#!/usr/bin/env python3
"""generate mutation-rate matched synthetic snvs for null variant sets.

generates all possible snvs within bed intervals, excludes observed variants,
and stratified samples to match the observed mutational spectrum.

args:
    bed (str): bed file with chr, start, end, gene_id
    observed_variants (str): tsv with observed variants (CHROM, POS, REF, ALT)
    output (str): output tsv path for synthetic variants
    tmp_dir (str): hail temporary directory
    spark_conf (str): optional spark conf json path
    mutation_rates (str): optional mutation rates table
    methylation_ht (str): methylation hail table path
    grch38_fasta (str): reference fasta path
    grch38_fai (str): reference fasta index path
    target_n (int): target number of synthetic variants
    hail_home (str): optional hail home for placeholder substitution
    gcs_connector_jar (str): optional gcs connector jar path

returns:
    None
"""

import argparse
import json
import os
import uuid
from pathlib import Path
from typing import Dict

import hail as hl
import pandas as pd


# methylation bin thresholds
METH_T1 = 0.3
METH_T2 = 0.6
REFERENCE = "GRCh38"
N_PARTS = 256


def _replace_placeholders(conf: Dict, placeholders: Dict[str, str]) -> Dict:
    updated = {}
    for key, value in conf.items():
        if isinstance(value, dict):
            updated[key] = _replace_placeholders(value, placeholders)
        elif isinstance(value, str):
            new_val = value
            for ph_key, ph_val in placeholders.items():
                new_val = new_val.replace(f"{{{ph_key}}}", ph_val)
            updated[key] = new_val
        else:
            updated[key] = value
    return updated


def load_spark_conf(path: Path, placeholders: Dict[str, str]) -> Dict:
    conf = json.loads(path.read_text())
    return _replace_placeholders(conf, placeholders)


def load_bed(path: Path) -> pd.DataFrame:
    """load bed file ensuring tab separation."""
    cols = ["chr", "start", "end", "gene_id"]
    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=cols,
        dtype={"chr": str, "start": int, "end": int, "gene_id": str},
    )
    if bed.empty:
        raise ValueError("bed file is empty")
    return bed


def bed_to_intervals(bed: pd.DataFrame):
    """convert bed to hail intervals list."""
    intervals = []
    for _, row in bed.iterrows():
        chrom = str(row["chr"])
        start = int(row["start"])
        end = int(row["end"])
        locus_interval = hl.parse_locus_interval(
            f"{chrom}:{start + 1}-{end}",
            reference_genome=REFERENCE,
        )
        intervals.append(locus_interval)
    return intervals


def bed_to_gene_mapping(bed: pd.DataFrame):
    """create interval to gene_id mapping table."""
    rows = []
    for _, row in bed.iterrows():
        chrom = str(row["chr"])
        start = int(row["start"])
        end = int(row["end"])
        gene_id = str(row["gene_id"])
        rows.append({
            "interval": hl.parse_locus_interval(
                f"{chrom}:{start + 1}-{end}",
                reference_genome=REFERENCE,
            ),
            "gene_id": gene_id,
        })
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            interval=hl.tinterval(hl.tlocus(REFERENCE)),
            gene_id=hl.tstr,
        ),
    ).key_by("interval")


def load_observed_variants(path: Path, min_parts: int = 256) -> hl.Table:
    """load observed variants tsv and key by locus/alleles."""
    ht = hl.import_table(
        str(path),
        impute=True,
        delimiter="\t",
        quote=None,
        min_partitions=min_parts,
    )

    ht = ht.annotate(
        locus=hl.locus(ht.CHROM, hl.int32(ht.POS), reference_genome=REFERENCE),
        alleles=[ht.REF, ht.ALT],
    )

    # keep only snvs
    ht = ht.filter((hl.len(ht.alleles[0]) == 1) & (hl.len(ht.alleles[1]) == 1))
    ht = ht.key_by("locus", "alleles")
    return ht


def bin_methylation(score):
    return (
        hl.case()
        .when(hl.is_missing(score), hl.int32(0))
        .when(score < METH_T1, hl.int32(0))
        .when(score < METH_T2, hl.int32(1))
        .default(hl.int32(2))
    )


def comp(b):
    return (
        hl.switch(b)
        .when("A", "T")
        .when("C", "G")
        .when("G", "C")
        .when("T", "A")
        .or_missing()
    )


def revcomp(trimer):
    return comp(trimer[2]) + comp(trimer[1]) + comp(trimer[0])


def canon_pyrimidine(ref, alt, ctx):
    """pyrimidine collapse: make ref in {C,T}."""
    flip = hl.set(["A", "G"]).contains(ref)
    c_ref = hl.if_else(flip, comp(ref), ref)
    c_alt = hl.if_else(flip, comp(alt), alt)
    c_ctx = hl.if_else(flip, revcomp(ctx), ctx)
    return c_ref, c_alt, c_ctx


def mut6_class(ref, alt, ctx):
    """compute 6-class mutation type."""
    c_ref, c_alt, _ = canon_pyrimidine(ref, alt, ctx)
    return hl.str(c_ref) + ">" + hl.str(c_alt)


def canon_AC(ref, alt, ctx):
    """a/c orient: make ref in {A,C}."""
    flip = hl.set(["G", "T"]).contains(ref)
    c_ref = hl.if_else(flip, comp(ref), ref)
    c_alt = hl.if_else(flip, comp(alt), alt)
    c_ctx = hl.if_else(flip, revcomp(ctx), ctx)
    return c_ref, c_alt, c_ctx


def load_mu_table(mu_path: str) -> hl.Table:
    """load karczewski mutation rate table."""
    mu = hl.import_table(
        mu_path,
        impute=True,
        delimiter="\t",
        types={"methylation_level": hl.tint32, "mu_snp": hl.tfloat64},
        force_bgz=True,
    )
    return mu.key_by("context", "ref", "alt", "methylation_level")


def annotate_mu_snp(ht: hl.Table, mu: hl.Table) -> hl.Table:
    """annotate mutation rate from karczewski table."""
    ref = ht.alleles[0]
    alt = ht.alleles[1]
    c_ref, c_alt, c_ctx = canon_AC(ref, alt, ht.context)

    # methylation only for CpG C>T
    c_meth = hl.if_else(
        (c_ref == "C") & (c_alt == "T") & (c_ctx[2] == "G"),
        ht.methyl_bin,
        hl.int32(0),
    )

    ht = ht.annotate(mu_snp=mu[c_ctx, c_ref, c_alt, c_meth].mu_snp)
    return ht


def _find_first_existing_field(ht: hl.Table, candidates):
    fields = set(ht.row_value.dtype.fields)
    for c in candidates:
        if c in fields:
            return c
    return None


def candidates_from_intervals(intervals, rg):
    """generate all possible snvs from reference within intervals."""
    ht = hl.Table.parallelize(
        [{"interval": iv} for iv in intervals],
        hl.tstruct(interval=hl.tinterval(hl.tlocus(REFERENCE))),
    )
    ht = ht.annotate(
        pos=hl.range(ht.interval.start.position, ht.interval.end.position)
    ).explode("pos")
    ht = ht.annotate(
        locus=hl.locus(ht.interval.start.contig, ht.pos, REFERENCE)
    ).key_by("locus")
    ht = ht.annotate(
        context=hl.get_sequence(
            ht.locus.contig, ht.locus.position, before=1, after=1, reference_genome=rg
        ),
        ref=hl.get_sequence(
            ht.locus.contig, ht.locus.position, before=0, after=0, reference_genome=rg
        ),
    )
    ht = ht.filter(
        hl.is_defined(ht.context)
        & (hl.len(ht.context) == 3)
        & hl.is_defined(ht.ref)
        & (hl.len(ht.ref) == 1)
        & (ht.context[1] == ht.ref)
    )
    ht = ht.annotate(
        alts=hl.array(["A", "C", "G", "T"]).filter(lambda b: b != ht.ref)
    ).explode("alts")
    ht = ht.annotate(alleles=[ht.ref, ht.alts]).key_by("locus", "alleles")
    return ht.select("context")


def compute_target_counts(observed: hl.Table, N: int):
    """compute per-stratum target counts from observed distribution."""
    dist = observed.aggregate(
        hl.agg.counter(hl.tuple([observed.mut6, observed.methyl_bin]))
    )
    tot = sum(dist.values())
    if tot == 0:
        raise ValueError("observed distribution is empty")

    target = {}
    running = 0
    items = list(dist.items())
    for i, (k, c) in enumerate(items):
        if i == len(items) - 1:
            n = N - running
        else:
            n = int(round(N * (c / tot)))
            running += n
        target[k] = max(0, n)

    # fix rounding drift
    drift = N - sum(target.values())
    if drift != 0:
        kmax = max(target.keys(), key=lambda kk: dist[kk])
        target[kmax] = max(0, target[kmax] + drift)
    return target


def cap_and_redistribute(target, available, N):
    """cap targets by available counts and redistribute deficit."""
    capped = {
        k: min(target.get(k, 0), available.get(k, 0))
        for k in set(target) | set(available)
    }
    deficit = N - sum(capped.values())
    if deficit <= 0:
        return capped

    for _ in range(50):
        if deficit <= 0:
            break
        rem = {k: available.get(k, 0) - capped.get(k, 0) for k in available}
        rem = {k: v for k, v in rem.items() if v > 0}
        if not rem:
            break
        kbest = max(rem.keys(), key=lambda kk: rem[kk])
        add = min(rem[kbest], deficit)
        capped[kbest] += add
        deficit -= add
    return capped


def stratified_sample_exact(cand: hl.Table, plan: dict):
    """perform exact stratified sampling per mut6/methyl_bin stratum."""
    parts = []
    for (mut6, mb), n in plan.items():
        if n <= 0:
            continue
        ht = cand.filter((cand.mut6 == mut6) & (cand.methyl_bin == mb))
        ht = ht.annotate(_r=hl.rand_unif()).order_by("_r").head(n).drop("_r")
        parts.append(ht)
    if not parts:
        return cand.head(0)
    out = parts[0]
    for ht in parts[1:]:
        out = out.union(ht)
    return out


def assign_gene_ids(ht: hl.Table, bed: pd.DataFrame) -> hl.Table:
    """assign gene_id to variants based on bed intervals."""
    # build interval table for gene assignment
    rows = []
    for _, row in bed.iterrows():
        chrom = str(row["chr"])
        start = int(row["start"])
        end = int(row["end"])
        gene_id = str(row["gene_id"])
        rows.append({"chrom": chrom, "start": start, "end": end, "gene_id": gene_id})
    
    gene_df = pd.DataFrame(rows)
    gene_ht = hl.Table.from_pandas(gene_df)
    gene_ht = gene_ht.annotate(
        interval=hl.locus_interval(
            gene_ht.chrom, gene_ht.start + 1, gene_ht.end, reference_genome=REFERENCE
        )
    ).key_by("interval")
    
    # annotate variants with gene_id
    ht = ht.annotate(
        _matches=hl.array([
            hl.struct(gene_id=g.gene_id)
            for g in hl.filter(
                lambda g: g.interval.contains(ht.locus),
                gene_ht.collect()
            )
        ])
    )
    
    # for simplicity, take first matching gene (intervals shouldn't overlap)
    ht = ht.annotate(gene_id=hl.if_else(
        hl.len(ht._matches) > 0,
        ht._matches[0].gene_id,
        "unknown"
    )).drop("_matches")
    
    return ht


def assign_gene_ids_simple(ht: hl.Table, bed: pd.DataFrame) -> hl.Table:
    """assign gene_id using pandas-based approach for efficiency."""
    # export to pandas, assign gene_id, reimport
    pdf = ht.to_pandas()
    
    # create interval lookup
    def find_gene(chrom, pos):
        for _, row in bed.iterrows():
            if row["chr"] == chrom and row["start"] <= pos < row["end"]:
                return row["gene_id"]
        return "unknown"
    
    # vectorized approach using interval tree would be better for large data
    # but for now use a simple merge approach
    pdf["_chrom"] = pdf["locus"].apply(lambda x: x.contig)
    pdf["_pos"] = pdf["locus"].apply(lambda x: x.position)
    
    # assign gene_id based on position
    gene_ids = []
    for _, row in pdf.iterrows():
        gene_ids.append(find_gene(row["_chrom"], row["_pos"]))
    pdf["gene_id"] = gene_ids
    pdf = pdf.drop(columns=["_chrom", "_pos"])
    
    return hl.Table.from_pandas(pdf).key_by("locus", "alleles")


def format_output_table(ht: hl.Table) -> hl.Table:
    """format output table with standard columns."""
    return ht.annotate(
        CHROM=ht.locus.contig,
        POS=ht.locus.position,
        REF=ht.alleles[0],
        ALT=ht.alleles[1],
        variant_id=hl.format(
            "%s:%d:%s>%s",
            ht.locus.contig,
            ht.locus.position,
            ht.alleles[0],
            ht.alleles[1],
        ),
    ).select("CHROM", "POS", "REF", "ALT", "gene_id", "variant_id")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="generate synthetic snvs")
    parser.add_argument("--bed", required=True, help="bed file with gene intervals")
    parser.add_argument("--observed-variants", required=True, help="observed variants tsv")
    parser.add_argument("--output", required=True, help="output tsv path")
    parser.add_argument("--tmp-dir", required=True, help="hail tmp dir")
    parser.add_argument("--spark-conf", help="spark conf json path")
    parser.add_argument("--mutation-rates", help="mutation rates table path")
    parser.add_argument("--methylation-ht", required=True, help="methylation hail table")
    parser.add_argument("--grch38-fasta", required=True, help="reference fasta")
    parser.add_argument("--grch38-fai", required=True, help="reference fasta index")
    parser.add_argument("--target-n", type=int, default=1_800_000, help="target synthetic count")
    parser.add_argument("--hail-home", help="hail home path")
    parser.add_argument("--gcs-connector-jar", help="gcs connector jar path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    bed_path = Path(args.bed)
    if not bed_path.exists():
        raise FileNotFoundError(f"bed not found: {bed_path}")

    observed_path = Path(args.observed_variants)
    if not observed_path.exists():
        raise FileNotFoundError(f"observed variants not found: {observed_path}")

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # setup spark conf
    placeholders = {
        "GCS_CONNECTOR_JAR": args.gcs_connector_jar or os.environ.get("GCS_CONNECTOR_JAR", ""),
        "HAIL_HOME": args.hail_home or os.environ.get("HAIL_HOME", ""),
        "TMP_DIR": args.tmp_dir,
    }
    spark_conf = None
    if args.spark_conf:
        spark_conf = load_spark_conf(Path(args.spark_conf), placeholders)
        spark_conf.setdefault("spark.sql.shuffle.partitions", str(N_PARTS))

    # init hail
    log_path = out_path.parent / "synthetic_generation.log"
    hl.init(
        app_name="generate_synthetic_snvs",
        spark_conf=spark_conf,
        tmp_dir=args.tmp_dir,
        log=str(log_path),
        default_reference=REFERENCE,
    )

    # add reference sequence
    rg = hl.get_reference(REFERENCE)
    if not rg.has_sequence():
        rg.add_sequence(args.grch38_fasta, args.grch38_fai)

    # load bed and observed variants
    print("loading bed file...")
    bed = load_bed(bed_path)
    intervals = bed_to_intervals(bed)

    print("loading observed variants...")
    observed = load_observed_variants(observed_path, min_parts=N_PARTS)

    # load methylation table
    print("loading methylation table...")
    meth = hl.read_table(args.methylation_ht)
    meth_field = _find_first_existing_field(
        meth,
        ["methylation", "methylation_level", "mean_methylation", "meth", "beta", "methylation_score"],
    )
    if meth_field is None:
        for f in meth.row_value.dtype.fields:
            if meth.row_value.dtype[f] in (hl.tfloat32, hl.tfloat64):
                meth_field = f
                break
    if meth_field is None:
        raise ValueError("could not detect methylation score field")

    # annotate observed with context and mutation class
    print("annotating observed variants...")
    observed = observed.annotate(
        context=hl.get_sequence(
            observed.locus.contig,
            observed.locus.position,
            before=1,
            after=1,
            reference_genome=rg,
        ),
        meth_score=meth[observed.locus][meth_field],
    )
    observed = observed.filter(
        hl.is_defined(observed.context)
        & (hl.len(observed.context) == 3)
        & (observed.context[1] == observed.alleles[0])
    )
    observed = observed.annotate(
        methyl_bin=bin_methylation(observed.meth_score),
        mut6=mut6_class(observed.alleles[0], observed.alleles[1], observed.context),
    )

    obs_ct = observed.count()
    print(f"observed snvs: {obs_ct}")
    if obs_ct == 0:
        raise RuntimeError("observed set empty")

    # generate candidates from reference
    print("generating candidate snvs from reference...")
    cand = candidates_from_intervals(intervals, rg)
    cand = cand.annotate(meth_score=meth[cand.locus][meth_field])
    cand = cand.annotate(
        methyl_bin=bin_methylation(cand.meth_score),
        mut6=mut6_class(cand.alleles[0], cand.alleles[1], cand.context),
    ).select("context", "methyl_bin", "mut6")

    cand0 = cand.count()
    print(f"candidate snvs: {cand0}")

    # exclude observed variants
    print("excluding observed variants...")
    excluded = observed.select().key_by("locus", "alleles")
    cand = cand.anti_join(excluded).persist()
    cand1 = cand.count()
    print(f"candidates after exclusion: {cand1}")

    if cand1 == 0:
        raise RuntimeError("all candidates excluded")

    # optionally annotate mutation rates
    mu = None
    if args.mutation_rates:
        print("loading mutation rates...")
        mu = load_mu_table(args.mutation_rates)
        cand = annotate_mu_snp(cand, mu)
        cand = cand.filter(hl.is_defined(cand.mu_snp)).persist()
        cand2 = cand.count()
        print(f"candidates with mu_snp: {cand2}")
        if cand2 == 0:
            raise RuntimeError("mu_snp missing for all candidates")
        cand1 = cand2

    # compute sampling plan
    print("computing sampling plan...")
    N = min(args.target_n, cand1)
    obs_dist = observed.select("mut6", "methyl_bin")
    avail = cand.aggregate(hl.agg.counter(hl.tuple([cand.mut6, cand.methyl_bin])))
    plan = compute_target_counts(obs_dist, N)
    plan = cap_and_redistribute(plan, avail, N)
    planned = sum(plan.values())
    print(f"planned sample: {planned}, target: {N}")

    # stratified sampling
    print("stratified sampling...")
    cand = cand.repartition(N_PARTS, shuffle=True).persist()
    sampled = stratified_sample_exact(cand, plan).persist()
    samp_ct = sampled.count()
    print(f"sampled: {samp_ct}")

    # assign gene_ids
    print("assigning gene ids...")
    # export sampled, assign gene_ids in pandas, then format
    sampled_pdf = sampled.to_pandas()
    sampled_pdf["CHROM"] = sampled_pdf["locus"].apply(lambda x: x.contig)
    sampled_pdf["POS"] = sampled_pdf["locus"].apply(lambda x: x.position)
    sampled_pdf["REF"] = sampled_pdf["alleles"].apply(lambda x: x[0])
    sampled_pdf["ALT"] = sampled_pdf["alleles"].apply(lambda x: x[1])
    
    # assign gene_id from bed
    def find_gene(chrom, pos):
        for _, row in bed.iterrows():
            if row["chr"] == chrom and row["start"] <= pos < row["end"]:
                return row["gene_id"]
        return "unknown"
    
    sampled_pdf["gene_id"] = sampled_pdf.apply(
        lambda r: find_gene(r["CHROM"], r["POS"]), axis=1
    )
    sampled_pdf["variant_id"] = sampled_pdf.apply(
        lambda r: f"{r['CHROM']}:{r['POS']}:{r['REF']}>{r['ALT']}", axis=1
    )

    # output columns
    out_cols = ["CHROM", "POS", "REF", "ALT", "gene_id", "variant_id"]
    out_pdf = sampled_pdf[out_cols]

    # write output
    print(f"writing output to {out_path}...")
    out_pdf.to_csv(out_path, sep="\t", index=False, compression="gzip")

    print(f"done. wrote {len(out_pdf)} synthetic variants")


if __name__ == "__main__":
    main()
