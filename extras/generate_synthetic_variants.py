"""unified workflow for generating mutation-rate matched synthetic variants.

this script:
1. automatically loads and keys raw variant tsvs if needed
2. generates mutation-rate matched synthetic variants
3. enforces cross-dataset disjointness
4. produces output tables and distribution plots

usage:
    python generate_synthetic_variants.py

configuration via environment variables:
    INTERMEDIATE_DIR: base directory for raw tsvs and keyed tables (default: /mnt/sdb/markus_files/gene_exp/intermediate)
    OUTPUT_DIR: output directory for final results (default: /mnt/sdb/markus_files/gene_exp/51_mutagenesis/ISM_variants)
    SPARK_CONF_JSON: path to spark configuration json
    HAIL_HOME: hail installation directory
    HAIL_TMP: temporary directory for hail
    HAIL_LOG: log file path

expected directory structure under INTERMEDIATE_DIR:
    {dataset}_variants.tsv          - raw variant file (will be auto-keyed)
    {dataset}_variants.keyed.ht/    - keyed hail table (created if missing)
    {dataset}_gene_set±10kb.tab.bed - target interval bed file

where {dataset} is ClinGen or Background
"""

import os
import json
import uuid
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import hail as hl

sys.path.insert(0, str(Path(__file__).parent.parent))
from helpers.path_manager import SyntheticVariantPaths


# config from env vars with defaults
BASE_DIR = Path(os.getenv("INTERMEDIATE_DIR", "/mnt/sdb/markus_files/gene_exp/intermediate"))
OUTPUT_DIR = Path(os.getenv("OUTPUT_DIR", "/mnt/sdb/markus_files/gene_exp/51_mutagenesis/ISM_variants"))
SPARK_CONF_JSON = os.getenv("SPARK_CONF_JSON", "/mnt/sdb/markus_files/gene_exp/11_prep_variants/spark_conf.json")
HAIL_HOME = os.getenv("HAIL_HOME", "/mnt/sdb/markus_files/hail")
TMP_DIR = Path(os.getenv("HAIL_TMP", "/mnt/sdb/markus_files/hail_tmp"))
LOG = Path(os.getenv("HAIL_LOG", "/mnt/sdb/markus_files/hail_logs/synthetic_variants.log"))

MUTATION_RATES = "/mnt/sdb/markus_files/gene_exp/51_mutagenesis/supplementary_dataset_10_mutation_rates.tsv.bgz"
METHYLATION_HT = "gs://gcp-public-data--gnomad/resources/grch38/methylation_sites/methylation.ht"

USE_GNOMAD_SITES = False
GNOMAD_SITES = "gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht/"

REFERENCE = "GRCh38"
TARGET_N_PER_SET = 1_800_000

GRCH38_FASTA = "gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz"
GRCH38_FAI = "gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai"

# methylation bin thresholds
METH_T1 = 0.3
METH_T2 = 0.6

N_PARTS = 256

# dataset names to process
DATASET_NAMES = ["ClinGen", "Background"]

os.makedirs(LOG.parent, exist_ok=True)
os.makedirs(TMP_DIR, exist_ok=True)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# -----------------------------
# Hail init
# -----------------------------
def init_hail():
    """initialize hail with configured spark settings."""
    with open(SPARK_CONF_JSON) as f:
        spark_conf = json.load(f)
    for k, v in list(spark_conf.items()):
        if isinstance(v, str):
            spark_conf[k] = v.replace("{hail_home}", HAIL_HOME)

    spark_conf.setdefault("spark.sql.shuffle.partitions", str(N_PARTS))

    hl.init(
        spark_conf=spark_conf,
        tmp_dir=str(TMP_DIR),
        log=str(LOG),
        default_reference=REFERENCE,
    )

    rg = hl.get_reference("GRCh38")
    if not rg.has_sequence():
        rg.add_sequence(GRCH38_FASTA, GRCH38_FAI)
    return rg


def load_and_key(tsv_in: str, ht_out: str, min_parts: int = 256) -> None:
    """import raw variant tsv and create keyed hail table.
    
    args:
        tsv_in (str): path to input tsv with CHROM, POS, REF, ALT columns.
        ht_out (str): output path for keyed hail table.
        min_parts (int): minimum partitions for import and repartitioning.
    """
    print(f"loading and keying: {tsv_in} -> {ht_out}")
    ht = hl.import_table(
        tsv_in, 
        impute=True, 
        delimiter="\t", 
        quote=None, 
        min_partitions=min_parts
    )

    ht = ht.annotate(
        locus=hl.locus(ht.CHROM, hl.int32(ht.POS), reference_genome=REFERENCE),
        alleles=[ht.REF, ht.ALT],
    )

    # keep only snvs
    ht = ht.filter((hl.len(ht.alleles[0]) == 1) & (hl.len(ht.alleles[1]) == 1))

    # key then repartition to parallelize sort
    ht = ht.key_by("locus", "alleles")
    ht = ht.repartition(min_parts, shuffle=True)

    ht = ht.checkpoint(ht_out, overwrite=True)
    n = ht.count()
    print(f"wrote {ht_out}, n={n}")


def ensure_keyed_table(paths: SyntheticVariantPaths) -> str:
    """ensure keyed hail table exists; create if missing.
    
    args:
        paths (SyntheticVariantPaths): path configuration for dataset.
    
    returns:
        str: path to keyed hail table.
    """
    if paths.keyed_ht.exists():
        print(f"keyed table exists: {paths.keyed_ht}")
        return str(paths.keyed_ht)
    
    if not paths.raw_variants_tsv.exists():
        raise FileNotFoundError(f"raw variants tsv not found: {paths.raw_variants_tsv}")
    
    load_and_key(str(paths.raw_variants_tsv), str(paths.keyed_ht), min_parts=N_PARTS)
    return str(paths.keyed_ht)


def _find_first_existing_field(ht: hl.Table, candidates):
    fields = set(ht.row_value.dtype.fields)
    for c in candidates:
        if c in fields:
            return c
    return None

def import_bed_safe(bed_path: str, reference_genome: str = "GRCh38") -> hl.Table:
    tmp_bed = os.path.join(TMP_DIR, f"bed_noheader_{uuid.uuid4().hex}.bed")
    with open(bed_path, "rt") as fin:
        first = fin.readline()
        rest = fin.read()

    parts = first.rstrip("\n").split("\t")
    has_header = (len(parts) < 3) or (not parts[1].isdigit()) or (not parts[2].isdigit())

    with open(tmp_bed, "wt") as fout:
        if not has_header:
            fout.write(first)
        fout.write(rest)

    contig_recoding = {
        **{str(i): f"chr{i}" for i in range(1, 23)},
        "X": "chrX",
        "Y": "chrY",
        "MT": "chrM",
    }
    return hl.import_bed(tmp_bed, reference_genome=reference_genome, contig_recoding=contig_recoding)

def load_target_intervals(bed_path: str, reference_genome: str = "GRCh38"):
    bed = import_bed_safe(bed_path, reference_genome=reference_genome)
    return bed.aggregate(hl.agg.collect(bed.interval))

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

# pyrimidine-collapse (for 6-class spectrum): make ref in {C,T} by flipping if ref in {A,G}
def canon_pyrimidine(ref, alt, ctx):
    flip = hl.set(["A","G"]).contains(ref)
    c_ref = hl.if_else(flip, comp(ref), ref)
    c_alt = hl.if_else(flip, comp(alt), alt)
    c_ctx = hl.if_else(flip, revcomp(ctx), ctx)
    return c_ref, c_alt, c_ctx

def mut6_class(ref, alt, ctx):
    c_ref, c_alt, _ = canon_pyrimidine(ref, alt, ctx)
    return hl.str(c_ref) + ">" + hl.str(c_alt)   # one of C>A,C>G,C>T,T>A,T>C,T>G

# A/C-orient (for Karczewski mu table): make ref in {A,C} by flipping if ref in {G,T}
def canon_AC(ref, alt, ctx):
    flip = hl.set(["G","T"]).contains(ref)
    c_ref = hl.if_else(flip, comp(ref), ref)
    c_alt = hl.if_else(flip, comp(alt), alt)
    c_ctx = hl.if_else(flip, revcomp(ctx), ctx)
    return c_ref, c_alt, c_ctx

def load_mu_table(mu_path: str) -> hl.Table:
    mu = hl.import_table(
        mu_path,
        impute=True,
        delimiter="\t",
        types={"methylation_level": hl.tint32, "mu_snp": hl.tfloat64},
        force_bgz=True,
    )
    return mu.key_by("context", "ref", "alt", "methylation_level")

def annotate_mu_snp(ht: hl.Table, mu: hl.Table) -> hl.Table:
    ref = ht.alleles[0]
    alt = ht.alleles[1]
    c_ref, c_alt, c_ctx = canon_AC(ref, alt, ht.context)

    # methylation only used for CpG C>T in Karczewski table, else 0
    c_meth = hl.if_else(
        (c_ref == "C") & (c_alt == "T") & (c_ctx[2] == "G"),
        ht.methyl_bin,
        hl.int32(0),
    )

    ht = ht.annotate(
        mu_snp=mu[c_ctx, c_ref, c_alt, c_meth].mu_snp
    )
    return ht

def compute_target_counts(observed: hl.Table, N: int):
    # observed has fields mut6, methyl_bin
    dist = observed.aggregate(hl.agg.counter(hl.tuple([observed.mut6, observed.methyl_bin])))
    tot = sum(dist.values())
    if tot == 0:
        raise ValueError("Observed distribution is empty.")
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
    capped = {k: min(target.get(k, 0), available.get(k, 0)) for k in set(target) | set(available)}
    deficit = N - sum(capped.values())
    if deficit <= 0:
        return capped
    # greedily fill deficit into buckets with remaining capacity
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

def format_output_table(ht: hl.Table) -> hl.Table:
    return ht.annotate(
        CHROM=ht.locus.contig,
        POS=ht.locus.position,
        REF=ht.alleles[0],
        ALT=ht.alleles[1],
        gene_tag=".",
        variant_id=hl.format("%s:%d:%s>%s", ht.locus.contig, ht.locus.position, ht.alleles[0], ht.alleles[1]),
    ).select("CHROM", "POS", "REF", "ALT", "gene_tag", "variant_id")

def _plot_distribution(counts, title, path):
    labels = [f"{k[0]}|{k[1]}" for k in counts]
    values = [counts[k] for k in counts]
    plt.figure(figsize=(8, 4))
    plt.bar(range(len(labels)), values)
    plt.xticks(range(len(labels)), labels, rotation=45, ha="right")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(path)
    plt.close()

def candidates_from_bed_fasta(intervals, reference_genome: str = "GRCh38"):
    ht = hl.Table.parallelize(
        [{"interval": iv} for iv in intervals],
        hl.tstruct(interval=hl.tinterval(hl.tlocus(reference_genome))),
    )
    ht = ht.annotate(pos=hl.range(ht.interval.start.position, ht.interval.end.position)).explode("pos")
    ht = ht.annotate(locus=hl.locus(ht.interval.start.contig, ht.pos, reference_genome)).key_by("locus")
    ht = ht.annotate(
        context=hl.get_sequence(ht.locus.contig, ht.locus.position, before=1, after=1, reference_genome=rg),
        ref=hl.get_sequence(ht.locus.contig, ht.locus.position, before=0, after=0, reference_genome=rg),
    )
    ht = ht.filter(
        hl.is_defined(ht.context)
        & (hl.len(ht.context) == 3)
        & hl.is_defined(ht.ref)
        & (hl.len(ht.ref) == 1)
        & (ht.context[1] == ht.ref)
    )
    ht = ht.annotate(alts=hl.array(["A", "C", "G", "T"]).filter(lambda b: b != ht.ref)).explode("alts")
    ht = ht.annotate(alleles=[ht.ref, ht.alts]).key_by("locus", "alleles")
    return ht.select("context")

# -----------------------------
# Main workflow
# -----------------------------
def main():
    """run complete synthetic variant generation workflow."""
    
    # init hail
    rg = init_hail()
    
    # build path configs for each dataset
    path_configs = {
        name: SyntheticVariantPaths(
            dataset_name=name,
            base_dir=BASE_DIR,
            output_dir=OUTPUT_DIR,
        )
        for name in DATASET_NAMES
    }
    
    # ensure all keyed tables exist (run load_and_key if needed)
    print("\n=== checking/creating keyed tables ===")
    keyed_tables = {}
    for name, paths in path_configs.items():
        keyed_tables[name] = ensure_keyed_table(paths)
    
    # load shared tables once
    print("\n=== loading shared reference tables ===")
    mu = load_mu_table(MUTATION_RATES)
    
    meth = hl.read_table(METHYLATION_HT)
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
        raise ValueError("could not detect methylation score field in methylation ht")
    
    gnomad_sites = None
    if USE_GNOMAD_SITES:
        gnomad_sites = hl.read_table(GNOMAD_SITES)
        if list(gnomad_sites.key) != ["locus", "alleles"]:
            gnomad_sites = gnomad_sites.key_by("locus", "alleles")
        gnomad_sites = gnomad_sites.select()
    
    # accumulates sampled keys to enforce disjointness across datasets
    global_excluded = None
    
    # process each dataset
    for name in DATASET_NAMES:
        print(f"\n=== processing {name} ===")
        
        paths = path_configs[name]
        
        intervals = load_target_intervals(str(paths.bed_file), reference_genome=REFERENCE)
        
        # (A) observed scored variants restricted to targets
        scored = hl.read_table(keyed_tables[name])
        scored_locus = scored.key_by("locus")
        scored_in_targets = hl.filter_intervals(scored_locus, intervals).key_by("locus", "alleles")

        scored_in_targets = scored_in_targets.annotate(
            context=hl.get_sequence(
                scored_in_targets.locus.contig,
                scored_in_targets.locus.position,
                before=1,
                after=1,
                reference_genome=rg,
            ),
            meth_score=meth[scored_in_targets.locus][meth_field]
        )
        scored_in_targets = scored_in_targets.filter(
            hl.is_defined(scored_in_targets.context) & (hl.len(scored_in_targets.context) == 3)
        )
        print(
            "context missing:",
            scored_in_targets.aggregate(hl.agg.count_where(hl.is_missing(scored_in_targets.context))),
        )
        print(
            "bad center base:",
            scored_in_targets.aggregate(
                hl.agg.count_where(
                    hl.is_defined(scored_in_targets.context)
                    & (hl.len(scored_in_targets.context) == 3)
                    & (scored_in_targets.context[1] != scored_in_targets.alleles[0])
                )
            ),
        )
        scored_in_targets = scored_in_targets.filter(
            scored_in_targets.context[1] == scored_in_targets.alleles[0]
        )
        scored_in_targets = scored_in_targets.annotate(
            methyl_bin=bin_methylation(scored_in_targets.meth_score),
            mut6=mut6_class(
                scored_in_targets.alleles[0], 
                scored_in_targets.alleles[1], 
                scored_in_targets.context
            ),
        )
        
        obs_dist = scored_in_targets.select("mut6", "methyl_bin")
        obs_ct = obs_dist.count()
        print("observed snvs in targets:", obs_ct)
        if obs_ct == 0:
            raise RuntimeError(f"{name}: observed set empty inside targets")
        obs_counts = obs_dist.aggregate(hl.agg.counter(hl.tuple([obs_dist.mut6, obs_dist.methyl_bin])))
        
        # (B) candidates from fasta within targets
        cand = candidates_from_bed_fasta(intervals, reference_genome=REFERENCE)
        cand = cand.annotate(meth_score=meth[cand.locus][meth_field])
        cand = cand.annotate(
            methyl_bin=bin_methylation(cand.meth_score),
            mut6=mut6_class(cand.alleles[0], cand.alleles[1], cand.context),
        ).select("context", "methyl_bin", "mut6")
        
        cand0 = cand.count()
        print("candidate snvs in targets:", cand0)
        
        # (C) exclude observed + cross-set + optional gnomad
        excluded = scored_in_targets.select().key_by("locus", "alleles")
        if USE_GNOMAD_SITES:
            gn_in = hl.filter_intervals(
                gnomad_sites.key_by("locus"), 
                intervals
            ).key_by("locus", "alleles").select()
            excluded = excluded.union(gn_in)
        if global_excluded is not None:
            excluded = excluded.union(global_excluded)
        
        excluded = excluded.distinct().persist()
        excl_ct = excluded.count()
        print("excluded keys:", excl_ct)
        
        cand = cand.anti_join(excluded).persist()
        cand1 = cand.count()
        print("candidates after exclusion:", cand1)
        if cand1 == 0:
            raise RuntimeError(f"{name}: all candidates excluded")
        
        # (D) annotate mu_snp and keep only defined
        cand = annotate_mu_snp(cand, mu)
        cand = cand.filter(hl.is_defined(cand.mu_snp)).persist()
        cand2 = cand.count()
        print("candidates with mu_snp defined:", cand2)
        if cand2 == 0:
            raise RuntimeError(f"{name}: mu_snp missing for all candidates")
        
        # (E) plan sample size from observed distribution
        N = min(TARGET_N_PER_SET, cand2)
        avail = cand.aggregate(hl.agg.counter(hl.tuple([cand.mut6, cand.methyl_bin])))
        plan = compute_target_counts(obs_dist, N)
        plan = cap_and_redistribute(plan, avail, N)
        planned = sum(plan.values())
        print("planned sample:", planned, "target:", N)
        
        # (F) stratified sampling
        cand = cand.repartition(N_PARTS, shuffle=True).persist()
        sampled = stratified_sample_exact(cand, plan).persist()
        samp_ct = sampled.count()
        print("sampled:", samp_ct)
        samp_counts = sampled.aggregate(hl.agg.counter(hl.tuple([sampled.mut6, sampled.methyl_bin])))
        
        # write outputs
        cand.write(str(paths.candidates_ht), overwrite=True)
        sampled.write(str(paths.sampled_ht), overwrite=True)
        tidy_sampled = format_output_table(sampled)
        tidy_sampled.export(str(paths.sampled_tsv))
        
        stats = sampled.aggregate(hl.struct(
            n=hl.agg.count(),
            mu_min=hl.agg.min(sampled.mu_snp),
            mu_max=hl.agg.max(sampled.mu_snp),
        ))
        print("final stats:", stats)
        print(sampled.key)
        sampled.describe()
        
        _plot_distribution(obs_counts, f"{name} observed mut6|methyl_bin", str(paths.observed_dist_plot))
        _plot_distribution(samp_counts, f"{name} sampled mut6|methyl_bin", str(paths.sampled_dist_plot))
        
        # update global exclusion for cross-set disjointness
        sampled_keys = sampled.key_by("locus", "alleles").select()
        global_excluded = (
            sampled_keys if global_excluded is None 
            else global_excluded.union(sampled_keys).distinct()
        )
    
    print("\ndone")


if __name__ == "__main__":
    main()
