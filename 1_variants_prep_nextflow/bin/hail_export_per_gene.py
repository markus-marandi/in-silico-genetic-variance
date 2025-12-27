#!/usr/bin/env python3
"""
export per-gene gnomad v4.1 variants for tss intervals.

args:
  bed (str): bed file with chr, start, end, name
  gnomad_ht (str): gnomad hail table path
  outdir (str): output directory for per-gene exports
  tmp_dir (str): hail temporary directory
  spark_conf (str): optional spark conf json path
  gcs_connector_jar (str): optional gcs connector jar path
  hail_home (str): optional hail home for placeholder substitution

returns:
  None
"""

import argparse
import os
import json
from pathlib import Path
from typing import Dict, Iterable, List

import hail as hl
from hail.utils import Interval
import pandas as pd


def load_bed(path: Path) -> pd.DataFrame:
    cols = ["chr", "start", "end", "name"]
    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=cols,
        dtype={"chr": str, "start": int, "end": int, "name": str},
    )
    if bed.empty:
        raise ValueError("bed file is empty")
    return bed


def bed_to_intervals(bed: pd.DataFrame) -> Dict[str, List[Interval]]:
    by_gene: Dict[str, List[Interval]] = {}
    for _, row in bed.iterrows():
        chrom = str(row["chr"])
        start = int(row["start"])
        end = int(row["end"])
        gene = str(row["name"])
        locus_interval = hl.parse_locus_interval(
            f"{chrom}:{start + 1}-{end}",
            reference_genome="GRCh38",
        )
        by_gene.setdefault(gene, []).append(locus_interval)
    return by_gene


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="export per-gene gnomAD v4.1 intervals with Hail")
    parser.add_argument("--bed", required=True, help="bed file with gene intervals")
    parser.add_argument("--gnomad-ht", required=True, help="gnomad hail table path")
    parser.add_argument("--outdir", required=True, help="directory for per-gene exports")
    parser.add_argument("--tmp-dir", required=True, help="hail tmp dir")
    parser.add_argument("--spark-conf", help="spark conf json path")
    parser.add_argument("--gcs-connector-jar", help="gcs connector jar path")
    parser.add_argument("--hail-home", help="hail home path for placeholder replacement")
    return parser.parse_args()


def export_gene_variants(ht: hl.Table, gene: str, intervals: Iterable[Interval], out_path: Path) -> None:
    filtered = hl.filter_intervals(ht, intervals)
    freq_struct = hl.or_missing(hl.len(filtered.freq) > 0, filtered.freq[0])
    export = filtered.select(
        CHROM=filtered.locus.contig,
        POS=filtered.locus.position,
        REF=filtered.alleles[0],
        ALT=filtered.alleles[1],
        AF=freq_struct.AF,
        AC=freq_struct.AC,
        AN=freq_struct.AN,
    )
    export = export.select_globals()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    hl.export_table(
        export,
        out_path.as_posix(),
        header=True,
        delimiter="\t",
        parallel="header_per_shard",
    )


def main() -> None:
    args = parse_args()
    bed_path = Path(args.bed)
    if not bed_path.exists():
        raise FileNotFoundError(f"bed not found: {bed_path}")
    placeholders = {
        "GCS_CONNECTOR_JAR": args.gcs_connector_jar or os.environ.get("GCS_CONNECTOR_JAR", ""),
        "HAIL_HOME": args.hail_home or os.environ.get("HAIL_HOME", ""),
        "TMP_DIR": args.tmp_dir,
    }
    spark_conf = None
    if args.spark_conf:
        spark_conf = load_spark_conf(Path(args.spark_conf), placeholders)
    conf_text_requires_jar = False
    if args.spark_conf:
        conf_text_requires_jar = "{GCS_CONNECTOR_JAR}" in Path(args.spark_conf).read_text()
    if conf_text_requires_jar and not placeholders["GCS_CONNECTOR_JAR"]:
        raise ValueError("spark_conf references {GCS_CONNECTOR_JAR} but no connector jar provided")

    bed = load_bed(bed_path)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log_path = outdir / "hail_export.log"
    hl.init(
        app_name="hail_export_per_gene",
        spark_conf=spark_conf,
        tmp_dir=args.tmp_dir,
        log=log_path.as_posix(),
    )
    ht = hl.read_table(args.gnomad_ht)
    ht = ht.filter(hl.len(ht.alleles) == 2)
    if "freq" not in ht.row:
        raise ValueError("gnomad table missing freq field required for AF/AC/AN")

    intervals = bed_to_intervals(bed)

    for gene, gene_intervals in intervals.items():
        out_path = outdir / f"{gene}.gnomad_v4.1.tsv.bgz"
        export_gene_variants(ht, gene, gene_intervals, out_path)


if __name__ == "__main__":
    main()


