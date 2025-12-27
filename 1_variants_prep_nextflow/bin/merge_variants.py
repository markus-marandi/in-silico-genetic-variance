#!/usr/bin/env python3
"""
merge per-gene hail exports into one gzipped tsv.

args:
  per_gene_dir (str): directory containing per-gene exports
  output (str): path to merged gzipped tsv

returns:
  None
"""

import argparse
from pathlib import Path
from typing import List

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="merge per-gene gnomAD exports")
    parser.add_argument("--per-gene-dir", required=True, help="directory with per-gene tsv.bgz files")
    parser.add_argument("--output", required=True, help="path to merged output .tsv.gz")
    return parser.parse_args()


def read_gene_table(path: Path) -> pd.DataFrame:
    gene_tag = path.name.replace(".gnomad_v4.1.tsv.bgz", "")
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        return df
    df["gene_tag"] = gene_tag
    df["variant_id"] = (
        df["CHROM"].astype(str)
        + ":"
        + df["POS"].astype(str)
        + ":"
        + df["REF"].astype(str)
        + ">"
        + df["ALT"].astype(str)
    )
    return df


def merge_tables(files: List[Path]) -> pd.DataFrame:
    frames = [read_gene_table(path) for path in files]
    frames = [f for f in frames if not f.empty]
    if not frames:
        raise ValueError("no per-gene variant files with data found to merge")
    merged = pd.concat(frames, ignore_index=True)
    merged = merged.drop_duplicates(subset=["variant_id", "gene_tag"])
    return merged


def main() -> None:
    args = parse_args()
    per_gene_dir = Path(args.per_gene_dir)
    if not per_gene_dir.exists():
        raise FileNotFoundError(f"per_gene directory not found: {per_gene_dir}")
    files = sorted(per_gene_dir.glob("*.gnomad_v4.1.tsv.bgz"))
    if not files:
        raise FileNotFoundError("no per-gene files found to merge")
    merged = merge_tables(files)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(out_path, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main()


