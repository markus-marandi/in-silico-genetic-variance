#!/usr/bin/env python3
"""
compare tss choices between longest-transcript heuristic and mane select-first.

args:
    mane (str): path to mane summary tsv.
    genelist (str): optional path to restrict to these ENSG IDs.
    diff_out (str): optional path to write per-gene differences tsv.

returns:
    None
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Optional

import pandas as pd

CANON_CHR_RE = re.compile(r"^chr([1-9]|1[0-9]|2[0-2]|X|Y)$")


def strip_version(value: str) -> str:
    return re.sub(r"\.\d+$", "", value)


def nc_to_ucsc_chr(value: str) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    s = str(value)
    if s.startswith("chr") and CANON_CHR_RE.match(s):
        return s
    match = re.match(r"^NC_0+(\d+)\.\d+$", s)
    if match:
        num = int(match.group(1))
        if num == 23:
            return "chrX"
        if num == 24:
            return "chrY"
        chrom = f"chr{num}"
        return chrom if CANON_CHR_RE.match(chrom) else None
    if s.isdigit():
        chrom = f"chr{s}"
        return chrom if CANON_CHR_RE.match(chrom) else None
    if s in {"X", "Y"}:
        return f"chr{s}"
    return None


def detect_columns(df: pd.DataFrame) -> dict[str, str | None]:
    col_map = {
        "ensg": ["Ensembl_Gene", "Ensembl gene", "Ensembl gene ID", "EnsemblGene"],
        "symbol": ["Approved_Symbol", "HGNC_symbol", "symbol", "Gene", "Approved symbol", "HGNC Symbol"],
        "chrom": ["GRCh38_chr", "chr", "GRCh38 chromosome", "Chromosome"],
        "start": ["chr_start", "GRCh38_start", "start"],
        "end": ["chr_end", "GRCh38_end", "end"],
        "strand": ["chr_strand", "strand", "GRCh38_strand"],
        "status": ["MANE_status", "mane_status", "Status"],
    }
    resolved: dict[str, str | None] = {k: None for k in col_map}
    for key, candidates in col_map.items():
        for name in candidates:
            if name in df.columns:
                resolved[key] = name
                break
    missing = [k for k in ("ensg", "chrom", "start", "end", "strand") if resolved[k] is None]
    if missing:
        raise ValueError(f"mane missing required columns: {missing}")
    return resolved


def load_mane(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    cols = detect_columns(df)
    use_cols = [
        cols["ensg"],
        cols["chrom"],
        cols["start"],
        cols["end"],
        cols["strand"],
    ]
    optional_symbol = cols.get("symbol")
    optional_status = cols.get("status")
    if optional_symbol:
        use_cols.insert(1, optional_symbol)
    if optional_status:
        use_cols.append(optional_status)
    mane = df[use_cols].copy()
    if optional_symbol and optional_status:
        mane.columns = [
            "Ensembl_Gene",
            "Symbol",
            "GRCh38_chr",
            "chr_start",
            "chr_end",
            "chr_strand",
            "MANE_status",
        ]
    elif optional_symbol:
        mane.columns = [
            "Ensembl_Gene",
            "Symbol",
            "GRCh38_chr",
            "chr_start",
            "chr_end",
            "chr_strand",
        ]
        mane["MANE_status"] = pd.NA
    elif optional_status:
        mane.columns = [
            "Ensembl_Gene",
            "GRCh38_chr",
            "chr_start",
            "chr_end",
            "chr_strand",
            "MANE_status",
        ]
        mane["Symbol"] = pd.NA
    else:
        mane.columns = [
            "Ensembl_Gene",
            "GRCh38_chr",
            "chr_start",
            "chr_end",
            "chr_strand",
        ]
        mane["Symbol"] = pd.NA
        mane["MANE_status"] = pd.NA

    mane["Ensembl_Gene"] = mane["Ensembl_Gene"].map(strip_version)
    mane["chrom"] = mane["GRCh38_chr"].map(nc_to_ucsc_chr)
    mane["chr_start"] = pd.to_numeric(mane["chr_start"], errors="coerce")
    mane["chr_end"] = pd.to_numeric(mane["chr_end"], errors="coerce")
    mane = mane.dropna(subset=["chrom", "chr_start", "chr_end", "chr_strand"])
    mane = mane[mane["chrom"].str.match(CANON_CHR_RE)].copy()
    mane["TSS"] = mane.apply(
        lambda r: int(r["chr_start"]) if r["chr_strand"] == "+" else int(r["chr_end"]),
        axis=1,
    )
    mane["span"] = mane["chr_end"] - mane["chr_start"]
    return mane


def choose_longest(mane: pd.DataFrame) -> pd.DataFrame:
    return (
        mane.sort_values(["Ensembl_Gene", "span"], ascending=[True, False])
        .drop_duplicates("Ensembl_Gene", keep="first")
    )


def choose_mane_curated(mane: pd.DataFrame) -> pd.DataFrame:
    def priority(status: str | float | None) -> int:
        if isinstance(status, float) and pd.isna(status):
            return 2
        if status == "MANE Select":
            return 0
        if status == "MANE Plus Clinical":
            return 1
        return 2

    mane = mane.copy()
    mane["status_rank"] = mane["MANE_status"].map(priority)
    return (
        mane.sort_values(
            ["Ensembl_Gene", "status_rank", "span"], ascending=[True, True, False]
        )
        .drop_duplicates("Ensembl_Gene", keep="first")
        .drop(columns=["status_rank"])
    )


def compare_choices(mane: pd.DataFrame) -> pd.DataFrame:
    longest = choose_longest(mane)
    curated = choose_mane_curated(mane)
    merged = longest.merge(
        curated,
        on="Ensembl_Gene",
        how="outer",
        suffixes=("_longest", "_curated"),
        indicator=True,
    )
    changed = merged[
        (merged["chrom_longest"] != merged["chrom_curated"])
        | (merged["TSS_longest"] != merged["TSS_curated"])
        | (merged["MANE_status_longest"] != merged["MANE_status_curated"])
    ].copy()
    changed["tss_delta"] = merged["TSS_curated"] - merged["TSS_longest"]
    return changed


def summarize(changed: pd.DataFrame, total_genes: int) -> str:
    n_changed = len(changed)
    pct = (n_changed / total_genes * 100) if total_genes else 0.0
    median_delta = (
        changed["tss_delta"].abs().median() if not changed["tss_delta"].isna().all() else 0.0
    )
    parts = [
        f"total genes: {total_genes}",
        f"different picks: {n_changed} ({pct:.2f}%)",
        f"median |delta TSS|: {int(median_delta)} bp",
    ]
    return "\n".join(parts)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="compare tss choices between longest transcript and mane select-first"
    )
    parser.add_argument("--mane", required=True, help="path to mane summary tsv")
    parser.add_argument(
        "--genelist",
        help="optional path to ensg list (headerless; versions ok) to restrict comparison",
    )
    parser.add_argument(
        "--diff-out",
        help="optional path to write differences tsv (genes where choice changes)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mane_path = Path(args.mane)
    if not mane_path.exists():
        raise FileNotFoundError(f"mane file not found: {mane_path}")

    mane = load_mane(mane_path)
    if args.genelist:
        genes_path = Path(args.genelist)
        if not genes_path.exists():
            raise FileNotFoundError(f"genelist not found: {genes_path}")
        gene_ids = {
            strip_version(line.strip())
            for line in genes_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        }
        mane = mane[mane["Ensembl_Gene"].isin(gene_ids)].copy()

    longest = choose_longest(mane)
    curated = choose_mane_curated(mane)
    changed = compare_choices(mane)

    print(summarize(changed, total_genes=len(curated)))
    if args.diff_out:
        out_path = Path(args.diff_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        changed.to_csv(out_path, sep="\t", index=False)
        print(f"wrote differences to {out_path}")


if __name__ == "__main__":
    main()

