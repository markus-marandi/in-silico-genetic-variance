#!/usr/bin/env python3
"""
build tss±pad bed from ensg list and mane summary.

args:
  ensg_list (str): path to ensg list
  mane (str): path to mane summary file
  pad_tss (int): padding around tss
  bed_out (str): path to write bed
  qc_out (str): path to write qc tsv
  missing_out (str): path to write missing ensg list

returns:
  None
"""

import argparse
import re
from pathlib import Path
from typing import Dict, Optional

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


def detect_columns(df: pd.DataFrame) -> Dict[str, str]:
    col_map = {
        "ensg": ["Ensembl_Gene", "Ensembl gene", "Ensembl gene ID", "EnsemblGene"],
        "symbol": ["Approved_Symbol", "HGNC_symbol", "symbol", "Gene", "Approved symbol", "HGNC Symbol"],
        "chrom": ["GRCh38_chr", "chr", "GRCh38 chromosome", "Chromosome"],
        "start": ["chr_start", "GRCh38_start", "start"],
        "end": ["chr_end", "GRCh38_end", "end"],
        "strand": ["chr_strand", "strand", "GRCh38_strand"],
    }
    resolved: Dict[str, Optional[str]] = {k: None for k in col_map}
    for key, candidates in col_map.items():
        for name in candidates:
            if name in df.columns:
                resolved[key] = name
                break
    missing = [k for k, v in resolved.items() if k in {"ensg", "chrom", "start", "end", "strand"} and v is None]
    if missing:
        raise ValueError(f"mane missing required columns: {missing}")
    return resolved


def load_mane(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    col_map = detect_columns(df)
    use_cols = [
        col_map["ensg"],
        col_map["chrom"],
        col_map["start"],
        col_map["end"],
        col_map["strand"],
    ]
    optional_symbol = col_map.get("symbol")
    if optional_symbol:
        use_cols.insert(1, optional_symbol)
    mane = df[use_cols].copy()
    if optional_symbol:
        mane.columns = [
            "Ensembl_Gene",
            "Symbol",
            "GRCh38_chr",
            "chr_start",
            "chr_end",
            "chr_strand",
        ]
    else:
        mane.columns = [
            "Ensembl_Gene",
            "GRCh38_chr",
            "chr_start",
            "chr_end",
            "chr_strand",
        ]
        mane["Symbol"] = pd.NA
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
    mane = (
        mane.sort_values(["Ensembl_Gene", "span"], ascending=[True, False])
        .drop_duplicates("Ensembl_Gene", keep="first")
    )
    return mane


def build_bed(mane: pd.DataFrame, ensg_list: list[str], pad_tss: int) -> pd.DataFrame:
    ensg_df = pd.DataFrame({"Ensembl_Gene": ensg_list})
    merged = ensg_df.merge(mane, on="Ensembl_Gene", how="left")
    mapped = merged.dropna(subset=["chrom"]).copy()
    if mapped.empty:
        raise ValueError("no ENSG IDs from the list were found in MANE")
    bed = mapped.apply(
        lambda r: pd.Series(
            {
                "chr": r["chrom"],
                "start": max(0, int(r["TSS"]) - pad_tss),
                "end": int(r["TSS"]) + pad_tss,
                "name": r["Ensembl_Gene"],
            }
        ),
        axis=1,
    )
    bed["start"] = pd.to_numeric(bed["start"], errors="coerce").fillna(0).astype(int)
    bed["end"] = pd.to_numeric(bed["end"], errors="coerce").fillna(0).astype(int)
    bed = bed[(bed["end"] > bed["start"]) & (bed["chr"].str.match(CANON_CHR_RE))]
    bed = bed.sort_values(["chr", "start", "end", "name"]).reset_index(drop=True)
    if bed.empty:
        raise ValueError("bed is empty after filtering")
    return bed


def write_bed(bed: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    bed.to_csv(path, sep="\t", header=False, index=False)


def write_qc(total_inputs_raw: int, total_inputs_unique: int, mapped: int, pad_tss: int, path: Path) -> None:
    missing = max(total_inputs_unique - mapped, 0)
    qc = pd.DataFrame(
        [
            {
                "total_inputs_raw": total_inputs_raw,
                "total_inputs_unique": total_inputs_unique,
                "mapped_unique": mapped,
                "missing_unique": missing,
                "pad_tss": pad_tss,
            }
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    qc.to_csv(path, sep="\t", index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="build TSS±pad BED from MANE")
    parser.add_argument("--ensg-list", required=True, help="path to ensg list")
    parser.add_argument("--mane", required=True, help="path to mane summary file")
    parser.add_argument("--pad-tss", required=True, type=int, help="padding around TSS in bp")
    parser.add_argument("--bed-out", required=True, help="path to write bed")
    parser.add_argument("--qc-out", required=True, help="path to write qc tsv")
    parser.add_argument("--missing-out", required=True, help="path to write missing ensg list")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensg_path = Path(args.ensg_list)
    if not ensg_path.exists():
        raise FileNotFoundError(f"ensg list not found: {ensg_path}")
    mane_path = Path(args.mane)
    if not mane_path.exists():
        raise FileNotFoundError(f"mane file not found: {mane_path}")

    ensgs_raw = [
        line.strip()
        for line in ensg_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if not ensgs_raw:
        raise ValueError("ensg list is empty")
    ensgs = list(dict.fromkeys(ensgs_raw))

    mane = load_mane(mane_path)
    bed = build_bed(mane, ensgs, args.pad_tss)
    write_bed(bed, Path(args.bed_out))
    write_qc(
        total_inputs_raw=len(ensgs_raw),
        total_inputs_unique=len(ensgs),
        mapped=bed["name"].nunique(),
        pad_tss=args.pad_tss,
        path=Path(args.qc_out),
    )
    missing = set(ensgs) - set(bed["name"].unique())
    missing_path = Path(args.missing_out)
    missing_path.parent.mkdir(parents=True, exist_ok=True)
    missing_path.write_text("\n".join(sorted(missing)), encoding="utf-8")


if __name__ == "__main__":
    main()


