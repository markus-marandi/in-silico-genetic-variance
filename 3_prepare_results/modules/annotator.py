from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import polars as pl


def annotate_af(lf: pl.LazyFrame, variants_path: Optional[Path]) -> Tuple[pl.LazyFrame, bool]:
    """join AF from initial variants file if available; return lf and flag indicating AF presence."""
    if variants_path is None or not variants_path.exists():
        return lf, False

    af_cols = ["CHROM", "POS", "REF", "ALT", "AF", "variant_id"]
    variants = pl.scan_csv(
        variants_path,
        separator="\t",
        has_header=True,
        null_values=["", "NA", "None"],
        infer_schema_length=1000,
    )
    keep = [c for c in variants.columns if c in af_cols]
    variants = variants.select(keep).with_columns(
        CHROM=pl.col("CHROM").cast(pl.Utf8),
        POS=pl.col("POS").cast(pl.Int64),
        REF=pl.col("REF").cast(pl.Utf8),
        ALT=pl.col("ALT").cast(pl.Utf8),
        variant_id=pl.col("variant_id").cast(pl.Utf8) if "variant_id" in variants.columns else pl.lit(None),
    )

    if "variant_id" in variants.columns and "variant_id" in lf.columns:
        joined = lf.join(variants, on="variant_id", how="left", suffix="_af")
    else:
        joined = lf.join(variants, on=["CHROM", "POS", "REF", "ALT"], how="left", suffix="_af")

    if "AF_af" in joined.columns:
        joined = joined.with_columns(
            pl.when(pl.col("AF_af").is_not_null()).then(pl.col("AF_af")).otherwise(pl.col("AF")).alias("AF")
        ).drop("AF_af", strict=False)

    any_af = joined.select(pl.col("AF").is_not_null().any()).collect().item()
    return joined, bool(any_af)


def annotate_gnomad(lf: pl.LazyFrame, gnomad_path: Optional[Path]) -> pl.LazyFrame:
    if gnomad_path is None:
        return lf
    if not gnomad_path.exists():
        raise FileNotFoundError(f"gnomad parquet not found: {gnomad_path}")
    gnomad = (
        pl.scan_parquet(gnomad_path)
        .select(
            pl.col("CHROM").cast(pl.Utf8),
            pl.col("POS").cast(pl.Int64),
            pl.col("REF").cast(pl.Utf8),
            pl.col("ALT").cast(pl.Utf8),
            pl.col("AF").alias("AF_gnomad"),
        )
    )
    return (
        lf.with_columns(
            CHROM=pl.col("CHROM").cast(pl.Utf8),
            POS=pl.col("POS").cast(pl.Int64),
            REF=pl.col("REF").cast(pl.Utf8),
            ALT=pl.col("ALT").cast(pl.Utf8),
        )
        .join(gnomad, on=["CHROM", "POS", "REF", "ALT"], how="left")
        .with_columns(
            pl.when(pl.col("AF").is_null()).then(pl.col("AF_gnomad")).otherwise(pl.col("AF")).alias("AF")
        )
        .drop("AF_gnomad", strict=False)
    )