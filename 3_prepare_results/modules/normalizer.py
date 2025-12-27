from __future__ import annotations

from pathlib import Path

import polars as pl

# Import functions directly
from .normalisation_helper import friendly_method_name, ag_variant_to_canonical


VARIANT_RE = r"^(?:chr)?(?P<chrom>[0-9XYM]+):(?P<pos>\d+):(?P<ref>[ACGTN]+)>(?P<alt>[ACGTN]+)$"


def normalize_and_backfill(lf: pl.LazyFrame) -> pl.LazyFrame:
    """ensure core columns and identifiers are present for downstream steps."""
    friendly_fn = friendly_method_name
    canon_fn = ag_variant_to_canonical

    # 1. Collect Schema ONCE (Schema-Safe)
    # This prevents the "ColumnNotFoundError" by checking what actually exists first.
    schema_cols = set(lf.collect_schema().names())

    has_vid = "variant_id" in schema_cols
    has_chrom = "CHROM" in schema_cols
    has_pos = "POS" in schema_cols
    has_ref = "REF" in schema_cols
    has_alt = "ALT" in schema_cols
    has_gene_tag = "gene_tag" in schema_cols
    has_gene_id = "gene_id" in schema_cols
    has_method = "method_friendly" in schema_cols
    has_scorer = "scorer_friendly" in schema_cols
    has_vid_canon = "variant_id_canonical" in schema_cols
    has_var_scorer = "variant_scorer" in schema_cols

    # 2. Ensure variant_id exists
    lf = lf.with_columns(
        pl.col("variant_id").cast(pl.Utf8) if has_vid else pl.lit(None).alias("variant_id")
    )

    # 3. Extract parts from variant_id (for backfilling)
    lf = lf.with_columns(
        pl.when(pl.col("variant_id").is_not_null())
        .then(pl.col("variant_id").str.extract(VARIANT_RE, group_index=1))
        .otherwise(None)
        .alias("_chrom_from_id"),
        pl.when(pl.col("variant_id").is_not_null())
        .then(pl.col("variant_id").str.extract(VARIANT_RE, group_index=2).cast(pl.Int64))
        .otherwise(None)
        .alias("_pos_from_id"),
        pl.when(pl.col("variant_id").is_not_null())
        .then(pl.col("variant_id").str.extract(VARIANT_RE, group_index=3))
        .otherwise(None)
        .alias("_ref_from_id"),
        pl.when(pl.col("variant_id").is_not_null())
        .then(pl.col("variant_id").str.extract(VARIANT_RE, group_index=4))
        .otherwise(None)
        .alias("_alt_from_id"),
    )

    # 4. Create or cast core columns safely
    lf = lf.with_columns(
        pl.col("CHROM").cast(pl.Utf8) if has_chrom else pl.lit(None).alias("CHROM"),
        pl.col("POS").cast(pl.Int64) if has_pos else pl.lit(None).alias("POS"),
        pl.col("REF").cast(pl.Utf8) if has_ref else pl.lit(None).alias("REF"),
        pl.col("ALT").cast(pl.Utf8) if has_alt else pl.lit(None).alias("ALT"),
    ).with_columns(
        CHROM=pl.col("CHROM").fill_null(pl.col("_chrom_from_id")),
        POS=pl.col("POS").fill_null(pl.col("_pos_from_id")),
        REF=pl.col("REF").fill_null(pl.col("_ref_from_id")),
        ALT=pl.col("ALT").fill_null(pl.col("_alt_from_id")),
    )

    # 5. Consolidate Gene Tag (SAFE LOGIC)
    # Only references columns that definitely exist
    if has_gene_tag and has_gene_id:
        tag_expr = pl.when(pl.col("gene_tag").is_null()).then(pl.col("gene_id")).otherwise(pl.col("gene_tag"))
    elif has_gene_tag:
        tag_expr = pl.col("gene_tag")
    elif has_gene_id:
        tag_expr = pl.col("gene_id")
    else:
        tag_expr = pl.lit(None)

    lf = lf.with_columns(gene_tag=tag_expr.alias("gene_tag"))

    # 6. Friendly Method Name
    if has_method:
        # If method_friendly exists, prioritize it, backfill from scorer
        fallback_col = "scorer_friendly" if has_scorer else ("variant_scorer" if has_var_scorer else None)
        if fallback_col:
             method_expr = (
                pl.when(pl.col("method_friendly").is_null())
                .then(pl.col(fallback_col).map_elements(friendly_fn, return_dtype=pl.Utf8))
                .otherwise(pl.col("method_friendly"))
             )
        else:
             method_expr = pl.col("method_friendly")
    else:
        # Create it from scorer
        src_col = "scorer_friendly" if has_scorer else ("variant_scorer" if has_var_scorer else None)
        if src_col:
            method_expr = pl.col(src_col).map_elements(friendly_fn, return_dtype=pl.Utf8)
        else:
            method_expr = pl.lit(None).cast(pl.Utf8)

    lf = lf.with_columns(method_friendly=method_expr)

    # 7. Canonical Variant ID
    if has_vid_canon:
        canon_base = pl.col("variant_id_canonical")
    else:
        canon_base = pl.lit(None)

    fallback_id = (
        pl.when(pl.col("variant_id").is_not_null())
        .then(pl.col("variant_id"))
        .otherwise(
            pl.concat_str(
                [
                    pl.col("CHROM"),
                    pl.lit(":"),
                    pl.col("POS").cast(pl.Utf8),
                    pl.lit(":"),
                    pl.col("REF"),
                    pl.lit(">"),
                    pl.col("ALT"),
                ]
            )
        )
    )

    final_canon = (
        pl.when(canon_base.is_null())
        .then(fallback_id)
        .otherwise(canon_base)
        .map_elements(canon_fn, return_dtype=pl.Utf8)
    )

    lf = lf.with_columns(variant_id_canonical=final_canon)

    # Cleanup temp columns
    lf = lf.drop("_chrom_from_id", "_pos_from_id", "_ref_from_id", "_alt_from_id")
    return lf