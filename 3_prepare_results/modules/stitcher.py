from __future__ import annotations

from pathlib import Path
from typing import Optional

import polars as pl

# Import helper directly from the neighbor module
from .normalisation_helper import strip_ensembl_version


def load_gene_list(path: Optional[Path]) -> Optional[pl.Series]:
    if path is None or not path.exists():
        return None
    df = pl.read_csv(path, has_header=False, new_columns=["gene_id"])
    series = df["gene_id"].str.strip_chars().str.split(".").list.get(0)
    return series.filter(series != "")


def stitch_variants(chunks_dir: Path, gene_list_path: Optional[Path] = None) -> pl.LazyFrame:
    """lazy scan chunk tsv.gz files, normalize gene ids, optional whitelist."""
    strip_fn = strip_ensembl_version

    pattern = str(chunks_dir / "chunk_*.tsv.gz")
    lf = pl.scan_csv(
        pattern,
        separator="\t",
        infer_schema_length=1000,
        null_values=["", "NA", "None"],
        low_memory=True,
    )

    # DYNAMIC SCHEMA CHECK (Safe)
    # This prevents the "ColumnNotFoundError: gene_tag" crash
    # We inspect the actual file columns before building the expressions
    schema_cols = set(lf.collect_schema().names())

    has_gene_id = "gene_id" in schema_cols
    has_gene_tag = "gene_tag" in schema_cols

    # Build the gene_norm expression based ONLY on existing columns
    if has_gene_id and has_gene_tag:
        norm_expr = (
            pl.when(pl.col("gene_id").is_not_null())
            .then(pl.col("gene_id").map_elements(strip_fn, return_dtype=pl.Utf8))
            .when(pl.col("gene_tag").is_not_null())
            .then(pl.col("gene_tag").map_elements(strip_fn, return_dtype=pl.Utf8))
            .otherwise(None)
        )
    elif has_gene_id:
        norm_expr = pl.col("gene_id").map_elements(strip_fn, return_dtype=pl.Utf8)
    elif has_gene_tag:
        norm_expr = pl.col("gene_tag").map_elements(strip_fn, return_dtype=pl.Utf8)
    else:
        # Fallback (should not happen based on your files)
        norm_expr = pl.lit(None).cast(pl.Utf8)

    lf = lf.with_columns(gene_norm=norm_expr)

    genes = load_gene_list(gene_list_path)
    if genes is not None:
        lf = lf.filter(pl.col("gene_norm").is_in(genes))

    return lf