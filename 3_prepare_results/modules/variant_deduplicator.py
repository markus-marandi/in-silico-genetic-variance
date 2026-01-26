"""deduplication utilities for variant data."""
from __future__ import annotations

from pathlib import Path

import polars as pl


def deduplicate_by_gene_and_variant(
    df: pl.DataFrame,
    verbose: bool = True,
) -> pl.DataFrame:
    """deduplicate by gene_id and variant_id, keeping row with max absolute score."""
    if verbose:
        print(f'  dedup input: {df.shape}')
    
    schema_cols = df.columns
    if 'raw_score' in schema_cols:
        score_col = 'raw_score'
    elif 'score' in schema_cols:
        score_col = 'score'
    else:
        raise ValueError('missing raw_score/score column')
    
    if score_col != 'raw_score':
        df = df.rename({score_col: 'raw_score'})
    
    df = df.filter(pl.col('raw_score').is_not_null())
    
    dedup = (
        df
        .with_columns(pl.col('raw_score').abs().alias('_abs_score'))
        .sort('_abs_score', descending=True)
        .unique(subset=['gene_id', 'variant_id'], keep='first', maintain_order=True)
        .drop('_abs_score')
    )

    dedup = dedup.sort(['gene_id', 'variant_id'])
    
    if verbose:
        print(f'  Deduplicated output: {dedup.shape}')
    
    return dedup


def deduplicate_from_parquet(
    path: Path,
    label: str | None = None,
    verbose: bool = True,
) -> pl.DataFrame:
    """load and deduplicate variants from parquet."""
    if verbose and label:
        print(f'  loading {label} from {path}')
    
    df = pl.read_parquet(path)
    return deduplicate_by_gene_and_variant(df, verbose=verbose)
