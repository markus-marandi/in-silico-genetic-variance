"""deduplication utilities for variant data."""
from __future__ import annotations

from pathlib import Path

import polars as pl

from .normalisation_helper import UNRESOLVED_TRACK_KEY, protocol_track_exprs


DEDUP_KEYS = ['gene_id', 'variant_id', 'track_key']


def _ensure_protocol_identity(df: pl.DataFrame) -> pl.DataFrame:
    schema_cols = set(df.columns)
    df = df.with_columns(*protocol_track_exprs(schema_cols))

    missing_track_rows = df.filter(pl.col('track_key') == UNRESOLVED_TRACK_KEY)
    if len(missing_track_rows) > 0:
        raise ValueError(
            'track_key could not be derived for all rows; rerun scoring or rebuild the '
            'variant parquet from chunk outputs that still contain track metadata'
        )

    return df


def _print_protocol_summary(df: pl.DataFrame, label: str) -> None:
    summary = (
        df.group_by('protocol_group')
        .agg(
            pl.len().alias('n_rows'),
            pl.col('track_key').n_unique().alias('n_track_keys'),
            pl.col('variant_id').n_unique().alias('n_variant_ids'),
            pl.col('gene_id').n_unique().alias('n_genes'),
        )
        .sort('protocol_group')
    )

    print(f'  protocol counts {label}:')
    for row in summary.iter_rows(named=True):
        print(
            '    '
            f"{row['protocol_group']}: "
            f"rows={row['n_rows']:,}, "
            f"track_keys={row['n_track_keys']:,}, "
            f"variant_ids={row['n_variant_ids']:,}, "
            f"genes={row['n_genes']:,}"
        )


def deduplicate_by_gene_and_variant(
    df: pl.DataFrame,
    verbose: bool = True,
) -> pl.DataFrame:
    """deduplicate within gene, variant, and track identity."""
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
    
    df = _ensure_protocol_identity(df)
    df = df.filter(pl.col('raw_score').is_not_null())

    if verbose:
        _print_protocol_summary(df, 'before dedup')
    
    dedup = (
        df
        .with_columns(pl.col('raw_score').abs().alias('_abs_score'))
        .sort(
            ['_abs_score', 'gene_id', 'variant_id', 'track_key'],
            descending=[True, False, False, False],
        )
        .unique(subset=DEDUP_KEYS, keep='first', maintain_order=True)
        .drop('_abs_score')
    )

    dedup = dedup.sort(['gene_id', 'protocol_group', 'variant_id', 'track_key'])
    
    if verbose:
        print(f'  deduplicated output: {dedup.shape}')
        _print_protocol_summary(dedup, 'after dedup')
    
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
