"""downsampling utilities for null variant data."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import polars as pl

from .normalisation_helper import UNRESOLVED_TRACK_KEY, protocol_track_exprs


def _group_columns(df: pl.DataFrame) -> list[str]:
    cols = ['gene_id']
    if 'protocol_group' in df.columns and (
        len(df) > 0 and df.filter(pl.col('protocol_group') != 'other').height > 0
    ):
        cols.append('protocol_group')
    if 'track_key' in df.columns and (
        len(df) > 0 and df.filter(pl.col('track_key') != UNRESOLVED_TRACK_KEY).height > 0
    ):
        cols.append('track_key')
    return cols


def _group_key(row: dict[str, object], group_cols: list[str]) -> str | tuple[object, ...]:
    if len(group_cols) == 1:
        return str(row[group_cols[0]])
    return tuple(row[col] for col in group_cols)


def load_real_variant_counts(
    real_path: Path,
    filter_snvs: bool = False,
) -> dict[str | tuple[object, ...], int]:
    """load per-gene variant counts from real data.
    
    auto-detects gene-level (has n_variants column) or variant-level (needs counting).
    
    args:
        real_path (Path): parquet with either gene_id + n_variants or variant-level data.
        filter_snvs (bool): if True, filter to SNVs only before counting
    
    returns:
        dict[str | tuple[object, ...], int]: group key -> n_variants mapping.
    """
    df = pl.read_parquet(real_path)
    df = df.with_columns(*protocol_track_exprs(df.columns))
    group_cols = _group_columns(df)
    
    if filter_snvs and 'REF' in df.columns and 'ALT' in df.columns:
        print(f'  Filtering real reference to SNVs only...')
        n_before = len(df)
        df = df.filter(
            (pl.col('REF').str.len_chars() == 1) &
            (pl.col('ALT').str.len_chars() == 1) &
            pl.col('REF').is_in(['A', 'C', 'G', 'T']) &
            pl.col('ALT').is_in(['A', 'C', 'G', 'T'])
        )
        n_after = len(df)
        print(f'    {n_before:,} → {n_after:,} variants ({n_before-n_after:,} removed)')
    
    if 'n_variants' in df.columns and 'gene_id' in df.columns:
        if filter_snvs:
            print("  Warning: SNV filtering requested but input is gene-level (already aggregated)")
        counts_df = df.select(group_cols + ['n_variants']).unique(group_cols)
    elif 'gene_id' in df.columns:
        counts_df = df.group_by(group_cols).agg(pl.len().alias('n_variants'))
    else:
        raise ValueError(f'missing gene_id column in {real_path}')

    return {
        _group_key(row, group_cols): int(row['n_variants'])
        for row in counts_df.iter_rows(named=True)
    }


def downsample_to_real_counts(
    null_df: pl.DataFrame,
    real_counts: dict[str | tuple[object, ...], int],
    seed: int = 42,
    verbose: bool = True,
    observed_df: pl.DataFrame | None = None,
    enforce_disjoint: bool = True,
) -> pl.DataFrame:
    """downsample null variants to match real variant counts per gene.
    
    Args:
        null_df: synthetic/null variants to downsample
        real_counts: target counts per gene from real data
        seed: random seed for reproducibility
        verbose: print progress
        observed_df: observed variants to exclude (optional)
        enforce_disjoint: remove overlaps with observed if provided
    """
    if verbose:
        print(f'  downsample input: {null_df.shape}')
        print(f'  real gene count: {len(real_counts)}')

    null_df = null_df.with_columns(*protocol_track_exprs(null_df.columns))
    group_cols = _group_columns(null_df)
    if verbose:
        print(f'  downsample groups: {group_cols}')
    
    if enforce_disjoint and observed_df is not None:
        observed_df = observed_df.with_columns(*protocol_track_exprs(observed_df.columns))
        if verbose:
            print(f'  enforcing disjointness with observed data...')
        
        key_cols = ['CHROM', 'POS', 'REF', 'ALT', 'gene_id']
        if 'protocol_group' in null_df.columns and 'protocol_group' in observed_df.columns:
            key_cols.append('protocol_group')
        if 'track_key' in null_df.columns and 'track_key' in observed_df.columns:
            key_cols.append('track_key')
        obs_keys = observed_df.select(key_cols).unique()
        
        null_before = len(null_df)
        null_df = null_df.join(obs_keys, on=key_cols, how='anti')
        null_after = len(null_df)
        
        if verbose:
            removed = null_before - null_after
            if removed > 0:
                print(f'removed {removed:,} overlapping variants ({removed/null_before*100:.2f}%)')
            else:
                print(f'no overlaps detected')
    
    pdf = null_df.to_pandas()
    sampled_parts = []
    
    for group_value, group in pdf.groupby(group_cols, dropna=False):
        group_key = group_value if isinstance(group_value, tuple) else str(group_value)
        n_req = int(real_counts.get(group_key, 0))
        
        if n_req <= 0:
            continue
        
        avail = len(group)
        if avail == 0:
            continue
        
        if avail <= n_req:
            sampled_parts.append(group)
        else:
            sort_cols = [col for col in ['variant_id', 'track_key', 'protocol_group'] if col in group.columns]
            group_sorted = group.sort_values(sort_cols) if sort_cols else group
            
            taken = group_sorted.sample(n=n_req, random_state=seed)
            sampled_parts.append(taken)
    
    if sampled_parts:
        sampled_pdf = pd.concat(sampled_parts, axis=0)
        result = pl.from_pandas(sampled_pdf)
        
        result = result.sort([col for col in ['gene_id', 'protocol_group', 'track_key', 'variant_id'] if col in result.columns])
        
        if verbose:
            print(f'  downsample output: {result.shape}')
        return result
    else:
        if verbose:
            print('  downsample output: empty')
        return pl.DataFrame(schema=null_df.schema)