"""downsampling utilities for null variant data."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import polars as pl


def load_real_variant_counts(real_path: Path, filter_snvs: bool = False) -> dict[str, int]:
    """load per-gene variant counts from real data.
    
    auto-detects gene-level (has n_variants column) or variant-level (needs counting).
    
    args:
        real_path (Path): parquet with either gene_id + n_variants or variant-level data.
        filter_snvs (bool): if True, filter to SNVs only before counting
    
    returns:
        dict[str, int]: gene_id -> n_variants mapping.
    """
    df = pl.read_parquet(real_path)
    
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
        return (
            df.select(['gene_id', 'n_variants'])
            .unique('gene_id')
            .to_pandas()
            .set_index('gene_id')['n_variants']
            .to_dict()
        )
    elif 'gene_id' in df.columns:
        counts = (
            df.group_by('gene_id')
            .agg(pl.len().alias('n_variants'))
            .to_pandas()
            .set_index('gene_id')['n_variants']
            .to_dict()
        )
        return counts
    else:
        raise ValueError(f'missing gene_id column in {real_path}')

def downsample_to_real_counts(
    null_df: pl.DataFrame,
    real_counts: dict[str, int],
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
    
    if enforce_disjoint and observed_df is not None:
        if verbose:
            print(f'  enforcing disjointness with observed data...')
        
        key_cols = ['CHROM', 'POS', 'REF', 'ALT', 'gene_id']
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
    
    for gene, group in pdf.groupby('gene_id'):
        n_req = int(real_counts.get(gene, 0))
        
        if n_req <= 0:
            continue
        
        avail = len(group)
        if avail == 0:
            continue
        
        if avail <= n_req:
            sampled_parts.append(group)
        else:
            group_sorted = group.sort_values('variant_id')
            
            taken = group_sorted.sample(n=n_req, random_state=seed)
            sampled_parts.append(taken)
    
    if sampled_parts:
        sampled_pdf = pd.concat(sampled_parts, axis=0)
        result = pl.from_pandas(sampled_pdf)
        
        result = result.sort(['gene_id', 'variant_id'])
        
        if verbose:
            print(f'  downsample output: {result.shape}')
        return result
    else:
        if verbose:
            print('  downsample output: empty')
        return pl.DataFrame(schema=null_df.schema)