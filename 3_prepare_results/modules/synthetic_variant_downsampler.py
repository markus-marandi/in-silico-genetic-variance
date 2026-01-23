"""downsampling utilities for null variant data."""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import polars as pl


def load_real_variant_counts(real_path: Path) -> dict[str, int]:
    """load per-gene variant counts from real data.
    
    auto-detects gene-level (has n_variants column) or variant-level (needs counting).
    
    args:
        real_path (Path): parquet with either gene_id + n_variants or variant-level data.
    
    returns:
        dict[str, int]: gene_id -> n_variants mapping.
    """
    df = pl.read_parquet(real_path)
    
    if 'n_variants' in df.columns and 'gene_id' in df.columns:
        # gene-level file with pre-computed counts
        return (
            df.select(['gene_id', 'n_variants'])
            .unique('gene_id')
            .to_pandas()
            .set_index('gene_id')['n_variants']
            .to_dict()
        )
    elif 'gene_id' in df.columns:
        # variant-level file, count per gene
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
) -> pl.DataFrame:
    """downsample null variants to match real variant counts per gene.
    
    args:
        null_df (pl.DataFrame): null variant dataframe with gene_id column.
        real_counts (dict[str, int]): gene_id -> n_variants mapping from real data.
        seed (int): random seed for reproducibility.
        verbose (bool): print diagnostic info.
    
    returns:
        pl.DataFrame: downsampled null variants matching real counts per gene.
    """
    if verbose:
        print(f'  downsample input: {null_df.shape}')
        print(f'  real gene count: {len(real_counts)}')
    
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
            # take all available
            sampled_parts.append(group)
        else:
            # sample without replacement
            taken = group.sample(n=n_req, random_state=seed)
            sampled_parts.append(taken)
    
    if sampled_parts:
        sampled_pdf = pd.concat(sampled_parts, axis=0)
        result = pl.from_pandas(sampled_pdf)
        if verbose:
            print(f'  downsample output: {result.shape}')
        return result
    else:
        if verbose:
            print('  downsample output: empty')
        return pl.DataFrame(schema=null_df.schema)
