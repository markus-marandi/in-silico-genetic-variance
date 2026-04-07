"""module for adding gene-aware permuted AFs to variants."""
from __future__ import annotations

import logging
from pathlib import Path
import numpy as np
import polars as pl

log = logging.getLogger(__name__)

def load_gene_af_pools(source: Path | pl.DataFrame) -> dict[str, np.ndarray]:
    """extract non-zero AFs per gene from source file or DataFrame."""
    
    if isinstance(source, (str, Path)):
        log.info('loading per-gene AF pools from file: %s', source)
        df = pl.read_parquet(source, columns=["gene_id", "AF"])
    else:
        log.info('loading per-gene AF pools from in-memory DataFrame')
        # Select only necessary columns to avoid overhead
        df = source.select(["gene_id", "AF"])
        
    # normalize ids and keep full vector lengths by mapping null AF -> 0.0
    df = (
        df.with_columns(
            pl.col('gene_id').cast(pl.Utf8).str.split('.').list.get(0).alias('gene_id'),
            pl.col('AF').cast(pl.Float64, strict=False).fill_null(0.0).alias('AF')
        )
    )
    
    if len(df) == 0:
        raise ValueError('no valid AFs found in source')
    
    # group by gene and collect AFs
    gene_af_pools = {}
    
    for gene_df in df.partition_by('gene_id', as_dict=False):
        gene_id = gene_df[0, "gene_id"]
        afs = gene_df['AF'].to_numpy()
        gene_af_pools[gene_id] = afs
    
    log.info('loaded AF pools for %d genes', len(gene_af_pools))
    return gene_af_pools


def add_perm_af_gene_aware(
    df: pl.DataFrame,
    gene_af_pools: dict[str, np.ndarray],
    seed: int = 42,
    strict_no_replacement: bool = True,
) -> pl.DataFrame:
    """add perm_AF column by sampling from the same gene's AF pool.
    
    uses replace=False (shuffle) when possible to preserve exact distribution.
    falls back to replace=True only when pool size < needed variants.
    """
    log.info('adding perm_AF to %d variants...', len(df))
    
    # clean up potential collision columns
    cols_to_drop = [c for c in ['AF_x', 'AF_y', 'perm_AF'] if c in df.columns]
    if cols_to_drop:
        df = df.drop(cols_to_drop)
    
    if 'gene_id' not in df.columns:
        raise ValueError('input dataframe must contain gene_id')

    # preserve original order while grouping by normalized gene id
    df = (
        df.with_row_index('_row_idx')
        .with_columns(pl.col('gene_id').cast(pl.Utf8).str.split('.').list.get(0).alias('_gene_norm'))
        .sort('_gene_norm')
    )
    
    rng = np.random.default_rng(seed)
    perm_parts: list[pl.DataFrame] = []
    
    # iterate over genes in the target dataframe
    for gene_df in df.partition_by('_gene_norm', maintain_order=True):
        gene_id = gene_df[0, '_gene_norm']
        n_needed = len(gene_df)
        
        # determine which pool to use
        pool = gene_af_pools.get(gene_id)
        if pool is None or len(pool) == 0:
            local_af = np.array([], dtype=float)
            if 'AF' in gene_df.columns:
                local_af = gene_df.with_columns(pl.col('AF').cast(pl.Float64, strict=False).fill_null(0.0))['AF'].to_numpy()
            if len(local_af) > 0:
                pool = local_af
                log.warning('gene %s missing in AF pools; using local AFs from target dataframe', gene_id)
            else:
                raise ValueError(f'gene {gene_id} has no AF pool and no local AF values to shuffle')
        
        n_available = len(pool)
        
        # prefer shuffling/subsetting without replacement to preserve gene-level AF distribution
        if n_available >= n_needed:
            sampled = rng.permutation(pool)[:n_needed]
        else:
            msg = f'gene {gene_id}: need {n_needed} AFs but pool has {n_available}'
            if strict_no_replacement:
                raise ValueError(msg + '; adjust downsampling or disable strict mode')
            log.warning(msg + '; sampling with replacement')
            sampled = rng.choice(pool, size=n_needed, replace=True)

        perm_parts.append(pl.DataFrame({'_row_idx': gene_df['_row_idx'], 'perm_AF': sampled}))

    perm_df = pl.concat(perm_parts) if perm_parts else pl.DataFrame({'_row_idx': [], 'perm_AF': []})

    return (
        df.join(perm_df, on='_row_idx', how='left')
        .sort('_row_idx')
        .drop(['_row_idx', '_gene_norm'])
    )