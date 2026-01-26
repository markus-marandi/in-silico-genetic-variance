"""module for adding gene-aware permuted AFs to variants."""
from __future__ import annotations

import logging
from pathlib import Path
import numpy as np
import polars as pl

log = logging.getLogger(__name__)

def load_gene_af_pools(perm_variants_path: Path) -> dict[str, np.ndarray]:
    """extract non-zero AFs per gene from source file."""
    log.info('loading per-gene AF pools from %s', perm_variants_path)
    df = pl.read_parquet(perm_variants_path, columns=["gene_id", "AF"])
    
    # filter to non-null AFs (include AF=0)
    df = df.filter(pl.col('AF').is_not_null())
    
    if len(df) == 0:
        raise ValueError(f'no valid AFs found in {perm_variants_path}')
    
    # group by gene and collect AFs
    # conversion to python dict of numpy arrays is efficient for sampling loop
    gene_af_pools = {}
    
    # Use partition_by for speed if dataset fits in memory, otherwise loop
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
) -> pl.DataFrame:
    """add perm_AF column by sampling from the same gene's AF pool."""
    log.info('adding perm_AF to %d variants...', len(df))
    
    # clean up potential collision columns
    cols_to_drop = [c for c in ['AF_x', 'AF_y', 'perm_AF'] if c in df.columns]
    if cols_to_drop:
        df = df.drop(cols_to_drop)
    
    # sort by gene_id to ensure group order matches loop
    df = df.sort('gene_id')
    
    # create global fallback pool (rare case where gene has variants but no AFs in pool)
    global_pool = np.concatenate(list(gene_af_pools.values()))
    
    rng = np.random.default_rng(seed)
    perm_afs = []
    
    # iterate over genes in the target dataframe
    for gene_df in df.partition_by('gene_id', maintain_order=True):
        gene_id = gene_df[0, "gene_id"]
        n_needed = len(gene_df)
        
        if gene_id in gene_af_pools:
            # Sample with replacement from specific gene pool
            gene_pool = gene_af_pools[gene_id]
            sampled = rng.choice(gene_pool, size=n_needed, replace=True)
        else:
            # Fallback to global pool
            sampled = rng.choice(global_pool, size=n_needed, replace=True)
        
        perm_afs.append(sampled)
    
    perm_afs_array = np.concatenate(perm_afs)
    
    df = df.with_columns(pl.Series('perm_AF', perm_afs_array))
    
    return df