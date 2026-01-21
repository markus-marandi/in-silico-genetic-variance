#!/usr/bin/env python3
"""downsample synthetic variants to match real variant counts per gene.

args:
    synthetic (str): synthetic variants tsv with gene_id column
    real (str): real variants tsv with gene_id column
    output (str): output path for downsampled variants
    seed (int): random seed for reproducibility

returns:
    None
"""

import argparse
from pathlib import Path
from typing import Dict

import pandas as pd


def count_variants_per_gene(df: pd.DataFrame, gene_col: str = "gene_id") -> Dict[str, int]:
    """count variants per gene."""
    return df[gene_col].value_counts().to_dict()


def downsample_to_real_counts(
    synthetic_df: pd.DataFrame,
    real_counts: Dict[str, int],
    gene_col: str = "gene_id",
    seed: int = 42,
) -> pd.DataFrame:
    """downsample synthetic variants to match real variant counts per gene.
    
    args:
        synthetic_df (pd.DataFrame): synthetic variants with gene_id column.
        real_counts (Dict[str, int]): count of real variants per gene.
        gene_col (str): column name for gene identifier.
        seed (int): random seed for reproducibility.
    
    returns:
        pd.DataFrame: downsampled synthetic variants.
    """
    sampled_parts = []
    
    for gene, group in synthetic_df.groupby(gene_col):
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
        return pd.concat(sampled_parts, axis=0).reset_index(drop=True)
    else:
        return pd.DataFrame(columns=synthetic_df.columns)


def detect_gene_column(df: pd.DataFrame) -> str:
    """detect gene id column name."""
    candidates = ["gene_id", "gene_tag", "gene", "name", "ensg", "ensembl_gene"]
    for c in candidates:
        if c in df.columns:
            return c
        # case insensitive
        for col in df.columns:
            if col.lower() == c.lower():
                return col
    raise ValueError(f"could not find gene id column in: {df.columns.tolist()}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="downsample synthetic variants to match real counts per gene"
    )
    parser.add_argument(
        "--synthetic", required=True, help="synthetic variants tsv path"
    )
    parser.add_argument(
        "--real", required=True, help="real variants tsv path"
    )
    parser.add_argument(
        "--output", required=True, help="output path for downsampled variants"
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="random seed (default: 42)"
    )
    parser.add_argument(
        "--gene-col", help="gene id column name (auto-detected if not specified)"
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    
    synthetic_path = Path(args.synthetic)
    real_path = Path(args.real)
    output_path = Path(args.output)
    
    if not synthetic_path.exists():
        raise FileNotFoundError(f"synthetic variants not found: {synthetic_path}")
    if not real_path.exists():
        raise FileNotFoundError(f"real variants not found: {real_path}")
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # detect compression from extension
    synthetic_compression = "gzip" if str(synthetic_path).endswith(".gz") else None
    real_compression = "gzip" if str(real_path).endswith(".gz") else None
    output_compression = "gzip" if str(output_path).endswith(".gz") else None
    
    print(f"loading synthetic variants from {synthetic_path}...")
    synthetic_df = pd.read_csv(synthetic_path, sep="\t", compression=synthetic_compression)
    print(f"loaded {len(synthetic_df)} synthetic variants")
    
    print(f"loading real variants from {real_path}...")
    real_df = pd.read_csv(real_path, sep="\t", compression=real_compression)
    print(f"loaded {len(real_df)} real variants")
    
    # detect gene column
    syn_gene_col = args.gene_col or detect_gene_column(synthetic_df)
    
    # for real variants, need to infer gene from the file structure
    # if real variants have gene_id column, use it
    # otherwise, if they came from per_gene export, we need to use the bed file
    # for now, assume real variants have gene_id or similar column
    try:
        real_gene_col = args.gene_col or detect_gene_column(real_df)
    except ValueError:
        # real variants might not have gene_id; need to add it from bed intervals
        # for now, create a dummy column if not present
        print("warning: real variants missing gene_id column, cannot match per-gene")
        print("outputting full synthetic set without downsampling")
        synthetic_df.to_csv(output_path, sep="\t", index=False, compression=output_compression)
        return
    
    print(f"using gene column: synthetic={syn_gene_col}, real={real_gene_col}")
    
    # count real variants per gene
    real_counts = count_variants_per_gene(real_df, gene_col=real_gene_col)
    n_genes = len(real_counts)
    total_real = sum(real_counts.values())
    print(f"real variants: {total_real} across {n_genes} genes")
    
    # downsample synthetic to match real counts
    print("downsampling synthetic variants...")
    downsampled_df = downsample_to_real_counts(
        synthetic_df,
        real_counts,
        gene_col=syn_gene_col,
        seed=args.seed,
    )
    
    print(f"downsampled to {len(downsampled_df)} variants")
    
    # write output
    print(f"writing output to {output_path}...")
    downsampled_df.to_csv(output_path, sep="\t", index=False, compression=output_compression)
    
    # summary stats
    syn_genes = synthetic_df[syn_gene_col].nunique()
    down_genes = downsampled_df[syn_gene_col].nunique() if len(downsampled_df) > 0 else 0
    
    print("\nsummary:")
    print(f"  synthetic variants: {len(synthetic_df)} ({syn_genes} genes)")
    print(f"  real variants: {total_real} ({n_genes} genes)")
    print(f"  downsampled: {len(downsampled_df)} ({down_genes} genes)")
    
    if len(downsampled_df) < total_real:
        deficit = total_real - len(downsampled_df)
        print(f"  note: {deficit} fewer variants due to insufficient synthetic coverage")


if __name__ == "__main__":
    main()
