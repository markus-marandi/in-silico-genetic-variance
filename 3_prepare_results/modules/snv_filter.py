"""Filter to keep only single nucleotide variants (SNVs).

This module removes:
- Insertions and deletions (indels)
- Multi-nucleotide variants (MNVs)
- Any variants where REF or ALT is not a single base

This ensures clean SNV-only datasets for downstream analysis.
"""

import polars as pl


def is_valid_snv(ref: str | None, alt: str | None) -> bool:
    """Check if a variant is a valid single nucleotide variant.
    
    Args:
        ref: Reference allele
        alt: Alternate allele
        
    Returns:
        True if both are single bases (A, C, G, T), False otherwise
    """
    if ref is None or alt is None:
        return False
    if len(ref) != 1 or len(alt) != 1:
        return False
    if ref not in ['A', 'C', 'G', 'T'] or alt not in ['A', 'C', 'G', 'T']:
        return False
    return True


def filter_to_snvs(df: pl.DataFrame, verbose: bool = True) -> pl.DataFrame:
    """Filter dataframe to keep only single nucleotide variants.
    
    Args:
        df: Input dataframe with REF and ALT columns
        verbose: If True, print filtering statistics
        
    Returns:
        Filtered dataframe containing only SNVs
        
    Raises:
        ValueError: If REF or ALT columns are missing
    """
    if 'REF' not in df.columns or 'ALT' not in df.columns:
        raise ValueError("DataFrame must contain 'REF' and 'ALT' columns")
    
    n_before = len(df)
    
    # Filter to single base substitutions
    df_filtered = df.filter(
        (pl.col('REF').str.len_chars() == 1) &
        (pl.col('ALT').str.len_chars() == 1) &
        pl.col('REF').is_in(['A', 'C', 'G', 'T']) &
        pl.col('ALT').is_in(['A', 'C', 'G', 'T'])
    )
    
    n_after = len(df_filtered)
    n_removed = n_before - n_after
    
    if verbose:
        print(f"  SNV filtering: {n_before:,} → {n_after:,} variants")
        print(f"  Removed {n_removed:,} non-SNV variants ({n_removed/n_before*100:.2f}%)")
        
        if n_removed > 0:
            # Show examples of removed variants
            removed = df.filter(
                ~(
                    (pl.col('REF').str.len_chars() == 1) &
                    (pl.col('ALT').str.len_chars() == 1) &
                    pl.col('REF').is_in(['A', 'C', 'G', 'T']) &
                    pl.col('ALT').is_in(['A', 'C', 'G', 'T'])
                )
            )
            
            # Count by type
            removed_types = removed.group_by(
                pl.struct(['REF', 'ALT'])
            ).agg(
                pl.count().alias('n')
            ).sort('n', descending=True).head(10)
            
            print(f"  Top 10 removed variant types:")
            for row in removed_types.iter_rows(named=True):
                ref = row['REF']
                alt = row['ALT']
                n = row['n']
                ref_len = len(ref) if ref else 0
                alt_len = len(alt) if alt else 0
                
                if ref_len > 1 or alt_len > 1:
                    vtype = "MNV/Indel"
                elif ref_len == 0 or alt_len == 0:
                    vtype = "Missing"
                else:
                    vtype = "Invalid base"
                    
                print(f"    {ref} → {alt} ({vtype}): {n:,}")
    
    return df_filtered


def filter_snvs_lazy(lf: pl.LazyFrame) -> pl.LazyFrame:
    """Filter lazy frame to keep only single nucleotide variants.
    
    Args:
        lf: Input lazy frame with REF and ALT columns
        
    Returns:
        Filtered lazy frame containing only SNVs
    """
    return lf.filter(
        (pl.col('REF').str.len_chars() == 1) &
        (pl.col('ALT').str.len_chars() == 1) &
        pl.col('REF').is_in(['A', 'C', 'G', 'T']) &
        pl.col('ALT').is_in(['A', 'C', 'G', 'T'])
    )
