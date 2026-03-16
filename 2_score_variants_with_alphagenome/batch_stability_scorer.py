"""
Batch stability scorer: query same 128 variants 50 times to assess model consistency.

Each iteration:
  - Scores the same batch of 128 variants
  - Deduplicates: keeps MAX abs score as raw_score
  - Collects all duplicates as within-batch confidence interval data
  - Saves iteration results
  - Rests 60 seconds before next iteration

After all 50 iterations:
  - Aggregates across runs to compute between-iteration statistics
  - Reports mean ± stdev, percentile CIs (2.5th, 97.5th for 95% CI), min/max, count
"""

from __future__ import annotations

import os
import sys
import time
import gzip
import shutil
import importlib.util
from pathlib import Path
from datetime import datetime
from typing import Any
import math

import pandas as pd
import numpy as np
import polars as pl
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome.models import dna_client as _dc


VAR_TSV = "/cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_variants.tsv"
TSS_BED = "/cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_gene_set±10kb.bed"
GENE_LIST_PATH = "/cfs/klemming/scratch/m/mmarandi/intermediate_input/dataset5_null/background_gene_set_380.tsv"

REPO_ROOT = Path(__file__).resolve().parent.parent
PATH_MANAGER = Path("/cfs/klemming/home/m/mmarandi/lab/users/mmarandi/alphagenome/in-silico-genetic-variance/helpers/path_manager.py")

BATCH_SIZE = 128
NUM_ITERATIONS = 50
REST_SECONDS = 60
SEQ_LEN = 1048576  # 1 MB context window
RENAME_FRIENDLY = {
    "GeneMaskLFCScorer": "gene_exonmask_delta_log2",
}

os.environ.setdefault("DATASET_ID", "dataset_4")
os.environ.setdefault("SAMPLE_ID", "background_batch_stability")

def _load_layout():
    if not PATH_MANAGER.exists():
        raise FileNotFoundError(f"path manager not found: {PATH_MANAGER}")
    
    spec = importlib.util.spec_from_file_location("ag_path_manager", PATH_MANAGER)
    if spec is None or spec.loader is None:
        raise ImportError(f"unable to load path manager from {PATH_MANAGER}")
    
    module = importlib.util.module_from_spec(spec)
    sys.modules["ag_path_manager"] = module 
    spec.loader.exec_module(module)
    return module.ProjectLayout.from_env()

layout = _load_layout()
layout.make_dirs()

STABILITY_DIR = layout.results_dir / "batch_stability"
STABILITY_DIR.mkdir(exist_ok=True, parents=True)
ITERATIONS_DIR = STABILITY_DIR / "iterations"
ITERATIONS_DIR.mkdir(exist_ok=True, parents=True)

print(f"Output directory: {STABILITY_DIR}")
print(f"Iterations directory: {ITERATIONS_DIR}")

load_dotenv()
API_KEY = os.getenv("API_KEY_PERSONAL")
assert API_KEY, "Set API_KEY_PERSONAL in your .env"

def ensg_core(x: str) -> str | None:
    """Extract ENSG ID without version."""
    if x is None:
        return None
    import re
    m = re.search(r'(ENSG\d+)', str(x))
    return m.group(1) if m else None

def _interval_to_str(x):
    """Convert interval to string representation."""
    if isinstance(x, str):
        return x
    if hasattr(x, "chromosome") and hasattr(x, "start") and hasattr(x, "end"):
        return f"{x.chromosome}:{int(x.start)}-{int(x.end)}:."
    return str(x)

def _normalize_label(s):
    """Normalize labels for matching."""
    import re
    s = '' if s is None else str(s)
    return re.sub(r'[^a-z0-9]+', '_', s.lower()).strip('_')

def _load_gene_whitelist(path: str | None) -> set[str]:
    """Load gene whitelist from TSV."""
    if not path:
        return set()
    genes = pd.read_csv(path, header=None, names=["gene_id"], dtype=str)
    genes["gene_id"] = genes["gene_id"].astype(str).str.strip()
    genes = genes[genes["gene_id"] != ""]
    genes["gene_id_core"] = genes["gene_id"].map(ensg_core)
    whitelist = set(genes["gene_id_core"].dropna())
    print(f"Loaded {len(whitelist)} genes from {path}")
    return whitelist

def make_intervals(seq_len: int, chroms, poses):
    """Create genome intervals for each variant."""
    half = seq_len // 2
    poses = np.asarray(poses, dtype=int)
    starts = np.maximum(poses - half, 1)
    ends = starts + seq_len
    return [
        genome.Interval(chromosome=str(c), start=int(s), end=int(e))
        for c, s, e in zip(chroms, starts, ends)
    ]

def make_variants(chroms, poses, refs, alts, names):
    """Create genome variant objects."""
    return [
        genome.Variant(
            chromosome=str(c),
            position=int(p),
            reference_bases=str(rf),
            alternate_bases=str(al),
            name=str(nm)
        )
        for c, p, rf, al, nm in zip(chroms, poses, refs, alts, names)
    ]

def normalize_tidy(tidy: pd.DataFrame) -> pd.DataFrame:
    """Normalize tidy scores output."""
    tidy = tidy.copy()
    if "variant_id" in tidy.columns:
        tidy["variant_id"] = tidy["variant_id"].map(lambda v: v.name if hasattr(v, "name") else str(v))
    elif "variant" in tidy.columns:
        tidy["variant_id"] = tidy["variant"].map(lambda v: v.name if hasattr(v, "name") else str(v))
    return tidy

def map_friendly(name: str) -> str:
    """Map scorer name to friendly version."""
    s = str(name)
    for k, v in RENAME_FRIENDLY.items():
        if k in s:
            return v
    return s

def inject_into_anndata(scores, meta_df):
    """Inject variant metadata into AnnData objects."""
    for per_var, (_, row) in zip(scores, meta_df.iterrows()):
        if not isinstance(per_var, (list, tuple)):
            per_var = [per_var]
        for ad in per_var:
            ad.obs["variant_id"] = str(row.variant_id)
            ad.obs["CHROM"] = str(row.CHROM)
            ad.obs["POS"] = int(row.POS)
            ad.obs["REF"] = str(row.REF)
            ad.obs["ALT"] = str(row.ALT)
            ad.obs["gene_tag"] = str(row.gene_tag) if pd.notna(row.gene_tag) else ""

def scores_to_df(scores, meta_df) -> pd.DataFrame:
    """Flatten scores into tidy rows."""
    rows: list[dict[str, object]] = []

    for per_variant, (_, mrow) in zip(scores, meta_df.iterrows()):
        if not isinstance(per_variant, (list, tuple)):
            per_variant = [per_variant]

        for ad_obj in per_variant:
            X = np.asarray(ad_obj.X)
            if X.ndim == 1:
                X = X[None, :]

            n_genes, n_tracks = X.shape
            interval = ad_obj.uns.get("scored_interval", None)
            iv_str = _interval_to_str(interval) if interval is not None else pd.NA
            scorer_name = str(ad_obj.uns.get("variant_scorer", ""))
            out_type = str(ad_obj.uns.get("output_type", ""))

            track_meta = getattr(ad_obj, "var", pd.DataFrame(index=range(n_tracks)))

            for g_idx in range(n_genes):
                gene_row = ad_obj.obs.iloc[g_idx] if g_idx < ad_obj.obs.shape[0] else {}
                for t_idx in range(n_tracks):
                    val = float(X[g_idx, t_idx])
                    trow = track_meta.iloc[t_idx] if t_idx < len(track_meta) else {}

                    rows.append({
                        "variant_id": str(mrow.variant_id),
                        "scored_interval_str": iv_str,
                        "output_type": out_type,
                        "variant_scorer": scorer_name,
                        "track_name": trow.get("name", pd.NA),
                        "track_strand": trow.get("strand", pd.NA),
                        "Assay_title": trow.get("assay_title", pd.NA),
                        "ontology_curie": trow.get("ontology_curie", pd.NA),
                        "biosample_name": trow.get("biosample_name", pd.NA),
                        "biosample_type": trow.get("biosample_type", pd.NA),
                        "gtex_tissue": trow.get("gtex_tissue", pd.NA),
                        "raw_score": val,
                        "gene_id": gene_row.get("gene_id", pd.NA),
                        "gene_name": gene_row.get("gene_name", pd.NA),
                        "gene_type": gene_row.get("gene_type", pd.NA),
                        "gene_strand": gene_row.get("gene_strand", pd.NA),
                        "CHROM": mrow.CHROM,
                        "POS": int(mrow.POS),
                        "REF": mrow.REF,
                        "ALT": mrow.ALT,
                        "gene_tag": mrow.gene_tag,
                    })

    return pd.DataFrame(rows)

def deduplicate_with_ci(
    df: pd.DataFrame,
    iteration: int,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Deduplicate by gene and variant, keeping MAX abs score as raw_score.
    Return both the deduplicated frame and a CI frame with all duplicates.
    
    Returns:
        (dedup_df, ci_df): deduplicated main results and raw CI data from duplicates
    """
    if verbose:
        print(f"  dedup input: {df.shape}")
    
    # Ensure we have a score column
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
    
    # Convert to polars for efficient deduplication
    df_pl = pl.from_pandas(df) if isinstance(df, pd.DataFrame) else df
    
    # Add abs score and sort by it (descending) so first is max
    df_pl = df_pl.with_columns(pl.col('raw_score').abs().alias('_abs_score'))
    df_pl = df_pl.sort('_abs_score', descending=True)
    
    # Keep first occurrence (max abs score) per gene-variant pair
    dedup = (
        df_pl
        .unique(subset=['gene_id', 'variant_id'], keep='first', maintain_order=True)
        .drop('_abs_score')
        .sort(['gene_id', 'variant_id'])
        .to_pandas()
    )
    
    # Collect ALL raw scores (including the ones we kept) as CI data per variant
    ci_data = df_pl.select([
        'gene_id', 'variant_id', 'raw_score', '_abs_score'
    ]).drop('_abs_score').to_pandas()
    ci_data['iteration'] = iteration
    
    if verbose:
        print(f"  deduplicated output: {dedup.shape}")
        print(f"  CI data rows (all scores): {ci_data.shape}")
    
    return dedup, ci_data

# ============================================================================
# MAIN EXECUTION
# ============================================================================

print("=" * 80)
print("BATCH STABILITY SCORER: 50x same-variant queries with CI estimation")
print("=" * 80)

# Load variants
print(f"\nLoading variants from {VAR_TSV}...")
df = pd.read_csv(VAR_TSV, sep="\t", compression="infer", low_memory=False)

need_base = {"CHROM", "POS", "REF", "ALT", "gene_tag"}
missing = need_base - set(df.columns)
if missing:
    raise ValueError(f"Missing columns in {VAR_TSV}: {missing}")

if "variant_id" not in df.columns:
    df["variant_id"] = (
        df["CHROM"].astype(str)
        + ":" + pd.to_numeric(df["POS"], errors="coerce").astype("Int64").astype(str)
        + ":" + df["REF"].astype(str) + ">" + df["ALT"].astype(str)
    )

df = df.loc[:, ["variant_id", "gene_tag", "CHROM", "POS", "REF", "ALT"]].copy()
df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
df = df.dropna(subset=["POS"]).copy()
df["POS"] = df["POS"].astype(int)

# Take only first BATCH_SIZE variants
query_batch = df.iloc[:BATCH_SIZE].copy()
print(f"Selected first {len(query_batch)} variants for repeated queries")

# Initialize API
print(f"\nInitializing AlphaGenome client...")
t0_client = time.time()
model = dna_client.create(api_key=API_KEY)
print(f"API client ready in {time.time() - t0_client:.1f}s")

RNA = dna_client.OutputType.RNA_SEQ
ORG = _dc.Organism.HOMO_SAPIENS

# Load gene whitelist
gene_whitelist = _load_gene_whitelist(GENE_LIST_PATH)

# Collect all CI data across iterations
all_ci_data = []


for iteration in range(1, NUM_ITERATIONS + 1):
    print(f"\n{'=' * 80}")
    print(f"ITERATION {iteration}/{NUM_ITERATIONS}")
    print(f"{'=' * 80}")
    
    start_iter = time.time()
    
    # Prepare batch
    meta = query_batch.copy()
    intervals = make_intervals(SEQ_LEN, meta["CHROM"].to_numpy(), meta["POS"].to_numpy())
    variants = make_variants(meta["CHROM"], meta["POS"], meta["REF"], meta["ALT"], meta["variant_id"])
    
    # Score
    t0_score = time.time()
    print(f"Scoring {len(meta)} variants...")
    scores = model.score_variants(
        intervals=intervals,
        variants=variants,
        variant_scorers=[variant_scorers.GeneMaskLFCScorer(requested_output=RNA)],
        organism=ORG,
        progress_bar=False,
    )
    elapsed_score = time.time() - t0_score
    print(f"Scoring done in {elapsed_score:.1f}s")
    
    # Convert to dataframe
    inject_into_anndata(scores, meta)
    tidy = variant_scorers.tidy_scores(scores, match_gene_strand=False)
    tidy = normalize_tidy(tidy)
    tidy["scorer_friendly"] = tidy["variant_scorer"].map(map_friendly)
    
    # Convert to polars for deduplication
    tidy_pl = pl.from_pandas(tidy)
    
    # Deduplicate with CI tracking
    dedup_df, ci_df = deduplicate_with_ci(tidy_pl, iteration=iteration, verbose=True)
    
    # Filter by whitelist
    if gene_whitelist:
        dedup_df["gene_id_core"] = dedup_df["gene_id"].map(ensg_core)
        dedup_df = dedup_df[dedup_df["gene_id_core"].isin(gene_whitelist)].copy()
        dedup_df = dedup_df.drop(columns=["gene_id_core"], errors="ignore")
    
    # Save iteration results
    iter_out = ITERATIONS_DIR / f"iteration_{iteration:02d}.tsv.gz"
    dedup_df.to_csv(iter_out, sep="\t", index=False, compression="gzip")
    print(f"Saved iteration {iteration} results to {iter_out.name}")
    
    # Accumulate CI data
    all_ci_data.append(ci_df)
    
    # Rest before next iteration (except after last)
    if iteration < NUM_ITERATIONS:
        print(f"Resting for {REST_SECONDS} seconds before next iteration...")
        time.sleep(REST_SECONDS)
    
    elapsed_iter = time.time() - start_iter
    print(f"Iteration {iteration} completed in {elapsed_iter:.1f}s")

print(f"\n{'=' * 80}")
print("ALL ITERATIONS COMPLETE")
print(f"{'=' * 80}")

print(f"\nAggregating results across {NUM_ITERATIONS} iterations...")

# Combine all CI data
all_ci = pd.concat(all_ci_data, ignore_index=True)
print(f"Total raw scores collected: {len(all_ci)}")

# Group by variant and gene to compute statistics
agg_stats = (
    all_ci
    .groupby(['variant_id', 'gene_id'])
    .agg({
        'raw_score': [
            'count',
            'mean',
            'std',
            'min',
            'max',
            ('p025', lambda x: x.quantile(0.025)),
            ('p25', lambda x: x.quantile(0.25)),
            ('p50', lambda x: x.quantile(0.50)),
            ('p75', lambda x: x.quantile(0.75)),
            ('p975', lambda x: x.quantile(0.975)),
        ]
    })
    .reset_index()
)

# Flatten column names
agg_stats.columns = ['_'.join(col).strip('_') if col[1] else col[0] 
                     for col in agg_stats.columns.values]
agg_stats = agg_stats.rename(columns={
    'raw_score_count': 'count',
    'raw_score_mean': 'mean_score',
    'raw_score_std': 'std_score',
    'raw_score_min': 'min_score',
    'raw_score_max': 'max_score',
    'raw_score_p025': 'ci95_lower',
    'raw_score_p25': 'ci50_lower',
    'raw_score_p50': 'median_score',
    'raw_score_p75': 'ci50_upper',
    'raw_score_p975': 'ci95_upper',
})

# Save aggregated statistics
agg_out = STABILITY_DIR / "aggregated_statistics.tsv.gz"
agg_stats.to_csv(agg_out, sep="\t", index=False, compression="gzip")
print(f"Saved aggregated statistics to {agg_out.name}")
print(f"Records: {len(agg_stats)}")

# Save raw CI data for post-hoc analysis
raw_ci_out = STABILITY_DIR / "all_raw_ci_data.tsv.gz"
all_ci.to_csv(raw_ci_out, sep="\t", index=False, compression="gzip")
print(f"Saved raw CI data to {raw_ci_out.name}")

# Print summary
print(f"Total unique variants scored: {len(agg_stats)}")
print(f"Total iterations: {NUM_ITERATIONS}")
print(f"Total raw scores collected: {len(all_ci)}")
print(f"Mean observations per variant: {agg_stats['count'].mean():.1f}")
print(f"Min observations per variant: {agg_stats['count'].min()}")
print(f"Max observations per variant: {agg_stats['count'].max()}")

print(f"\nMean score across all variants:")
print(f"  Mean:   {agg_stats['mean_score'].mean():.6f}")
print(f"  Std:    {agg_stats['mean_score'].std():.6f}")
print(f"  Range:  [{agg_stats['mean_score'].min():.6f}, {agg_stats['mean_score'].max():.6f}]")

print(f"\nScore variability (std within variant):")
print(f"  Mean std: {agg_stats['std_score'].mean():.6f}")
print(f"  Median:   {agg_stats['std_score'].median():.6f}")
print(f"  Max:      {agg_stats['std_score'].max():.6f}")

print("Output files:")
print(f"  Iterations: {ITERATIONS_DIR}")
print(f"  Aggregated: {agg_out}")
print(f"  Raw CI data: {raw_ci_out}")
