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
    - Reports mean +/- stdev, percentile CIs (2.5th, 97.5th for 95% CI), min/max, count
"""

from __future__ import annotations

import os
import sys
import time
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.models import dna_client, variant_scorers
from alphagenome.models import dna_client as _dc

# ============================================================================
# CONFIGURATION & SETUP
# ============================================================================

VAR_TSV = "/cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_variants.tsv"
TSS_BED = "/cfs/klemming/scratch/m/mmarandi/experiments/dataset4/background/01_inputs/background_gene_set±10kb.bed"
GENE_LIST_PATH = "/cfs/klemming/scratch/m/mmarandi/intermediate_input/dataset5_null/background_gene_set_380.tsv"

REPO_ROOT = Path(__file__).resolve().parent.parent
PATH_MANAGER = Path("/cfs/klemming/home/m/mmarandi/lab/users/mmarandi/alphagenome/in-silico-genetic-variance/helpers/path_manager.py")
PREPARE_RESULTS_DIR = REPO_ROOT / "3_prepare_results"

if PREPARE_RESULTS_DIR.as_posix() not in sys.path:
    sys.path.insert(0, PREPARE_RESULTS_DIR.as_posix())

from modules.normalisation_helper import (  # noqa: E402
    TRACK_METADATA_COLUMNS,
    add_protocol_track_columns,
)

BATCH_SIZE = 128
NUM_ITERATIONS = 50
REST_SECONDS = 60
SEQ_LEN = 1048576  # 1 MB context window
RENAME_FRIENDLY = {
    "GeneMaskLFCScorer": "gene_exonmask_delta_log2",
}
DEDUP_KEYS = ["gene_id", "variant_id", "track_key"]

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

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def ensg_core(x: str) -> str | None:
    """Extract ENSG ID without version."""
    if x is None or pd.isna(x):
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


def print_protocol_counts(df: pd.DataFrame, label: str) -> None:
    """print counts by protocol group for quick validation.

    args:
        df (pd.DataFrame): dataframe to summarize.
        label (str): summary label.
    """
    if df.empty or "protocol_group" not in df.columns:
        print(f"  protocol counts {label}: none")
        return

    summary = (
        df.groupby("protocol_group", dropna=False)
        .agg(
            n_rows=("variant_id", "size"),
            n_track_keys=("track_key", "nunique"),
            n_variant_ids=("variant_id", "nunique"),
            n_genes=("gene_id", "nunique"),
        )
        .reset_index()
        .sort_values("protocol_group")
    )

    print(f"  protocol counts {label}:")
    for row in summary.itertuples(index=False):
        print(
            "    "
            f"{row.protocol_group}: "
            f"rows={row.n_rows:,}, "
            f"track_keys={row.n_track_keys:,}, "
            f"variant_ids={row.n_variant_ids:,}, "
            f"genes={row.n_genes:,}"
        )


def build_protocol_qc_frame(
    df: pd.DataFrame,
    iteration: int,
    stage: str,
) -> pd.DataFrame:
    """summarize protocol counts for one iteration and stage.

    args:
        df (pd.DataFrame): dataframe to summarize.
        iteration (int): one-based iteration index.
        stage (str): stage label.

    returns:
        pd.DataFrame: one row per protocol group.
    """
    if df.empty:
        return pd.DataFrame(
            columns=[
                "iteration",
                "stage",
                "protocol_group",
                "n_rows",
                "n_track_keys",
                "n_variant_ids",
                "n_genes",
            ]
        )

    summary = (
        df.groupby("protocol_group", dropna=False)
        .agg(
            n_rows=("variant_id", "size"),
            n_track_keys=("track_key", "nunique"),
            n_variant_ids=("variant_id", "nunique"),
            n_genes=("gene_id", "nunique"),
        )
        .reset_index()
    )
    summary.insert(0, "stage", stage)
    summary.insert(0, "iteration", iteration)
    return summary

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


def _first_present_value(series: pd.Series) -> object:
    non_null = series.dropna()
    if non_null.empty:
        return pd.NA
    return non_null.iloc[0]

def deduplicate_with_ci(
    df: pd.DataFrame,
    iteration: int,
    verbose: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Deduplicate by gene, variant, and track using pure Pandas.

    Keep the row with the largest absolute score within each exact track.
    Return both the deduplicated frame and a CI frame with all rows.
    """
    if verbose:
        print(f"  dedup input: {df.shape}")

    if 'raw_score' not in df.columns:
        if 'score' in df.columns:
            df = df.rename(columns={'score': 'raw_score'})
        else:
            raise ValueError('missing raw_score/score column')

    df = add_protocol_track_columns(df)
    df = df.dropna(subset=['raw_score']).copy()

    if df['track_key'].isna().any() or (df['track_key'] == '').any():
        raise ValueError('track_key missing in stability scoring output')

    if verbose:
        print_protocol_counts(df, 'before dedup')

    df['_abs_score'] = df['raw_score'].abs()
    df['_row_idx'] = np.arange(len(df))
    df = df.sort_values(
        ['_abs_score', 'gene_id', 'variant_id', 'track_key'],
        ascending=[False, True, True, True],
        kind='mergesort',
    )

    dedup = df.drop_duplicates(subset=DEDUP_KEYS, keep='first').copy()
    winner_idx = set(dedup['_row_idx'].tolist())

    ci_data = df.drop(columns=['_abs_score']).copy()
    ci_data['iteration'] = iteration
    ci_data['is_iteration_winner'] = ci_data['_row_idx'].isin(winner_idx)

    dedup = (
        dedup.drop(columns=['_abs_score', '_row_idx'])
        .sort_values(['gene_id', 'protocol_group', 'variant_id', 'track_key'])
    )
    ci_data = ci_data.drop(columns=['_row_idx']).sort_values(
        ['iteration', 'gene_id', 'protocol_group', 'variant_id', 'track_key']
    )

    if verbose:
        print(f"  deduplicated output: {dedup.shape}")
        print(f"  raw CI rows (all scores): {ci_data.shape}")
        print_protocol_counts(dedup, 'after dedup')

    return dedup, ci_data

# ============================================================================
# MAIN EXECUTION
# ============================================================================

if __name__ == "__main__":
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
    protocol_qc_frames = []

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
        tidy = add_protocol_track_columns(tidy)
        tidy["scorer_friendly"] = tidy["variant_scorer"].map(map_friendly)
        protocol_qc_frames.append(build_protocol_qc_frame(tidy, iteration, "before_dedup"))
        
        # Deduplicate with CI tracking (pure Pandas)
        dedup_df, ci_df = deduplicate_with_ci(tidy, iteration=iteration, verbose=True)
        
        # Filter BOTH dataframes by whitelist
        if gene_whitelist:
            dedup_df["gene_id_core"] = dedup_df["gene_id"].map(ensg_core)
            dedup_df = dedup_df[dedup_df["gene_id_core"].isin(gene_whitelist)].copy()
            dedup_df = dedup_df.drop(columns=["gene_id_core"], errors="ignore")

            ci_df["gene_id_core"] = ci_df["gene_id"].map(ensg_core)
            ci_df = ci_df[ci_df["gene_id_core"].isin(gene_whitelist)].copy()
            ci_df = ci_df.drop(columns=["gene_id_core"], errors="ignore")

        print_protocol_counts(dedup_df, "after whitelist")
        protocol_qc_frames.append(build_protocol_qc_frame(dedup_df, iteration, "after_dedup"))

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
    print_protocol_counts(all_ci, "across all raw CI rows")

    protocol_qc = pd.concat(protocol_qc_frames, ignore_index=True)
    protocol_qc_out = STABILITY_DIR / "protocol_qc_by_iteration.tsv.gz"
    protocol_qc.to_csv(protocol_qc_out, sep="\t", index=False, compression="gzip")
    print(f"Saved protocol QC summary to {protocol_qc_out.name}")

    # Group by variant, gene, protocol, and track to preserve exact RNA identity
    metadata_cols = [col for col in TRACK_METADATA_COLUMNS if col in all_ci.columns]
    agg_stats = (
        all_ci
        .groupby(["variant_id", "gene_id", "protocol_group", "track_key"], dropna=False)
        .agg(
            count=("raw_score", "count"),
            mean_score=("raw_score", "mean"),
            std_score=("raw_score", "std"),
            min_score=("raw_score", "min"),
            max_score=("raw_score", "max"),
            ci95_lower=("raw_score", lambda x: x.quantile(0.025)),
            ci50_lower=("raw_score", lambda x: x.quantile(0.25)),
            median_score=("raw_score", "median"),
            ci50_upper=("raw_score", lambda x: x.quantile(0.75)),
            ci95_upper=("raw_score", lambda x: x.quantile(0.975)),
            presence_iterations=("iteration", "nunique"),
            win_iterations=("is_iteration_winner", "sum"),
            **{
                col: (col, _first_present_value)
                for col in metadata_cols
            },
        )
        .reset_index()
    )
    agg_stats["presence_fraction"] = agg_stats["presence_iterations"] / NUM_ITERATIONS
    agg_stats["win_fraction"] = agg_stats["win_iterations"] / NUM_ITERATIONS

    protocol_selection_summary = (
        agg_stats
        .groupby("protocol_group", dropna=False)
        .agg(
            n_variant_gene_tracks=("track_key", "size"),
            total_presence_iterations=("presence_iterations", "sum"),
            total_win_iterations=("win_iterations", "sum"),
            mean_presence_fraction=("presence_fraction", "mean"),
            mean_win_fraction=("win_fraction", "mean"),
        )
        .reset_index()
        .sort_values("protocol_group")
    )

    # Save aggregated statistics
    agg_out = STABILITY_DIR / "aggregated_statistics.tsv.gz"
    agg_stats.to_csv(agg_out, sep="\t", index=False, compression="gzip")
    print(f"Saved aggregated statistics to {agg_out.name}")
    print(f"Records: {len(agg_stats)}")

    protocol_selection_out = STABILITY_DIR / "protocol_selection_summary.tsv.gz"
    protocol_selection_summary.to_csv(
        protocol_selection_out,
        sep="\t",
        index=False,
        compression="gzip",
    )
    print(f"Saved protocol selection summary to {protocol_selection_out.name}")

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

    print("\nMean score across all variants:")
    print(f"  Mean:   {agg_stats['mean_score'].mean():.6f}")
    print(f"  Std:    {agg_stats['mean_score'].std():.6f}")
    print(f"  Range:  [{agg_stats['mean_score'].min():.6f}, {agg_stats['mean_score'].max():.6f}]")

    print("\nScore variability (std within variant):")
    print(f"  Mean std: {agg_stats['std_score'].mean():.6f}")
    print(f"  Median:   {agg_stats['std_score'].median():.6f}")
    print(f"  Max:      {agg_stats['std_score'].max():.6f}")

    print("\nOutput files:")
    print(f"  Iterations: {ITERATIONS_DIR}")
    print(f"  Aggregated: {agg_out}")
    print(f"  Raw CI data:   {raw_ci_out}")
    print(f"  Protocol QC:   {protocol_qc_out}")
    print(f"  Protocol wins: {protocol_selection_out}")