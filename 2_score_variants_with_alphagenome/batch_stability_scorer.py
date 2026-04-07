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

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from alphagenome_shared import build_layout, resolve_legacy_input_paths

def _env_int(name: str, default: int) -> int:
    """Parse an integer environment variable with fallback."""
    raw = os.getenv(name)
    if raw is None:
        return default
    try:
        return int(raw)
    except ValueError as exc:
        raise ValueError(f"{name} must be an integer, got {raw!r}") from exc

def _env_flag(name: str, default: bool = False) -> bool:
    """Parse a boolean environment variable with fallback."""
    raw = os.getenv(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "y", "on"}

BATCH_SIZE = _env_int("AG_BATCH_STABILITY_BATCH_SIZE", 128)
NUM_ITERATIONS = _env_int("AG_BATCH_STABILITY_NUM_ITERATIONS", 50)
REST_SECONDS = _env_int("AG_BATCH_STABILITY_REST_SECONDS", 60)
SEQ_LEN = _env_int("AG_BATCH_STABILITY_SEQ_LEN", 1048576)  # 1 MB context window
QUERY_START = _env_int("AG_BATCH_STABILITY_START_INDEX", 0)
RESET_OUTPUT = _env_flag("AG_BATCH_STABILITY_RESET", default=False)
RENAME_FRIENDLY = {
    "GeneMaskLFCScorer": "gene_exonmask_delta_log2",
}

os.environ.setdefault("DATASET_ID", "dataset_4")
os.environ.setdefault("SAMPLE_ID", "background_batch_stability")

layout = build_layout()
layout.make_dirs()
VAR_TSV, _, GENE_LIST_PATH = resolve_legacy_input_paths(layout)

STABILITY_DIR = layout.results_dir / "batch_stability"
if RESET_OUTPUT and STABILITY_DIR.exists():
    archived_dir = STABILITY_DIR.parent / f"{STABILITY_DIR.name}.reset_{datetime.now():%Y%m%d_%H%M%S}"
    shutil.move(STABILITY_DIR.as_posix(), archived_dir.as_posix())
    print(f"Archived existing batch stability outputs to {archived_dir}")

STABILITY_DIR.mkdir(exist_ok=True, parents=True)
ITERATIONS_DIR = STABILITY_DIR / "iterations"
ITERATIONS_DIR.mkdir(exist_ok=True, parents=True)
PROTOCOL_ITERATIONS_DIR = STABILITY_DIR / "iterations_by_protocol"
PROTOCOL_ITERATIONS_DIR.mkdir(exist_ok=True, parents=True)
RAW_CI_DIR = STABILITY_DIR / "raw_ci"
RAW_CI_DIR.mkdir(exist_ok=True, parents=True)

print(f"Output directory: {STABILITY_DIR}")
print(f"Iterations directory: {ITERATIONS_DIR}")
print(f"Protocol iterations directory: {PROTOCOL_ITERATIONS_DIR}")
print(f"Raw CI directory: {RAW_CI_DIR}")

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


def _add_protocol_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Annotate rows with protocol_group and stable track identity."""
    df = df.copy()

    assay_col = "Assay_title" if "Assay_title" in df.columns else (
        "Assay title" if "Assay title" in df.columns else None
    )
    assay_source = df[assay_col] if assay_col else pd.Series(pd.NA, index=df.index)
    track_source = df["track_name"] if "track_name" in df.columns else pd.Series(pd.NA, index=df.index)
    strand_source = df["track_strand"] if "track_strand" in df.columns else pd.Series(pd.NA, index=df.index)
    ontology_source = df["ontology_curie"] if "ontology_curie" in df.columns else pd.Series(pd.NA, index=df.index)

    protocol_source = assay_source.fillna(track_source).fillna(ontology_source)
    protocol_norm = protocol_source.map(_normalize_label)

    df["protocol_group"] = np.where(
        protocol_norm.str.contains("polya") & protocol_norm.str.contains("rna_seq"),
        "polyA_plus_rna_seq",
        np.where(
            protocol_norm.str.contains("total") & protocol_norm.str.contains("rna_seq"),
            "total_rna_seq",
            "other",
        ),
    )

    df["assay_label"] = assay_source.astype("string").fillna("")
    df["track_label"] = track_source.astype("string").fillna("")
    df["track_key"] = (
        df["track_label"].replace("", pd.NA)
        .fillna(df["assay_label"].replace("", pd.NA))
        .fillna(ontology_source.astype("string").replace("", pd.NA))
        .fillna("unknown_track")
    )
    df["track_key"] = (
        df["protocol_group"].astype(str)
        + "||"
        + df["track_key"].astype(str)
        + "||"
        + strand_source.astype("string").fillna("")
    )
    return df

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
        tidy = tidy.drop(columns=["variant"])
    if "scored_interval" in tidy.columns:
        tidy["scored_interval_str"] = tidy["scored_interval"].map(_interval_to_str)
        tidy = tidy.drop(columns=["scored_interval"])
    else:
        tidy["scored_interval_str"] = pd.NA
    return tidy

def _ensure_variant_ids(df: pd.DataFrame, interval_to_varid: dict[str, str]) -> None:
    """Fill missing variant IDs using the scored interval lookup."""
    if "variant_id" not in df.columns:
        df["variant_id"] = pd.NA
    mask = df["variant_id"].isna() | (df["variant_id"].astype(str).str.strip() == "")
    if mask.any():
        df.loc[mask, "variant_id"] = df.loc[mask, "scored_interval_str"].map(interval_to_varid)

def _iteration_result_path(iteration: int) -> Path:
    """Path for per-iteration deduplicated results."""
    return ITERATIONS_DIR / f"iteration_{iteration:02d}.tsv.gz"

def _iteration_ci_path(iteration: int) -> Path:
    """Path for per-iteration raw CI rows."""
    return RAW_CI_DIR / f"iteration_{iteration:02d}.raw_ci.tsv.gz"

def _iteration_protocol_path(iteration: int) -> Path:
    """Path for per-iteration protocol-separated selected rows."""
    return PROTOCOL_ITERATIONS_DIR / f"iteration_{iteration:02d}.protocol.tsv.gz"

def _load_iteration_score_frames(
    max_iterations: int,
) -> tuple[pd.DataFrame, list[int], list[int]]:
    """Load usable iteration result rows from disk for aggregation."""
    frames: list[pd.DataFrame] = []
    completed_iterations: list[int] = []
    invalid_iterations: list[int] = []

    for iteration in range(1, max_iterations + 1):
        iter_path = _iteration_result_path(iteration)
        if not iter_path.exists():
            continue

        completed_iterations.append(iteration)
        iter_df = pd.read_csv(iter_path, sep="\t", compression="gzip", low_memory=False)
        required = {"variant_id", "gene_id", "raw_score"}
        missing = required - set(iter_df.columns)
        if missing:
            invalid_iterations.append(iteration)
            print(
                f"Warning: iteration {iteration} is missing required columns {sorted(missing)}; "
                "excluding it from aggregation"
            )
            continue

        iter_df["variant_id"] = iter_df["variant_id"].astype("string").str.strip().replace("", pd.NA)
        iter_df["gene_id"] = iter_df["gene_id"].astype("string").str.strip().replace("", pd.NA)
        iter_df["raw_score"] = pd.to_numeric(iter_df["raw_score"], errors="coerce")
        usable = iter_df.dropna(subset=["variant_id", "gene_id", "raw_score"]).copy()
        if usable.empty:
            invalid_iterations.append(iteration)
            print(
                f"Warning: iteration {iteration} has no usable rows with variant_id/gene_id/raw_score; "
                "excluding it from aggregation"
            )
            continue

        usable["iteration"] = iteration
        frames.append(usable[["iteration", "variant_id", "gene_id", "raw_score"]])

    if not frames:
        empty = pd.DataFrame(columns=["iteration", "variant_id", "gene_id", "raw_score"])
        return empty, completed_iterations, invalid_iterations

    return pd.concat(frames, ignore_index=True), completed_iterations, invalid_iterations


def _load_protocol_iteration_frames(
    max_iterations: int,
) -> tuple[pd.DataFrame, list[int], list[int]]:
    """Load usable protocol-separated iteration rows from disk for aggregation."""
    frames: list[pd.DataFrame] = []
    completed_iterations: list[int] = []
    invalid_iterations: list[int] = []

    for iteration in range(1, max_iterations + 1):
        iter_path = _iteration_protocol_path(iteration)
        if not iter_path.exists():
            continue

        completed_iterations.append(iteration)
        iter_df = pd.read_csv(iter_path, sep="\t", compression="gzip", low_memory=False)
        required = {"variant_id", "gene_id", "protocol_group", "raw_score"}
        missing = required - set(iter_df.columns)
        if missing:
            invalid_iterations.append(iteration)
            print(
                f"Warning: protocol iteration {iteration} is missing required columns {sorted(missing)}; "
                "excluding it from protocol aggregation"
            )
            continue

        iter_df["variant_id"] = iter_df["variant_id"].astype("string").str.strip().replace("", pd.NA)
        iter_df["gene_id"] = iter_df["gene_id"].astype("string").str.strip().replace("", pd.NA)
        iter_df["protocol_group"] = iter_df["protocol_group"].astype("string").str.strip().replace("", pd.NA)
        iter_df["raw_score"] = pd.to_numeric(iter_df["raw_score"], errors="coerce")
        usable = iter_df.dropna(subset=["variant_id", "gene_id", "protocol_group", "raw_score"]).copy()
        if usable.empty:
            invalid_iterations.append(iteration)
            print(
                f"Warning: protocol iteration {iteration} has no usable rows; excluding it from protocol aggregation"
            )
            continue

        usable["iteration"] = iteration
        frames.append(usable)

    if not frames:
        empty = pd.DataFrame(columns=["iteration", "variant_id", "gene_id", "protocol_group", "raw_score"])
        return empty, completed_iterations, invalid_iterations

    return pd.concat(frames, ignore_index=True), completed_iterations, invalid_iterations

def _load_raw_ci_frames(max_iterations: int) -> tuple[pd.DataFrame, list[int], list[int]]:
    """Load available per-iteration raw CI rows from disk."""
    frames: list[pd.DataFrame] = []
    loaded_iterations: list[int] = []
    missing_iterations: list[int] = []

    for iteration in range(1, max_iterations + 1):
        iter_path = _iteration_result_path(iteration)
        ci_path = _iteration_ci_path(iteration)
        if not iter_path.exists():
            continue
        if not ci_path.exists():
            missing_iterations.append(iteration)
            continue

        ci_df = pd.read_csv(ci_path, sep="\t", compression="gzip", low_memory=False)
        if "iteration" not in ci_df.columns:
            ci_df["iteration"] = iteration
        frames.append(ci_df)
        loaded_iterations.append(iteration)

    if not frames:
        empty = pd.DataFrame(columns=["gene_id", "variant_id", "raw_score", "iteration"])
        return empty, loaded_iterations, missing_iterations

    return pd.concat(frames, ignore_index=True), loaded_iterations, missing_iterations

def _existing_iteration_numbers(max_iterations: int) -> list[int]:
    """List completed iteration numbers already on disk."""
    return [iteration for iteration in range(1, max_iterations + 1) if _iteration_result_path(iteration).exists()]

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
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build three iteration outputs:
    - global selected rows: max abs score across all tracks
    - protocol-selected rows: max abs score within each protocol only
    - raw track rows with protocol metadata for QC / CI analysis

    Returns:
        (dedup_df, protocol_df, ci_df)
    """
    if verbose:
        print(f"  dedup input: {df.shape}")

    df_pd = df.to_pandas() if isinstance(df, pl.DataFrame) else df.copy()
    schema_cols = df_pd.columns
    if 'raw_score' in schema_cols:
        score_col = 'raw_score'
    elif 'score' in schema_cols:
        score_col = 'score'
    else:
        raise ValueError('missing raw_score/score column')

    if score_col != 'raw_score':
        df_pd = df_pd.rename(columns={score_col: 'raw_score'})

    df_pd["raw_score"] = pd.to_numeric(df_pd["raw_score"], errors="coerce")
    df_pd = df_pd[df_pd["raw_score"].notna()].copy()
    df_pd = _add_protocol_columns(df_pd)

    df_pl = pl.from_pandas(df_pd)
    df_pl = df_pl.with_columns(pl.col('raw_score').abs().alias('_abs_score'))
    df_pl = df_pl.sort('_abs_score', descending=True)

    dedup = (
        df_pl
        .unique(subset=['gene_id', 'variant_id'], keep='first', maintain_order=True)
        .drop('_abs_score')
        .sort(['gene_id', 'variant_id'])
        .to_pandas()
    )

    protocol_counts = (
        df_pl
        .group_by(['gene_id', 'variant_id', 'protocol_group'])
        .agg(
            pl.len().alias('n_candidate_tracks_in_protocol'),
            pl.col('track_key').n_unique().alias('n_distinct_tracks_in_protocol'),
        )
    )
    protocol_dedup = (
        df_pl
        .unique(subset=['gene_id', 'variant_id', 'protocol_group'], keep='first', maintain_order=True)
        .join(protocol_counts, on=['gene_id', 'variant_id', 'protocol_group'], how='left')
        .drop('_abs_score')
        .sort(['gene_id', 'variant_id', 'protocol_group'])
        .to_pandas()
    )

    ci_data = df_pl.drop('_abs_score').to_pandas()
    ci_data['iteration'] = iteration

    if verbose:
        print(f"  globally selected output: {dedup.shape}")
        print(f"  protocol-selected output: {protocol_dedup.shape}")
        print(f"  raw track rows: {ci_data.shape}")

    return dedup, protocol_dedup, ci_data


def _aggregate_score_statistics(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    """Aggregate repeated score observations into summary statistics."""
    agg_stats = (
        df
        .groupby(group_cols)
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
    agg_stats["std_score"] = agg_stats["std_score"].fillna(0.0)
    return agg_stats

# ============================================================================
# MAIN EXECUTION
# ============================================================================

print("=" * 80)
print("BATCH STABILITY SCORER: repeated same-variant queries with CI estimation")
print("=" * 80)
print(
    f"Config: batch_size={BATCH_SIZE}, iterations={NUM_ITERATIONS}, "
    f"rest_seconds={REST_SECONDS}, seq_len={SEQ_LEN}, "
    f"start_index={QUERY_START}"
)

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

# Take only the configured slice of variants
if QUERY_START < 0:
    raise ValueError(f"AG_BATCH_STABILITY_START_INDEX must be >= 0, got {QUERY_START}")

query_batch = df.iloc[QUERY_START:QUERY_START + BATCH_SIZE].copy()
if query_batch.empty:
    raise ValueError(
        f"No variants selected from start index {QUERY_START} with batch size {BATCH_SIZE}"
    )
if len(query_batch) < BATCH_SIZE:
    print(
        f"Warning: requested {BATCH_SIZE} variants starting at index {QUERY_START}, "
        f"but only found {len(query_batch)}"
    )
print(
    f"Selected variants [{QUERY_START}:{QUERY_START + len(query_batch)}) "
    f"for repeated queries"
)

existing_iterations = _existing_iteration_numbers(NUM_ITERATIONS)
if existing_iterations:
    print(
        f"\nFound {len(existing_iterations)} existing iteration result files: "
        f"{', '.join(f'{it:02d}' for it in existing_iterations[:10])}"
        + (" ..." if len(existing_iterations) > 10 else "")
    )
    _, _, invalid_existing_iterations = _load_iteration_score_frames(NUM_ITERATIONS)
    if invalid_existing_iterations:
        preview = ", ".join(f"{it:02d}" for it in invalid_existing_iterations[:10])
        more = "" if len(invalid_existing_iterations) <= 10 else f" ... (+{len(invalid_existing_iterations) - 10} more)"
        raise RuntimeError(
            "Found invalid existing iteration outputs that should not be resumed: "
            f"{preview}{more}. Start fresh with AG_BATCH_STABILITY_RESET=1 or use a new SAMPLE_ID."
        )

pending_iterations = [iteration for iteration in range(1, NUM_ITERATIONS + 1) if iteration not in existing_iterations]
if pending_iterations:
    print(
        f"Pending iterations to score: {len(pending_iterations)} "
        f"(starting at {pending_iterations[0]:02d})"
    )
else:
    print("All requested iteration files already exist; skipping new scoring and rebuilding summaries from disk.")

RNA = dna_client.OutputType.RNA_SEQ
ORG = _dc.Organism.HOMO_SAPIENS

if pending_iterations:
    # Initialize API only when there is scoring left to do.
    RAW_CI_DIR.mkdir(exist_ok=True, parents=True)
    print(f"\nInitializing AlphaGenome client...")
    t0_client = time.time()
    model = dna_client.create(api_key=API_KEY)
    print(f"API client ready in {time.time() - t0_client:.1f}s")

    # Load gene whitelist only when we are about to score.
    gene_whitelist = _load_gene_whitelist(GENE_LIST_PATH)

    for pending_idx, iteration in enumerate(pending_iterations, start=1):
        print(f"\n{'=' * 80}")
        print(f"ITERATION {iteration}/{NUM_ITERATIONS}")
        print(f"{'=' * 80}")

        start_iter = time.time()

        # Prepare batch and a lookup that can recover missing variant IDs.
        meta = query_batch.copy()
        intervals = make_intervals(SEQ_LEN, meta["CHROM"].to_numpy(), meta["POS"].to_numpy())
        interval_to_varid = {_interval_to_str(iv): vid for iv, vid in zip(intervals, meta["variant_id"])}
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

        # Convert to dataframe and force stable variant IDs before deduplication.
        inject_into_anndata(scores, meta)
        tidy = variant_scorers.tidy_scores(scores, match_gene_strand=False)
        tidy = normalize_tidy(tidy)
        _ensure_variant_ids(tidy, interval_to_varid)
        missing_variant_ids = tidy["variant_id"].isna() | (tidy["variant_id"].astype(str).str.strip() == "")
        if missing_variant_ids.any():
            raise RuntimeError(
                f"Iteration {iteration} still has {int(missing_variant_ids.sum())} rows without variant_id "
                "after interval-based backfill"
            )
        tidy["scorer_friendly"] = tidy["variant_scorer"].map(map_friendly)

        # Convert to polars for deduplication
        tidy_pl = pl.from_pandas(tidy)

        # Build global and protocol-separated outputs with full track-level rows.
        dedup_df, protocol_df, ci_df = deduplicate_with_ci(tidy_pl, iteration=iteration, verbose=True)

        # Filter by whitelist
        if gene_whitelist:
            dedup_df["gene_id_core"] = dedup_df["gene_id"].map(ensg_core)
            dedup_df = dedup_df[dedup_df["gene_id_core"].isin(gene_whitelist)].copy()
            dedup_df = dedup_df.drop(columns=["gene_id_core"], errors="ignore")

            protocol_df["gene_id_core"] = protocol_df["gene_id"].map(ensg_core)
            protocol_df = protocol_df[protocol_df["gene_id_core"].isin(gene_whitelist)].copy()
            protocol_df = protocol_df.drop(columns=["gene_id_core"], errors="ignore")

            ci_df["gene_id_core"] = ci_df["gene_id"].map(ensg_core)
            ci_df = ci_df[ci_df["gene_id_core"].isin(gene_whitelist)].copy()
            ci_df = ci_df.drop(columns=["gene_id_core"], errors="ignore")

        iter_out = _iteration_result_path(iteration)
        protocol_out = _iteration_protocol_path(iteration)
        ci_out = _iteration_ci_path(iteration)

        # Save global selected rows, protocol-separated selected rows, and full raw track rows.
        dedup_df.to_csv(iter_out, sep="\t", index=False, compression="gzip")
        protocol_df.to_csv(protocol_out, sep="\t", index=False, compression="gzip")
        ci_df.to_csv(ci_out, sep="\t", index=False, compression="gzip")
        print(f"Saved iteration {iteration} results to {iter_out.name}")
        print(f"Saved iteration {iteration} protocol-separated results to {protocol_out.name}")
        print(f"Saved iteration {iteration} raw CI rows to {ci_out.name}")

        if pending_idx < len(pending_iterations):
            print(f"Resting for {REST_SECONDS} seconds before next iteration...")
            time.sleep(REST_SECONDS)

        elapsed_iter = time.time() - start_iter
        print(f"Iteration {iteration} completed in {elapsed_iter:.1f}s")

print(f"\n{'=' * 80}")
print("SCORING PHASE COMPLETE")
print(f"{'=' * 80}")

iteration_scores, completed_iterations, invalid_iterations = _load_iteration_score_frames(NUM_ITERATIONS)
if invalid_iterations:
    preview = ", ".join(f"{it:02d}" for it in invalid_iterations[:10])
    more = "" if len(invalid_iterations) <= 10 else f" ... (+{len(invalid_iterations) - 10} more)"
    raise RuntimeError(
        "Cannot aggregate because some iteration outputs are invalid: "
        f"{preview}{more}. Restart cleanly with AG_BATCH_STABILITY_RESET=1 or use a new SAMPLE_ID."
    )
if iteration_scores.empty:
    raise RuntimeError("No usable iteration outputs were found for aggregation.")

print(f"\nAggregating results across {len(completed_iterations)} completed iterations...")

agg_stats = _aggregate_score_statistics(iteration_scores, ['variant_id', 'gene_id'])

agg_out = STABILITY_DIR / "aggregated_statistics.tsv.gz"
agg_stats.to_csv(agg_out, sep="\t", index=False, compression="gzip")
print(f"Saved aggregated statistics to {agg_out.name}")
print(f"Records: {len(agg_stats)}")

protocol_iteration_scores, protocol_completed_iterations, invalid_protocol_iterations = _load_protocol_iteration_frames(NUM_ITERATIONS)
if invalid_protocol_iterations:
    preview = ", ".join(f"{it:02d}" for it in invalid_protocol_iterations[:10])
    more = "" if len(invalid_protocol_iterations) <= 10 else f" ... (+{len(invalid_protocol_iterations) - 10} more)"
    raise RuntimeError(
        "Cannot aggregate protocol-separated outputs because some iteration protocol files are invalid: "
        f"{preview}{more}."
    )

protocol_scores_out = STABILITY_DIR / "protocol_iteration_scores.tsv.gz"
protocol_agg_out = STABILITY_DIR / "aggregated_statistics_by_protocol.tsv.gz"
if protocol_iteration_scores.empty:
    print("No protocol-separated iteration outputs found; skipping protocol aggregation")
else:
    protocol_iteration_scores.to_csv(protocol_scores_out, sep="\t", index=False, compression="gzip")

    protocol_agg = _aggregate_score_statistics(
        protocol_iteration_scores,
        ['variant_id', 'gene_id', 'protocol_group'],
    )

    if 'track_key' in protocol_iteration_scores.columns:
        ordered_protocol = protocol_iteration_scores.sort_values(
            ['variant_id', 'gene_id', 'protocol_group', 'iteration']
        ).copy()
        track_summary = (
            ordered_protocol
            .groupby(['variant_id', 'gene_id', 'protocol_group'])
            .agg(
                n_distinct_selected_tracks=('track_key', 'nunique'),
                track_switches=('track_key', lambda s: int((s != s.shift()).sum() - 1)),
                dominant_track_key=('track_key', lambda s: s.mode().iloc[0] if not s.mode().empty else s.iloc[0]),
                dominant_track_name=('track_name', lambda s: s.mode().iloc[0] if not s.mode().empty else s.iloc[0]),
            )
            .reset_index()
        )
        protocol_agg = protocol_agg.merge(
            track_summary,
            on=['variant_id', 'gene_id', 'protocol_group'],
            how='left',
        )

    protocol_agg.to_csv(protocol_agg_out, sep="\t", index=False, compression="gzip")
    print(f"Saved protocol iteration scores to {protocol_scores_out.name}")
    print(f"Saved protocol-separated aggregated statistics to {protocol_agg_out.name}")
    print(f"Protocol records: {len(protocol_agg)}")

all_ci, ci_iterations, missing_ci_iterations = _load_raw_ci_frames(NUM_ITERATIONS)
raw_ci_out = STABILITY_DIR / "all_raw_ci_data.tsv.gz"
if all_ci.empty:
    print("No per-iteration raw CI files were found; skipping combined raw CI export")
else:
    all_ci.to_csv(raw_ci_out, sep="\t", index=False, compression="gzip")
    print(f"Saved raw CI data to {raw_ci_out.name}")
if missing_ci_iterations:
    print(
        f"Warning: raw CI companion files are missing for {len(missing_ci_iterations)} iteration(s): "
        f"{', '.join(f'{it:02d}' for it in missing_ci_iterations[:10])}"
        + (" ..." if len(missing_ci_iterations) > 10 else "")
    )

print(f"Total unique variant/gene pairs scored: {len(agg_stats)}")
print(f"Requested iterations: {NUM_ITERATIONS}")
print(f"Completed iterations on disk: {len(completed_iterations)}")
print(f"Total raw CI rows collected: {len(all_ci)}")
print(f"Mean observations per variant/gene pair: {agg_stats['count'].mean():.1f}")
print(f"Min observations per variant/gene pair: {agg_stats['count'].min()}")
print(f"Max observations per variant/gene pair: {agg_stats['count'].max()}")

print(f"\nMean score across all variant/gene pairs:")
print(f"  Mean:   {agg_stats['mean_score'].mean():.6f}")
print(f"  Std:    {agg_stats['mean_score'].std():.6f}")
print(f"  Range:  [{agg_stats['mean_score'].min():.6f}, {agg_stats['mean_score'].max():.6f}]")

print(f"\nScore variability (std within variant/gene pair):")
print(f"  Mean std: {agg_stats['std_score'].mean():.6f}")
print(f"  Median:   {agg_stats['std_score'].median():.6f}")
print(f"  Max:      {agg_stats['std_score'].max():.6f}")

print("Output files:")
print(f"  Iterations: {ITERATIONS_DIR}")
print(f"  Protocol iterations: {PROTOCOL_ITERATIONS_DIR}")
print(f"  Raw CI dir: {RAW_CI_DIR}")
print(f"  Aggregated: {agg_out}")
if not protocol_iteration_scores.empty:
    print(f"  Protocol iteration scores: {protocol_scores_out}")
    print(f"  Protocol aggregated: {protocol_agg_out}")
if not all_ci.empty:
    print(f"  Raw CI data: {raw_ci_out}")
