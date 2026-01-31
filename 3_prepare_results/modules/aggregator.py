from __future__ import annotations

from pathlib import Path
import numpy as np
import polars as pl

from .external_data_loader import ExternalDataLoader
from .normalisation_helper import strip_ensembl_version
from .permuted_af import load_gene_af_pools


def _load_gene_meta(base: Path) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """load and normalize all reference metadata (returns eager dataframes)."""
    loader = ExternalDataLoader(base)
    mane = loader.load_mane()
    gtf = loader.load_gtf_genes()
    tpm = loader.load_tpm()
    vgh = loader.load_vgh()

    def _norm(df: pl.DataFrame) -> pl.DataFrame:
        if "gene_id" in df.columns:
            return df.with_columns(pl.col("gene_id").str.split(".").list.get(0))
        return df

    return _norm(mane), _norm(gtf), _norm(tpm), _norm(vgh)

def calculate_vg_ci_struct(
        scores: pl.Series,
        gene_id: str,
        af_pools: dict[str, np.ndarray],
        n_iter: int = 1000
) -> dict[str, float]:
    """Monte Carlo simulation for Vg CI."""
    if len(scores) == 0:
        return {"vg_perm_mean": 0.0, "vg_perm_p05": 0.0, "vg_perm_p95": 0.0}

    pool = af_pools.get(gene_id)
    if pool is None or len(pool) == 0:
        return {"vg_perm_mean": 0.0, "vg_perm_p05": 0.0, "vg_perm_p95": 0.0}

    scores_arr = scores.to_numpy()
    n_variants = len(scores_arr)

    rng = np.random.default_rng()
    random_afs = rng.choice(pool, size=(n_iter, n_variants), replace=True)

    term_p = 2.0 * random_afs * (1.0 - random_afs)
    term_beta = scores_arr ** 2

    vg_sims = term_p @ term_beta

    return {
        "vg_perm_mean": float(np.mean(vg_sims)),
        "vg_perm_p05": float(np.percentile(vg_sims, 5)),
        "vg_perm_p95": float(np.percentile(vg_sims, 95)),
    }


def calc_architecture_stats(struct_list: list[dict]) -> dict:
    """
    Computes N90, N85, and properties of the 'driver' variants.
    Expects list of dicts: [{'v': vg_contribution, 'e': abs_score}, ...]
    """
    if not struct_list:
        return {
            "N90": 0, "N85": 0,
            "variance_N90": 0.0, "cv_effect_N90": 0.0,
            "mean_effect_N90": 0.0
        }

    # Extract arrays
    v = np.array([x['v'] for x in struct_list])
    e = np.array([x['e'] for x in struct_list])

    # Sort descending by Variance Contribution (v)
    sort_idx = np.argsort(v)[::-1]
    v_sorted = v[sort_idx]
    e_sorted = e[sort_idx]

    total_vg = np.sum(v_sorted)
    if total_vg == 0:
        return {
            "N90": 0, "N85": 0,
            "variance_N90": 0.0, "cv_effect_N90": 0.0,
            "mean_effect_N90": 0.0
        }

    cumsum = np.cumsum(v_sorted)

    # --- Logic for N90 (90% of variance) ---
    idx_90 = np.searchsorted(cumsum, 0.90 * total_vg)

    # The subset of variants driving this 90%
    subset_e_90 = e_sorted[:idx_90 + 1]

    # Stats for N90 subset
    n90_count = idx_90 + 1
    mean_90 = np.mean(subset_e_90) if len(subset_e_90) > 0 else 0.0
    std_90 = np.std(subset_e_90) if len(subset_e_90) > 0 else 0.0
    cv_90 = (std_90 / mean_90) if mean_90 > 0 else 0.0

    # --- Logic for N85 ---
    idx_85 = np.searchsorted(cumsum, 0.85 * total_vg)
    n85_count = idx_85 + 1

    return {
        "N90": int(n90_count),
        "N85": int(n85_count),
        "variance_N90": float(cumsum[idx_90]),  # Actual Vg sum of the top 90%
        "cv_effect_N90": float(cv_90),
        "mean_effect_N90": float(mean_90)
    }


# --- 4. Spatial Window Expression Generator ---
def get_window_exprs(windows: dict[str, tuple[int, int]]) -> list[pl.Expr]:
    """Generates agg expressions for N, Mean Abs, and Sum Vg for spatial windows."""
    exprs = []
    for name, (start, end) in windows.items():
        # Condition: start <= dist_signed < end
        cond = (pl.col("dist_signed") >= start) & (pl.col("dist_signed") < end)

        exprs.extend([
            pl.col("abs_score").filter(cond).mean().alias(f"mean_abs_{name}"),
            pl.col("vg_contribution").filter(cond).sum().alias(f"vg_{name}"),  # Localized Variance
            cond.sum().alias(f"n_variants_{name}")
        ])
    return exprs


# --- 5. Main Aggregation Function ---
def aggregate_genes(
        variants_path: Path,
        out_path: Path,
        base_ref: Path,
        is_ism: bool,
        gene_list_path: Path | None = None,
        calculate_ci: bool = False,
        real_reference_path: Path | None = None,
        n_permutations: int = 1000,
) -> None:
    """aggregate variant parquet to gene metrics."""

    # 1. Load metadata
    mane, gtf, tpm, vgh = _load_gene_meta(base_ref)

    # 2. Scan variants
    lf = pl.scan_parquet(variants_path)
    schema_cols = set(lf.collect_schema().names())
    has_perm_af = "perm_AF" in schema_cols

    # Normalize 'AF' column source
    af_candidates = ["AF", "AF_x", "AF_y", "af", "af_x", "af_y"]
    af_source = next((c for c in af_candidates if c in schema_cols), None)

    # Ensure we have a clean 'AF' column
    if af_source and af_source != "AF":
        lf = lf.with_columns(pl.col(af_source).alias("AF"))
        schema_cols.add("AF")
    elif "AF" not in schema_cols:
        lf = lf.with_columns(pl.lit(0.0).alias("AF"))
        schema_cols.add("AF")

    # Clean aliases
    extras = [c for c in af_candidates if c != "AF" and c in schema_cols]
    if extras: lf = lf.drop(extras)

    # 3. Gene ID Normalization
    if "gene_norm" in schema_cols:
        gene_col = "gene_norm"
    elif "gene_id" in schema_cols:
        gene_col = "gene_id"
        lf = lf.with_columns(gene_norm=pl.col("gene_id").str.split(".").list.get(0))
        gene_col = "gene_norm"
    else:
        gene_col = "gene_id"

    # 4. Whitelist Filter
    if gene_list_path and gene_list_path.exists():
        gene_whitelist = {
            strip_ensembl_version(line.strip())
            for line in gene_list_path.read_text().splitlines()
            if line.strip()
        }
        if gene_whitelist:
            lf = lf.filter(pl.col(gene_col).is_in(gene_whitelist))

    # 5. Attach Spatial Info
    gtf_spatial = gtf.lazy().select(
        [pl.col("gene_id"), "tss", "strand", "start", "end"]
    ).unique("gene_id")
    lf = lf.join(gtf_spatial, left_on=gene_col, right_on="gene_id", how="left")

    # 6. Deduplication (Max Impact)
    lf = lf.with_columns(abs_score=pl.col("raw_score").abs())
    lf = (
        lf.sort("abs_score", descending=True, nulls_last=True)
        .unique(subset=["variant_id"], keep="first", maintain_order=True)
    )

    # 7. Pre-Calculations (Distances & Vg components)
    lf = lf.with_columns(
        dist_to_tss=pl.when((pl.col("POS").is_not_null()) & (pl.col("tss").is_not_null()))
        .then((pl.col("POS") - pl.col("tss")).abs())
        .otherwise(None),

        # Signed Distance: - = Upstream, + = Downstream
        dist_signed=pl.when(
            (pl.col("POS").is_not_null()) & (pl.col("tss").is_not_null()) & (pl.col("strand").is_not_null())
        )
        .then(
            pl.when(pl.col("strand") == "+")
            .then(pl.col("POS") - pl.col("tss"))
            .otherwise(pl.col("tss") - pl.col("POS"))
        )
        .otherwise(None),

        # Real Vg
        vg_contribution=pl.when(
            pl.col("raw_score").is_not_null() & pl.col("AF").is_not_null()
        )
        .then(2.0 * pl.col("AF") * (1.0 - pl.col("AF")) * pl.col("raw_score").pow(2))
        .otherwise(0.0),

        # Permuted Vg
        vg_contribution_perm=pl.when(
            pl.col("raw_score").is_not_null() & pl.col("perm_AF").is_not_null()
        )
        .then(2.0 * pl.col("perm_AF") * (1.0 - pl.col("perm_AF")) * pl.col("raw_score").pow(2))
        .otherwise(0.0) if has_perm_af else pl.lit(0.0),
    ).fill_nan(0.0).fill_null(0.0)

    spatial_windows = {
        # 1. Distal Upstream (-10kb to -2kb)
        # Regulatory region (Enhancers)
        "distal_upstream": (-10000, -2000),

        # 2. Proximal Upstream (-2kb to -200bp)
        # Immediate regulatory environment
        "proximal_upstream": (-2000, -200),

        # 3. Core Promoter (-200bp to +200bp)
        # The physical landing pad for RNA Pol II (Crossing 0)
        "promoter_core": (-200, 200),

        # 4. Downstream Proximal (+200bp to +2kb)
        # 5' UTR and First Intron (Splicing/Translation regulation)
        "down_proximal": (200, 2000),

        # 5. Downstream Distal (+2kb to +10kb)
        # The rest of the gene body you scored
        "down_distal": (2000, 10000)
    }

    # --- 9. Build Aggregation Expressions ---
    agg_exprs = [
        pl.count().alias("n_variants"),

        pl.col("vg_contribution").sum().alias("vg_predicted"),  # From AF
        pl.col("vg_contribution_perm").sum().alias("vg_predicted_perm"),  # From perm_AF

        # Pack data for Architecture Analysis (Step 11)
        pl.struct(["vg_contribution", "abs_score"].alias("v", "e")).alias("arch_data"),

        pl.col("raw_score").pow(2).sum().alias("sum_sq_raw_score"),
        pl.col("raw_score").mean().alias("mean_raw_score"),

        # Global stats
        pl.col("abs_score").mean().alias("mean_abs_effect"),
        pl.col("abs_score").median().alias("median_abs_effect"),
        pl.col("abs_score").std().alias("std_abs_effect"),
        pl.col("abs_score").min().alias("min_abs_effect"),
        pl.col("abs_score").max().alias("max_abs_effect"),
        pl.col("abs_score").skew().alias("skewness_effect"),
        pl.col("abs_score").quantile(0.99).alias("q99_abs_effect"),

        # ID tracking
        pl.col("variant_id").sort_by("raw_score").first().alias("min_variant_id"),
        pl.col("raw_score").min().alias("min_variant_score"),
        pl.col("variant_id").sort_by("raw_score").last().alias("max_variant_id"),
        pl.col("raw_score").max().alias("max_variant_score"),

        # Distance stats
        pl.col("dist_to_tss").mean().alias("mean_dist_to_tss"),
        pl.col("dist_to_tss").median().alias("median_dist_to_tss"),
        pl.col("dist_to_tss").min().alias("min_dist_to_tss"),
        pl.col("dist_to_tss").max().alias("max_dist_to_tss"),

        # High impact counts
        (pl.col("abs_score") > 0.5).sum().alias("n_high_impact_gt05"),
        (pl.col("abs_score") > 1.0).sum().alias("n_high_impact_gt1"),
    ]

    # Add Dynamic Spatial Window Expressions
    agg_exprs.extend(get_window_exprs(spatial_windows))

    # --- 10. Perform Aggregation ---
    print("Collecting and Aggregating Genes...")
    df_agg = lf.group_by(gene_col).agg(agg_exprs).collect()

    # --- 11. Post-Aggregation: Compute Architecture (Python UDF) ---
    print("Computing Architecture Metrics (N90, CV_N90)...")

    arch_schema = pl.Struct({
        "N90": pl.Int64, "N85": pl.Int64,
        "variance_N90": pl.Float64, "cv_effect_N90": pl.Float64,
        "mean_effect_N90": pl.Float64
    })

    df_agg = df_agg.with_columns(
        pl.col("arch_data").map_elements(
            lambda x: calc_architecture_stats(x),
            return_dtype=arch_schema
        ).alias("arch_metrics")
    ).unnest("arch_metrics").drop("arch_data")

    # --- 12. Monte Carlo CI (optional with CLI flag) ---
    if calculate_ci:
        print(f"Starting Monte Carlo CI simulation ({n_permutations} iterations)...")
        ref_path_for_pools = real_reference_path if real_reference_path else variants_path
        af_pools = load_gene_af_pools(ref_path_for_pools)

        raw_lf = pl.scan_parquet(variants_path).select([gene_col, "raw_score"])

        # Normalize Gene ID again
        if gene_col == "gene_norm":
            raw_lf = raw_lf.with_columns(gene_norm=pl.col("gene_id").str.split(".").list.get(0))

        raw_df = raw_lf.collect()
        partitioned = raw_df.partition_by(gene_col, as_dict=True)

        result_rows = []
        for gid, sub_df in partitioned.items():
            if gid is None: continue
            stats = calculate_vg_ci_struct(sub_df["raw_score"], gid, af_pools, n_iter=n_permutations)
            stats[gene_col] = gid
            result_rows.append(stats)

        ci_results = pl.DataFrame(result_rows)
        df_agg = df_agg.join(ci_results, left_on=gene_col, right_on=gene_col, how="left")
        df_agg = df_agg.with_columns([
            pl.col("vg_perm_mean").fill_null(0.0),
            pl.col("vg_perm_p05").fill_null(0.0),
            pl.col("vg_perm_p95").fill_null(0.0),
        ])

    gtf_meta = gtf.select(
        pl.all().exclude(["tss", "strand", "start", "end"])
    ).unique("gene_id")

    # Handle MANE
    if "gene_id" in mane.columns:
        mane_meta = mane.select(
            [pl.col("gene_id"), pl.col("mane_transcript_id"), pl.lit(True).alias("is_mane")]
        ).drop_nulls("gene_id").unique("gene_id")
    else:
        mane_meta = mane.select(
            [pl.col("Ensembl_Gene").str.split(".").list.get(0).alias("gene_id"), pl.lit(True).alias("is_mane")]
        ).drop_nulls("gene_id").unique("gene_id")

    tpm_meta = tpm.select(["gene_id", "tpm_muscle"]).unique("gene_id")
    vgh_meta = vgh.unique("gene_id")

    enriched = (
        df_agg.rename({gene_col: "gene_id"})
        .join(gtf_meta, on="gene_id", how="left")
        .join(mane_meta, on="gene_id", how="left")
        .join(tpm_meta, on="gene_id", how="left")
        .join(vgh_meta, on="gene_id", how="left")
    )

    enriched = enriched.with_columns(
        cv_effect=pl.col("std_abs_effect") / pl.col("mean_abs_effect"),
        frac_high_impact_05=pl.col("n_high_impact_gt05") / pl.col("n_variants"),
        frac_high_impact_10=pl.col("n_high_impact_gt1") / pl.col("n_variants"),
        variants_per_kb=(pl.col("n_variants") / (pl.col("genomic_length") / 1000.0)).fill_nan(0.0),

        # Spatial Enrichment: Promoter Density vs Global Density
        # (Avoid div by zero)
        enrich_promoter_vg=(
                (pl.col("vg_promoter_core") / pl.col("n_variants_promoter_core").clip(1)) /
                (pl.col("vg_predicted") / pl.col("n_variants").clip(1))
        ).fill_nan(0.0)
    )

    # 14. Write Output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.write_parquet(out_path, compression="zstd")
    print(f"Saved aggregated gene metrics to {out_path}")