from __future__ import annotations

from pathlib import Path
import numpy as np
import polars as pl

from modules.external_data_loader import ExternalDataLoader
from modules.normalisation_helper import strip_ensembl_version
from modules.permuted_af import load_gene_af_pools


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


def calc_architecture_stats(struct_list: list[dict] | pl.Series, suffix: str = "") -> dict:
    """computes N90, N85, and properties of the 'driver' variants.
    
    args:
        struct_list: list of dicts or polars series [{'v': vg_contribution, 'e': abs_score}, ...]
        suffix: suffix to append to output keys (e.g., "_perm")
    """
    # handle polars series ambiguity
    if struct_list is None or len(struct_list) == 0:
        return {
            f"N90{suffix}": 0, f"N85{suffix}": 0,
            f"variance_N90{suffix}": 0.0, f"cv_effect_N90{suffix}": 0.0,
            f"mean_effect_N90{suffix}": 0.0
        }

    # extract arrays
    # note: iterating over a series of structs yields dictionaries, so this works for both list and series
    v = np.array([x['v'] for x in struct_list])
    e = np.array([x['e'] for x in struct_list])

    # sort descending by variance contribution (v)
    sort_idx = np.argsort(v)[::-1]
    v_sorted = v[sort_idx]
    e_sorted = e[sort_idx]

    total_vg = np.sum(v_sorted)
    if total_vg == 0:
        return {
            f"N90{suffix}": 0, f"N85{suffix}": 0,
            f"variance_N90{suffix}": 0.0, f"cv_effect_N90{suffix}": 0.0,
            f"mean_effect_N90{suffix}": 0.0
        }

    cumsum = np.cumsum(v_sorted)

    # logic for N90 (90% of variance)
    idx_90 = np.searchsorted(cumsum, 0.90 * total_vg)

    # subset of variants driving this 90%
    subset_e_90 = e_sorted[:idx_90 + 1]

    # stats for N90 subset
    n90_count = idx_90 + 1
    mean_90 = np.mean(subset_e_90) if len(subset_e_90) > 0 else 0.0
    std_90 = np.std(subset_e_90) if len(subset_e_90) > 0 else 0.0
    cv_90 = (std_90 / mean_90) if mean_90 > 0 else 0.0

    # logic for N85
    idx_85 = np.searchsorted(cumsum, 0.85 * total_vg)
    n85_count = idx_85 + 1

    return {
        f"N90{suffix}": int(n90_count),
        f"N85{suffix}": int(n85_count),
        f"variance_N90{suffix}": float(cumsum[idx_90]),
        f"cv_effect_N90{suffix}": float(cv_90),
        f"mean_effect_N90{suffix}": float(mean_90)
    }


def _calculate_entropy(scores: pl.Series) -> float:
    """Calculate Shannon entropy of effect sizes."""
    scores_arr = scores.to_numpy()
    scores_arr = scores_arr[scores_arr > 0]  # Remove zeros
    if len(scores_arr) == 0:
        return 0.0
    p = scores_arr / np.sum(scores_arr)
    return float(-np.sum(p * np.log(p + 1e-10)))  # Add epsilon for numerical stability


def _calculate_af_gradient(struct_series) -> float:
    """Calculate Spearman correlation between distance and AF."""
    from scipy.stats import spearmanr
    
    df = struct_series.struct.unnest()
    dist = df["dist_to_tss"].to_numpy()
    af = df["AF"].to_numpy()
    
    # Remove nulls
    valid_mask = ~np.isnan(dist) & ~np.isnan(af)
    dist = dist[valid_mask]
    af = af[valid_mask]
    
    if len(dist) < 3:  # Need at least 3 points for correlation
        return 0.0
    
    corr, _ = spearmanr(np.abs(dist), af)
    return float(corr) if not np.isnan(corr) else 0.0


# --- 4. Spatial Window Expression Generator ---
def get_window_exprs(windows: dict[str, tuple[int, int]], vg_col: str = "vg_contribution", suffix: str = "") -> list[pl.Expr]:
    """generates agg expressions for N, Mean Abs, and Sum Vg for spatial windows.
    
    args:
        windows: dict mapping window names to (start, end) ranges
        vg_col: which vg column to use for variance calculations
        suffix: suffix to append to vg column names (e.g., "_perm")
    """
    exprs = []
    for name, (start, end) in windows.items():
        # condition: start <= dist_signed < end
        cond = (pl.col("dist_signed") >= start) & (pl.col("dist_signed") < end)

        exprs.extend([
            pl.col("abs_score").filter(cond).mean().alias(f"mean_abs_{name}"),
            pl.col(vg_col).filter(cond).sum().alias(f"vg_{name}{suffix}"),
            cond.sum().alias(f"n_variants_{name}")
        ])
    return exprs


def get_window_vg_exprs(windows: dict[str, tuple[int, int]], vg_col: str = "vg_contribution", suffix: str = "") -> list[pl.Expr]:
    """generates only vg sum expressions for spatial windows (for perm metrics)."""
    exprs = []
    for name, (start, end) in windows.items():
        cond = (pl.col("dist_signed") >= start) & (pl.col("dist_signed") < end)
        exprs.append(pl.col(vg_col).filter(cond).sum().alias(f"vg_{name}{suffix}"))
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
        is_synthetic: bool = False,
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

    # --- 8. Determine variance mode ---
    # in synthetic mode, use only perm_AF-based metrics with _perm suffix
    if is_synthetic:
        vg_col = "vg_contribution_perm"
        vg_suffix = "_perm"
        vg_label = "vg_predicted_perm"
    else:
        vg_col = "vg_contribution"
        vg_suffix = ""
        vg_label = "vg_predicted"

    # --- 9. Build Aggregation Expressions ---
    agg_exprs = [
        pl.len().alias("n_variants"),

        pl.col(vg_col).sum().alias(vg_label),
        
        # always include perm if available (for sanity checks in real data)
        pl.col("vg_contribution_perm").sum().alias("vg_predicted_perm") if has_perm_af and not is_synthetic else None,

        # AF cutoff metrics: variant counts
        pl.col("AF").filter(pl.col("AF") >= 0.05).count().alias("n_variants_common"),
        pl.col("AF").filter(pl.col("AF") < 0.05).count().alias("n_variants_rare"),
        
        # AF cutoff metrics: variance contributions
        pl.col("vg_contribution").filter(pl.col("AF") >= 0.05).sum().alias("vg_common"),
        pl.col("vg_contribution").filter(pl.col("AF") < 0.05).sum().alias("vg_rare"),
    ]
    
    # Add permuted AF cutoff metrics if available
    if has_perm_af:
        agg_exprs.extend([
            pl.col("perm_AF").filter(pl.col("perm_AF") >= 0.05).count().alias("n_variants_common_perm"),
            pl.col("perm_AF").filter(pl.col("perm_AF") < 0.05).count().alias("n_variants_rare_perm"),
            pl.col("vg_contribution_perm").filter(pl.col("perm_AF") >= 0.05).sum().alias("vg_common_perm"),
            pl.col("vg_contribution_perm").filter(pl.col("perm_AF") < 0.05).sum().alias("vg_rare_perm"),
        ])
    
    agg_exprs.extend([
        # architecture data using the appropriate vg column
        pl.struct([
            pl.col(vg_col).alias("v"), 
            pl.col("abs_score").alias("e")
        ]).alias("arch_data"),

        pl.col("raw_score").pow(2).sum().alias("sum_sq_raw_score"),
        pl.col("raw_score").mean().alias("mean_raw_score"),

        # global stats (unchanged, based on raw_score)
        pl.col("abs_score").mean().alias("mean_abs_effect"),
        pl.col("abs_score").median().alias("median_abs_effect"),
        pl.col("abs_score").std().alias("std_abs_effect"),
        pl.col("abs_score").min().alias("min_abs_effect"),
        pl.col("abs_score").max().alias("max_abs_effect"),
        pl.col("abs_score").skew().alias("skewness_effect"),
        pl.col("abs_score").quantile(0.99).alias("q99_abs_effect"),

        # id tracking
        pl.col("variant_id").sort_by("raw_score").first().alias("min_variant_id"),
        pl.col("raw_score").min().alias("min_variant_score"),
        pl.col("variant_id").sort_by("raw_score").last().alias("max_variant_id"),
        pl.col("raw_score").max().alias("max_variant_score"),

        # distance stats
        pl.col("dist_to_tss").mean().alias("mean_dist_to_tss"),
        pl.col("dist_to_tss").median().alias("median_dist_to_tss"),
        pl.col("dist_to_tss").min().alias("min_dist_to_tss"),
        pl.col("dist_to_tss").max().alias("max_dist_to_tss"),

        # high impact counts
        (pl.col("abs_score") > 0.5).sum().alias("n_high_impact_gt05"),
        (pl.col("abs_score") > 1.0).sum().alias("n_high_impact_gt1"),
        
        # Spatial distribution: median distance to TSS for high-impact variants
        pl.col("dist_to_tss").filter(pl.col("abs_score") > 0.1).median().alias("median_dist_tss_high_abs"),
        pl.col("dist_to_tss").filter(pl.col("abs_score") >= pl.col("abs_score").quantile(0.90)).median().alias("median_dist_tss_high_rel"),
        
        # Selection signature: weighted mean effect
        (pl.col("AF") * pl.col("abs_score")).sum().alias("_af_times_abs"),
        pl.col("AF").sum().alias("_af_sum"),
        
        # Collect data for post-aggregation metrics
        pl.col("abs_score").alias("scores_all"),
        pl.struct(["dist_to_tss", "AF"]).alias("_dist_af_struct"),
        
        # Collect abs_score lists per window for entropy calculation (will be processed after aggregation)
        pl.col("abs_score").filter((pl.col("dist_signed") >= -200) & (pl.col("dist_signed") < 200)).alias("scores_promoter_core"),
        pl.col("abs_score").filter((pl.col("dist_signed") >= -2000) & (pl.col("dist_signed") < -200)).alias("scores_proximal_upstream"),
        pl.col("abs_score").filter((pl.col("dist_signed") >= -10000) & (pl.col("dist_signed") < -2000)).alias("scores_distal_upstream"),
        pl.col("abs_score").filter((pl.col("dist_signed") >= 200) & (pl.col("dist_signed") < 2000)).alias("scores_down_proximal"),
        pl.col("abs_score").filter((pl.col("dist_signed") >= 2000) & (pl.col("dist_signed") < 10000)).alias("scores_down_distal"),
    ])
    
    # Add permuted AF-specific metrics if available
    if has_perm_af:
        agg_exprs.extend([
            (pl.col("perm_AF") * pl.col("abs_score")).sum().alias("_perm_af_times_abs"),
            pl.col("perm_AF").sum().alias("_perm_af_sum"),
            
            # Depletion scores
            (pl.col("AF").filter(pl.col("abs_score") >= pl.col("abs_score").quantile(0.90)).mean()).alias("_af_high_impact"),
            (pl.col("perm_AF").filter(pl.col("abs_score") >= pl.col("abs_score").quantile(0.90)).mean()).alias("_perm_af_high_impact"),
            
            (pl.col("AF").filter(pl.col("AF") >= 0.05).mean()).alias("_af_common"),
            (pl.col("perm_AF").filter(pl.col("perm_AF") >= 0.05).mean()).alias("_perm_af_common"),
            
            (pl.col("AF").filter(pl.col("AF") < 0.05).mean()).alias("_af_rare"),
            (pl.col("perm_AF").filter(pl.col("perm_AF") < 0.05).mean()).alias("_perm_af_rare"),
            
            # Depletion scores per spatial window
            (pl.col("AF").filter((pl.col("dist_signed") >= -200) & (pl.col("dist_signed") < 200)).mean()).alias("_af_promoter_core"),
            (pl.col("perm_AF").filter((pl.col("dist_signed") >= -200) & (pl.col("dist_signed") < 200)).mean()).alias("_perm_af_promoter_core"),
            
            (pl.col("AF").filter((pl.col("dist_signed") >= -2000) & (pl.col("dist_signed") < -200)).mean()).alias("_af_proximal_upstream"),
            (pl.col("perm_AF").filter((pl.col("dist_signed") >= -2000) & (pl.col("dist_signed") < -200)).mean()).alias("_perm_af_proximal_upstream"),
            
            (pl.col("AF").filter((pl.col("dist_signed") >= -10000) & (pl.col("dist_signed") < -2000)).mean()).alias("_af_distal_upstream"),
            (pl.col("perm_AF").filter((pl.col("dist_signed") >= -10000) & (pl.col("dist_signed") < -2000)).mean()).alias("_perm_af_distal_upstream"),
            
            (pl.col("AF").filter((pl.col("dist_signed") >= 200) & (pl.col("dist_signed") < 2000)).mean()).alias("_af_down_proximal"),
            (pl.col("perm_AF").filter((pl.col("dist_signed") >= 200) & (pl.col("dist_signed") < 2000)).mean()).alias("_perm_af_down_proximal"),
            
            (pl.col("AF").filter((pl.col("dist_signed") >= 2000) & (pl.col("dist_signed") < 10000)).mean()).alias("_af_down_distal"),
            (pl.col("perm_AF").filter((pl.col("dist_signed") >= 2000) & (pl.col("dist_signed") < 10000)).mean()).alias("_perm_af_down_distal"),
        ])
    
    # filter out None entries
    agg_exprs = [e for e in agg_exprs if e is not None]

    # add dynamic spatial window expressions
    agg_exprs.extend(get_window_exprs(spatial_windows, vg_col=vg_col, suffix=vg_suffix))
    
    # for real data with perm_AF, also create _perm vg window metrics (just vg, not counts/means)
    if has_perm_af and not is_synthetic:
        agg_exprs.extend(get_window_vg_exprs(spatial_windows, vg_col="vg_contribution_perm", suffix="_perm"))

    # --- 10. Perform Aggregation ---
    print("Collecting and Aggregating Genes...")
    df_agg = lf.group_by(gene_col).agg(agg_exprs).collect()

    # --- 11. Post-Aggregation: Compute Architecture (Python UDF) ---
    print("Computing Architecture Metrics (N90, CV_N90)...")

    # architecture metrics with appropriate suffix
    arch_col_suffix = "_perm" if is_synthetic else ""
    
    arch_schema = pl.Struct({
        f"N90{arch_col_suffix}": pl.Int64, 
        f"N85{arch_col_suffix}": pl.Int64,
        f"variance_N90{arch_col_suffix}": pl.Float64, 
        f"cv_effect_N90{arch_col_suffix}": pl.Float64,
        f"mean_effect_N90{arch_col_suffix}": pl.Float64
    })

    df_agg = df_agg.with_columns(
        pl.col("arch_data").map_elements(
            lambda x: calc_architecture_stats(x, suffix=arch_col_suffix),
            return_dtype=arch_schema
        ).alias("arch_metrics")
    ).unnest("arch_metrics").drop("arch_data")

    # --- 12. Monte Carlo CI (optional with CLI flag) ---
    if calculate_ci:
        print(f"Starting Monte Carlo CI simulation ({n_permutations} iterations)...")
        ref_path_for_pools = real_reference_path if real_reference_path else variants_path
        af_pools = load_gene_af_pools(ref_path_for_pools)
        
        # ensure pool keys match normalized gene_ids
        af_pools = {k.split('.')[0]: v for k, v in af_pools.items()}

        # always read gene_id from file, then normalize if needed
        raw_lf = pl.scan_parquet(variants_path).select(["gene_id", "raw_score"])

        # create gene_norm column if that's what we're using
        if gene_col == "gene_norm":
            raw_lf = raw_lf.with_columns(gene_norm=pl.col("gene_id").str.split(".").list.get(0))

        raw_df = raw_lf.collect()
        partitioned = raw_df.partition_by(gene_col, maintain_order=True)

        result_rows = []
        for sub_df in partitioned:
            # extract gene_id from first row to ensure proper type
            gid = sub_df[0, gene_col]
            if gid is None: 
                continue
            stats = calculate_vg_ci_struct(sub_df["raw_score"], gid, af_pools, n_iter=n_permutations)
            stats[gene_col] = gid
            result_rows.append(stats)

        ci_results = pl.DataFrame(result_rows)
        df_agg = df_agg.join(ci_results, on=gene_col, how="left")
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

    # build enrichment columns with appropriate suffix
    promoter_col = f"vg_promoter_core{vg_suffix}"
    vg_global_col = vg_label
    enrich_col = f"enrich_promoter_vg{vg_suffix}"
    
    # Build post-aggregation computed columns
    computed_cols = [
        # Existing metrics
        pl.col("n_high_impact_gt05").truediv(pl.col("n_variants")).alias("frac_high_impact_05"),
        pl.col("n_high_impact_gt1").truediv(pl.col("n_variants")).alias("frac_high_impact_10"),
        (pl.col("n_variants") / (pl.col("genomic_length") / 1000.0)).fill_nan(0.0).alias("variants_per_kb"),
        
        # Weighted mean effect
        (pl.col("_af_times_abs") / pl.col("_af_sum")).fill_nan(0.0).alias("weighted_mean_effect"),
        
        # Gene-wide entropy and AF gradient
        pl.col("scores_all").map_elements(lambda s: _calculate_entropy(s), return_dtype=pl.Float64).alias("entropy_effect"),
        pl.col("_dist_af_struct").map_elements(lambda s: _calculate_af_gradient(s), return_dtype=pl.Float64).alias("af_gradient"),
        
        # Entropy per window
        pl.col("scores_promoter_core").map_elements(lambda s: _calculate_entropy(s), return_dtype=pl.Float64).alias("entropy_promoter_core"),
        pl.col("scores_proximal_upstream").map_elements(lambda s: _calculate_entropy(s), return_dtype=pl.Float64).alias("entropy_proximal_upstream"),
        pl.col("scores_distal_upstream").map_elements(lambda s: _calculate_entropy(s), return_dtype=pl.Float64).alias("entropy_distal_upstream"),
        pl.col("scores_down_proximal").map_elements(lambda s: _calculate_entropy(s), return_dtype=pl.Float64).alias("entropy_down_proximal"),
        pl.col("scores_down_distal").map_elements(lambda s: _calculate_entropy(s), return_dtype=pl.Float64).alias("entropy_down_distal"),
        
        # Proportion of vg in each window
        (pl.col(f"vg_promoter_core{vg_suffix}") / pl.col(vg_global_col)).fill_nan(0.0).alias(f"prop_vg_promoter_core{vg_suffix}"),
        (pl.col(f"vg_proximal_upstream{vg_suffix}") / pl.col(vg_global_col)).fill_nan(0.0).alias(f"prop_vg_proximal_upstream{vg_suffix}"),
        (pl.col(f"vg_distal_upstream{vg_suffix}") / pl.col(vg_global_col)).fill_nan(0.0).alias(f"prop_vg_distal_upstream{vg_suffix}"),
        (pl.col(f"vg_down_proximal{vg_suffix}") / pl.col(vg_global_col)).fill_nan(0.0).alias(f"prop_vg_down_proximal{vg_suffix}"),
        (pl.col(f"vg_down_distal{vg_suffix}") / pl.col(vg_global_col)).fill_nan(0.0).alias(f"prop_vg_down_distal{vg_suffix}"),
        
        # Enrichment scores per window
        ((pl.col(f"vg_promoter_core{vg_suffix}") / pl.col("n_variants_promoter_core").clip(1)) /
         (pl.col(vg_global_col) / pl.col("n_variants").clip(1))).fill_nan(0.0).alias(f"enrich_vg_promoter_core{vg_suffix}"),
        
        ((pl.col(f"vg_proximal_upstream{vg_suffix}") / pl.col("n_variants_proximal_upstream").clip(1)) /
         (pl.col(vg_global_col) / pl.col("n_variants").clip(1))).fill_nan(0.0).alias(f"enrich_vg_proximal_upstream{vg_suffix}"),
        
        ((pl.col(f"vg_distal_upstream{vg_suffix}") / pl.col("n_variants_distal_upstream").clip(1)) /
         (pl.col(vg_global_col) / pl.col("n_variants").clip(1))).fill_nan(0.0).alias(f"enrich_vg_distal_upstream{vg_suffix}"),
        
        ((pl.col(f"vg_down_proximal{vg_suffix}") / pl.col("n_variants_down_proximal").clip(1)) /
         (pl.col(vg_global_col) / pl.col("n_variants").clip(1))).fill_nan(0.0).alias(f"enrich_vg_down_proximal{vg_suffix}"),
        
        ((pl.col(f"vg_down_distal{vg_suffix}") / pl.col("n_variants_down_distal").clip(1)) /
         (pl.col(vg_global_col) / pl.col("n_variants").clip(1))).fill_nan(0.0).alias(f"enrich_vg_down_distal{vg_suffix}"),
    ]
    
    # Add permuted AF-specific computed columns if available (only in non-synthetic mode)
    # In synthetic mode, these are already created in the main computed_cols with _perm suffix
    if has_perm_af and not is_synthetic:
        computed_cols.extend([
            # Weighted mean effect for perm_AF
            (pl.col("_perm_af_times_abs") / pl.col("_perm_af_sum")).fill_nan(0.0).alias("weighted_mean_effect_perm"),
            
            # Depletion scores
            (pl.col("_af_high_impact") / pl.col("_perm_af_high_impact")).fill_nan(0.0).alias("depletion_high_impact"),
            (pl.col("_af_common") / pl.col("_perm_af_common")).fill_nan(0.0).alias("depletion_common"),
            (pl.col("_af_rare") / pl.col("_perm_af_rare")).fill_nan(0.0).alias("depletion_rare"),
            
            # Depletion scores per window
            (pl.col("_af_promoter_core") / pl.col("_perm_af_promoter_core")).fill_nan(0.0).alias("depletion_promoter_core"),
            (pl.col("_af_proximal_upstream") / pl.col("_perm_af_proximal_upstream")).fill_nan(0.0).alias("depletion_proximal_upstream"),
            (pl.col("_af_distal_upstream") / pl.col("_perm_af_distal_upstream")).fill_nan(0.0).alias("depletion_distal_upstream"),
            (pl.col("_af_down_proximal") / pl.col("_perm_af_down_proximal")).fill_nan(0.0).alias("depletion_down_proximal"),
            (pl.col("_af_down_distal") / pl.col("_perm_af_down_distal")).fill_nan(0.0).alias("depletion_down_distal"),
            
            # Proportions and enrichments for _perm
            (pl.col("vg_promoter_core_perm") / pl.col("vg_predicted_perm")).fill_nan(0.0).alias("prop_vg_promoter_core_perm"),
            (pl.col("vg_proximal_upstream_perm") / pl.col("vg_predicted_perm")).fill_nan(0.0).alias("prop_vg_proximal_upstream_perm"),
            (pl.col("vg_distal_upstream_perm") / pl.col("vg_predicted_perm")).fill_nan(0.0).alias("prop_vg_distal_upstream_perm"),
            (pl.col("vg_down_proximal_perm") / pl.col("vg_predicted_perm")).fill_nan(0.0).alias("prop_vg_down_proximal_perm"),
            (pl.col("vg_down_distal_perm") / pl.col("vg_predicted_perm")).fill_nan(0.0).alias("prop_vg_down_distal_perm"),
            
            ((pl.col("vg_promoter_core_perm") / pl.col("n_variants_promoter_core").clip(1)) /
             (pl.col("vg_predicted_perm") / pl.col("n_variants").clip(1))).fill_nan(0.0).alias("enrich_vg_promoter_core_perm"),
            
            ((pl.col("vg_proximal_upstream_perm") / pl.col("n_variants_proximal_upstream").clip(1)) /
             (pl.col("vg_predicted_perm") / pl.col("n_variants").clip(1))).fill_nan(0.0).alias("enrich_vg_proximal_upstream_perm"),
            
            ((pl.col("vg_distal_upstream_perm") / pl.col("n_variants_distal_upstream").clip(1)) /
             (pl.col("vg_predicted_perm") / pl.col("n_variants").clip(1))).fill_nan(0.0).alias("enrich_vg_distal_upstream_perm"),
            
            ((pl.col("vg_down_proximal_perm") / pl.col("n_variants_down_proximal").clip(1)) /
             (pl.col("vg_predicted_perm") / pl.col("n_variants").clip(1))).fill_nan(0.0).alias("enrich_vg_down_proximal_perm"),
            
            ((pl.col("vg_down_distal_perm") / pl.col("n_variants_down_distal").clip(1)) /
             (pl.col("vg_predicted_perm") / pl.col("n_variants").clip(1))).fill_nan(0.0).alias("enrich_vg_down_distal_perm"),
        ])
    
    enriched = enriched.with_columns(computed_cols)
    
    # Drop intermediate columns
    cols_to_drop = [
        "_af_times_abs", "_af_sum", "scores_all", "_dist_af_struct",
        "scores_promoter_core", "scores_proximal_upstream", "scores_distal_upstream",
        "scores_down_proximal", "scores_down_distal"
    ]
    
    if has_perm_af:
        cols_to_drop.extend([
            "_perm_af_times_abs", "_perm_af_sum",
            "_af_high_impact", "_perm_af_high_impact",
            "_af_common", "_perm_af_common",
            "_af_rare", "_perm_af_rare",
            "_af_promoter_core", "_perm_af_promoter_core",
            "_af_proximal_upstream", "_perm_af_proximal_upstream",
            "_af_distal_upstream", "_perm_af_distal_upstream",
            "_af_down_proximal", "_perm_af_down_proximal",
            "_af_down_distal", "_perm_af_down_distal",
        ])
    
    enriched = enriched.drop([c for c in cols_to_drop if c in enriched.columns])

    # 14. Write Output
    out_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.write_parquet(out_path, compression="zstd")
    print(f"Saved aggregated gene metrics to {out_path}")


def main() -> None:
    import argparse
    import os
    from modules.variant_deduplicator import deduplicate_by_gene_and_variant
    from modules.synthetic_variant_downsampler import load_real_variant_counts, downsample_to_real_counts
    from modules.permuted_af import add_perm_af_gene_aware
    from modules.snv_filter import filter_to_snvs

    parser = argparse.ArgumentParser(description='process and aggregate variant scores')
    parser.add_argument('--variants-parquet', type=Path, required=True, help='input variants parquet')
    parser.add_argument('--gene-list', type=Path, help='optional gene whitelist')
    parser.add_argument('--variant-out', type=Path, help='output variants parquet')
    parser.add_argument('--gene-out', type=Path, required=True, help='output genes parquet')
    parser.add_argument('--root-dir', type=Path, help='root directory for reference data')
    parser.add_argument('--real-reference', type=Path, help='real dataset reference for downsampling and CI')
    parser.add_argument('--deduplicate', action='store_true', help='deduplicate variants')
    parser.add_argument('--permute-af', action='store_true', help='generate perm_AF')
    parser.add_argument('--calc-ci', action='store_true', help='calculate confidence intervals')
    parser.add_argument('--synthetic', action='store_true', help='synthetic/null mode: use only perm_AF for all Vg metrics')
    parser.add_argument('--filter-snvs', action='store_true', help='filter to keep only SNVs (removes indels and MNVs)')
    
    args = parser.parse_args()

    base_root = Path(args.root_dir) if args.root_dir else Path(
        os.getenv('ROOT_DIR') or os.getenv('PDC_TMP') or '/cfs/klemming/scratch/m/mmarandi'
    )
    variant_path = args.variants_parquet.resolve()
    
    # infer dataset type
    is_ism = 'ism' in str(variant_path).lower() or 'null' in str(variant_path).lower()
    
    # process variants if requested
    if args.deduplicate:
        import sys
        print(f'Loading variants from {variant_path}...', flush=True)
        sys.stdout.flush()
        df = pl.read_parquet(variant_path)
        print(f'  Loaded {len(df):,} rows, {len(df.columns)} columns', flush=True)
        sys.stdout.flush()

        # filter to SNVs only (CRITICAL: do this before deduplication)
        if args.filter_snvs:
            print('Filtering to SNVs only (removing indels and MNVs)...')
            df = filter_to_snvs(df, verbose=True)
        
        # filter by gene list
        if args.gene_list:
            print(f'Filtering variants using whitelist: {args.gene_list}')
            whitelist_df = pl.read_csv(
                args.gene_list, 
                has_header=False, 
                new_columns=['gene_id'], 
                separator='\t'
            )
            whitelist = whitelist_df['gene_id'].to_list()
            df = df.filter(pl.col('gene_id').is_in(whitelist))
            print(f'  Rows after gene list filtering: {len(df):,}')
        
        print('Deduplicating by gene_id and variant_id...')
        df = deduplicate_by_gene_and_variant(df, verbose=True)
        
        # downsample if ism/null
        if is_ism:
            if not args.real_reference:
                raise ValueError('Flag --real-reference is required for ISM/NULL processing')
            
            real_ref = args.real_reference.resolve()
            if not real_ref.exists():
                raise FileNotFoundError(f'Real reference not found: {real_ref}')
            
            print(f'Loading real variant counts from {real_ref}...')
            real_counts = load_real_variant_counts(real_ref, filter_snvs=args.filter_snvs)
            
            print('Downsampling to match real counts...')
            df = downsample_to_real_counts(df, real_counts, seed=42, verbose=True)
            
            suffix = '_downsampled'
        else:
            suffix = '_dedup'
        
        # add snv filter suffix if applied
        if args.filter_snvs:
            suffix += '_snv'

        # permute af
        if args.permute_af:
            print('--- Starting AF Permutation ---')
            perm_ref_path = None
            
            if args.real_reference:
                perm_ref_path = args.real_reference.resolve()
            elif not is_ism: 
                print('  No external reference provided for real data; using input variants as AF source')
                perm_ref_path = variant_path
            
            if not perm_ref_path:
                raise ValueError('Cannot permute AF: No reference provided (use --real-reference)')

            if not perm_ref_path.exists():
                raise FileNotFoundError(f'Permutation reference not found: {perm_ref_path}')

            print(f'  Loading AF pools from {perm_ref_path}...')
            af_pools = load_gene_af_pools(perm_ref_path)

            print('  Sampling and assigning perm_AF...')
            df = add_perm_af_gene_aware(df, af_pools, seed=42)
            
            suffix += '_perm'

        # write processed variants
        if args.variant_out:
            processed_path = args.variant_out.resolve()
        else:
            stem = variant_path.stem
            if not stem.endswith(suffix):
                stem += suffix
            processed_path = variant_path.with_name(f'{stem}.parquet')

        print(f'Writing processed variants to {processed_path}...')
        df.write_parquet(processed_path)
        
        # update variant_path for aggregation
        variant_path = processed_path
    
    # aggregate genes
    gene_path = args.gene_out.resolve()
    gene_list = args.gene_list.resolve() if args.gene_list else None

    # determine reference for ci
    ci_ref = None
    if args.real_reference:
        ci_ref = args.real_reference.resolve()
    elif not is_ism:
        ci_ref = variant_path

    print(f'Aggregating genes from {variant_path}...')
    aggregate_genes(
        variant_path,
        gene_path,
        base_ref=base_root,
        is_ism=is_ism,
        gene_list_path=gene_list,
        calculate_ci=args.calc_ci,
        real_reference_path=ci_ref,
        n_permutations=1000,
        is_synthetic=args.synthetic,
    )
    print('Done.')


if __name__ == '__main__':
    main()