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
        # If missing (e.g. pure ISM without annotation), fill 0
        lf = lf.with_columns(pl.lit(0.0).alias("AF"))
        schema_cols.add("AF")
    
    # Drop extra AF aliases but KEEP perm_AF
    extras = [c for c in af_candidates if c != "AF" and c in schema_cols]
    if extras:
        lf = lf.drop(extras)

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

    # 5. Attach Spatial
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

    # 7. Calculations (Calculate BOTH Vg metrics)
    lf = lf.with_columns(
        dist_to_tss=pl.when((pl.col("POS").is_not_null()) & (pl.col("tss").is_not_null()))
        .then((pl.col("POS") - pl.col("tss")).abs())
        .otherwise(None),
        
        dist_signed=pl.when(
            (pl.col("POS").is_not_null()) & (pl.col("tss").is_not_null()) & (pl.col("strand").is_not_null())
        )
        .then(
            pl.when(pl.col("strand") == "+")
            .then(pl.col("POS") - pl.col("tss"))
            .otherwise(pl.col("tss") - pl.col("POS"))
        )
        .otherwise(None),
        
        # 1. Real Vg (from AF)
        vg_contribution=pl.when(
            pl.col("raw_score").is_not_null() & pl.col("AF").is_not_null()
        )
        .then(2.0 * pl.col("AF") * (1.0 - pl.col("AF")) * pl.col("raw_score").pow(2))
        .otherwise(0.0),
        
        # 2. Permuted Vg (from perm_AF, if present)
        vg_contribution_perm=pl.when(
            pl.col("raw_score").is_not_null() & pl.col("perm_AF").is_not_null()
        )
        .then(2.0 * pl.col("perm_AF") * (1.0 - pl.col("perm_AF")) * pl.col("raw_score").pow(2))
        .otherwise(0.0) if has_perm_af else pl.lit(0.0),
    )

    # 8. Collect for Aggregation
    print("Collecting dataframe for aggregation...")
    df = lf.collect()

    ci_results = None
    if calculate_ci:
        print(f"Starting Monte Carlo CI simulation ({n_permutations} iterations)...")
        ref_path_for_pools = real_reference_path if real_reference_path else variants_path
        
        af_pools = load_gene_af_pools(ref_path_for_pools)

        if gene_col not in df.columns:
             pass 

        partitioned = df.partition_by(gene_col, as_dict=True)
        result_rows = []
        
        for gid, sub_df in partitioned.items():
            if gid is None: continue
            
            stats = calculate_vg_ci_struct(
                sub_df["raw_score"], 
                gid, 
                af_pools, 
                n_iter=n_permutations
            )
            stats[gene_col] = gid
            result_rows.append(stats)
            
        ci_results = pl.DataFrame(result_rows)

    # 9. Standard Aggregations
    
    agg_exprs = [
        # Basic counts
        pl.count().alias("n_variants"),
        
        # --- The Two Vg Metrics ---
        pl.col("vg_contribution").sum().alias("vg_predicted"),            # From AF
        pl.col("vg_contribution_perm").sum().alias("vg_predicted_perm"),  # From perm_AF
        
        pl.col("raw_score").pow(2).sum().alias("sum_sq_raw_score"),
        pl.col("raw_score").mean().alias("mean_raw_score"),
        
        # Global stats
        pl.col("abs_score").mean().alias("mean_abs_effect"),
        pl.col("abs_score").median().alias("median_abs_effect"),
        pl.col("abs_score").std().alias("std_abs_effect"),
        pl.col("abs_score").min().alias("min_abs_effect"),
        pl.col("abs_score").max().alias("max_abs_effect"),
        pl.col("abs_score").skew().alias("skewness_effect"),
        pl.col("abs_score").quantile(0.9).alias("q90_abs_effect"),
        
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

        # Spatial bins - Promoter
        pl.col("abs_score").filter(pl.col("dist_to_tss") <= 2000).mean().alias("mean_abs_promoter"),
        (pl.col("dist_to_tss") <= 2000).sum().alias("n_variants_promoter"),

        # Upstream 2kb
        pl.col("abs_score").filter((pl.col("dist_signed") >= -2000) & (pl.col("dist_signed") < 0)).mean().alias("mean_abs_up2kb"),
        ((pl.col("dist_signed") >= -2000) & (pl.col("dist_signed") < 0)).sum().alias("n_variants_up2kb"),

        # Upstream 10kb
        pl.col("abs_score").filter((pl.col("dist_signed") >= -10000) & (pl.col("dist_signed") < 0)).mean().alias("mean_abs_up10kb"),
        ((pl.col("dist_signed") >= -10000) & (pl.col("dist_signed") < 0)).sum().alias("n_variants_up10kb"),

        # Upstream 100kb
        pl.col("abs_score").filter((pl.col("dist_signed") >= -100000) & (pl.col("dist_signed") < 0)).mean().alias("mean_abs_up100kb"),
        ((pl.col("dist_signed") >= -100000) & (pl.col("dist_signed") < 0)).sum().alias("n_variants_up100kb"),

        # Downstream 2kb
        pl.col("abs_score").filter((pl.col("dist_signed") > 0) & (pl.col("dist_signed") <= 2000)).mean().alias("mean_abs_down2kb"),
        ((pl.col("dist_signed") > 0) & (pl.col("dist_signed") <= 2000)).sum().alias("n_variants_down2kb"),

        # Gene body
        pl.col("abs_score").filter(
            (pl.col("POS") >= pl.col("start")) & (pl.col("POS") <= pl.col("end")) & (pl.col("dist_to_tss") > 2000)
        ).mean().alias("mean_abs_gene_body"),
        ((pl.col("POS") >= pl.col("start")) & (pl.col("POS") <= pl.col("end")) & (pl.col("dist_to_tss") > 2000)).sum().alias("n_variants_gene_body"),
    ]

    # Perform aggregation
    genes = df.group_by(gene_col).agg(agg_exprs).rename({gene_col: "gene_id"})

    # 10. Merge CI Results
    if calculate_ci and ci_results is not None:
        genes = genes.join(ci_results, left_on="gene_id", right_on=gene_col, how="left")
        genes = genes.with_columns([
            pl.col("vg_perm_mean").fill_null(0.0),
            pl.col("vg_perm_p05").fill_null(0.0),
            pl.col("vg_perm_p95").fill_null(0.0),
        ])

    # 11. Enrich with Metadata 
    gtf_meta = gtf.select(
        pl.all().exclude(["tss", "strand", "start", "end"])
    ).unique("gene_id")

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
        genes.join(gtf_meta, on="gene_id", how="left")
        .join(mane_meta, on="gene_id", how="left")
        .join(tpm_meta, on="gene_id", how="left")
        .join(vgh_meta, on="gene_id", how="left")
    )

    # 11. post-aggregation calculations
    enriched = enriched.with_columns(
        cv_effect = pl.col("std_abs_effect") / pl.col("mean_abs_effect"), 
        frac_high_impact_05 = pl.col("n_high_impact_gt05") / pl.col("n_variants"),
        frac_high_impact_10 = pl.col("n_high_impact_gt1") / pl.col("n_variants"),
        enrich_up_vs_down_2kb = (pl.col("mean_abs_up2kb") / pl.col("mean_abs_down2kb")).fill_nan(0.0),
        enrich_up_vs_body = (pl.col("mean_abs_up2kb") / pl.col("mean_abs_gene_body")).fill_nan(0.0),
        variants_per_kb = (pl.col("n_variants") / (pl.col("genomic_length") / 1000.0)).fill_nan(0.0)
    )

    # 12. Write
    out_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.write_parquet(out_path, compression="zstd")