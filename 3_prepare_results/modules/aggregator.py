from __future__ import annotations

from pathlib import Path
import polars as pl

from .external_data_loader import ExternalDataLoader
from .normalisation_helper import strip_ensembl_version


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


def aggregate_genes(
    variants_path: Path,
    out_path: Path,
    base_ref: Path,
    is_ism: bool,
    gene_list_path: Path | None = None,
) -> None:
    """aggregate variant parquet to gene metrics with predicted vg and descriptive stats."""
    # 1. load metadata (eager)
    mane, gtf, tpm, vgh = _load_gene_meta(base_ref)

    # 2. scan variants (lazy)
    lf = pl.scan_parquet(variants_path)
    
    # collect schema once safely
    schema_cols = set(lf.collect_schema().names())

    # normalize/select AF source once (case-insensitive), then drop extras
    af_candidates = ["AF", "AF_x", "AF_y", "af", "af_x", "af_y"]
    af_source = next((c for c in af_candidates if c in schema_cols), None)

    if af_source and af_source != "AF":
        lf = lf.with_columns(pl.col(af_source).alias("AF"))
        schema_cols.add("AF")

    # allow missing AF in ism/null runs by injecting a zero column
    if "AF" not in schema_cols and is_ism:
        lf = lf.with_columns(pl.lit(0.0).alias("AF"))
        schema_cols.add("AF")

    # drop secondary AF columns to avoid duplicates later
    extras = [c for c in af_candidates if c != "AF" and c in schema_cols]
    if extras:
        lf = lf.drop(extras)
        schema_cols.difference_update(extras)

    # allow missing AF in ism/null runs by injecting a zero column
    if "AF" not in schema_cols and is_ism:
        lf = lf.with_columns(pl.lit(0.0).alias("AF"))
        schema_cols.add("AF")

    # 3. validation
    required_cols = {"raw_score"} if is_ism else {"raw_score", "AF"}
    missing = required_cols.difference(schema_cols)
    if missing:
        missing_cols = ", ".join(sorted(missing))
        raise ValueError(f"calculated vg requires columns: {missing_cols}. input has: {sorted(schema_cols)}")

    # 4. normalize gene id
    if "gene_norm" in schema_cols:
        gene_col = "gene_norm"
    elif "gene_id" in schema_cols:
        gene_col = "gene_id"
        lf = lf.with_columns(gene_norm=pl.col("gene_id").str.split(".").list.get(0))
        gene_col = "gene_norm"
    elif "gene_tag" in schema_cols:
        lf = lf.with_columns(gene_norm=pl.col("gene_tag").str.split(".").list.get(0))
        gene_col = "gene_norm"
    else:
        lf = lf.with_columns(gene_norm=pl.lit(None))
        gene_col = "gene_norm"

    # optional gene whitelist
    if gene_list_path and gene_list_path.exists():
        gene_whitelist = {
            strip_ensembl_version(line.strip())
            for line in gene_list_path.read_text().splitlines()
            if line.strip()
        }
        if gene_whitelist:
            lf = lf.filter(pl.col(gene_col).is_in(gene_whitelist))

    # 5. attach spatial info (lazy join)
    gtf_spatial = gtf.lazy().select(
        [pl.col("gene_id"), "tss", "strand", "start", "end"]
    ).unique("gene_id")
    
    lf = lf.join(gtf_spatial, left_on=gene_col, right_on="gene_id", how="left")

    # 6. deduplication strategy (select max impact track per variant)
    # calculate abs_score first to determine the winner
    lf = lf.with_columns(abs_score=pl.col("raw_score").abs())

    # dedupe variant rows keeping highest abs_score per variant
    # this reduces row count before expensive spatial logic
    lf = (
        lf.sort("abs_score", descending=True, nulls_last=True)
        .unique(subset=["variant_id"], keep="first", maintain_order=True)
    )

    # 7. post-dedupe calculations (spatial & variance)
    lf = lf.with_columns(
        # absolute distance to tss
        dist_to_tss=pl.when((pl.col("POS").is_not_null()) & (pl.col("tss").is_not_null()))
        .then((pl.col("POS") - pl.col("tss")).abs())
        .otherwise(None),
        
        # signed distance (negative = upstream, positive = downstream/body)
        dist_signed=pl.when(
            (pl.col("POS").is_not_null()) & (pl.col("tss").is_not_null()) & (pl.col("strand").is_not_null())
        )
        .then(
            pl.when(pl.col("strand") == "+")
            .then(pl.col("POS") - pl.col("tss"))
            .otherwise(pl.col("tss") - pl.col("POS"))
        )
        .otherwise(None),
        
        # genetic variance contribution (2pq beta^2)
        # calculated on the "winner" score only
        vg_contribution=pl.when(
            pl.col("raw_score").is_not_null() & pl.col("AF").is_not_null()
        )
        .then(2.0 * pl.col("AF") * (1.0 - pl.col("AF")) * pl.col("raw_score").pow(2))
        .otherwise(0.0),
    )

    # 8. define aggregations (comprehensive list)
    agg_exprs = [
        #  basic counts & vg
        pl.count().alias("n_variants"),
        pl.col("vg_contribution").sum().alias("vg_predicted"),
        pl.col("raw_score").pow(2).sum().alias("sum_sq_raw_score"),
        pl.col("raw_score").mean().alias("mean_raw_score"),
        
        #  global stats (exon/window)
        pl.col("abs_score").mean().alias("mean_abs_effect"),
        pl.col("abs_score").median().alias("median_abs_effect"),
        pl.col("abs_score").std().alias("std_abs_effect"),
        pl.col("abs_score").min().alias("min_abs_effect"),
        pl.col("abs_score").max().alias("max_abs_effect"),
        pl.col("abs_score").skew().alias("skewness_effect"),
        pl.col("abs_score").quantile(0.9).alias("q90_abs_effect"),
        
        #  id tracking (min/max score)
        pl.col("variant_id").sort_by("raw_score").first().alias("min_variant_id"),
        pl.col("raw_score").min().alias("min_variant_score"),
        pl.col("variant_id").sort_by("raw_score").last().alias("max_variant_id"),
        pl.col("raw_score").max().alias("max_variant_score"),

        #  distance stats
        pl.col("dist_to_tss").mean().alias("mean_dist_to_tss"),
        pl.col("dist_to_tss").median().alias("median_dist_to_tss"),
        pl.col("dist_to_tss").min().alias("min_dist_to_tss"),
        pl.col("dist_to_tss").max().alias("max_dist_to_tss"),

        #  high impact counts
        (pl.col("abs_score") > 0.5).sum().alias("n_high_impact_gt05"),
        (pl.col("abs_score") > 1.0).sum().alias("n_high_impact_gt1"),

        #  spatial bins
        # promoter (<= 2kb from tss)
        pl.col("abs_score").filter(pl.col("dist_to_tss") <= 2000).mean().alias("mean_abs_promoter"),
        (pl.col("dist_to_tss") <= 2000).sum().alias("n_variants_promoter"),

        # upstream 2kb
        pl.col("abs_score").filter((pl.col("dist_signed") >= -2000) & (pl.col("dist_signed") < 0)).mean().alias("mean_abs_up2kb"),
        ((pl.col("dist_signed") >= -2000) & (pl.col("dist_signed") < 0)).sum().alias("n_variants_up2kb"),

        # upstream 10kb
        pl.col("abs_score").filter((pl.col("dist_signed") >= -10000) & (pl.col("dist_signed") < 0)).mean().alias("mean_abs_up10kb"),
        ((pl.col("dist_signed") >= -10000) & (pl.col("dist_signed") < 0)).sum().alias("n_variants_up10kb"),

        # upstream 100kb
        pl.col("abs_score").filter((pl.col("dist_signed") >= -100000) & (pl.col("dist_signed") < 0)).mean().alias("mean_abs_up100kb"),
        ((pl.col("dist_signed") >= -100000) & (pl.col("dist_signed") < 0)).sum().alias("n_variants_up100kb"),

        # downstream 2kb
        pl.col("abs_score").filter((pl.col("dist_signed") > 0) & (pl.col("dist_signed") <= 2000)).mean().alias("mean_abs_down2kb"),
        ((pl.col("dist_signed") > 0) & (pl.col("dist_signed") <= 2000)).sum().alias("n_variants_down2kb"),

        # gene body (inside start/end, excluding promoter area > 2kb)
        pl.col("abs_score").filter(
            (pl.col("POS") >= pl.col("start")) & (pl.col("POS") <= pl.col("end")) & (pl.col("dist_to_tss") > 2000)
        ).mean().alias("mean_abs_gene_body"),
        ((pl.col("POS") >= pl.col("start")) & (pl.col("POS") <= pl.col("end")) & (pl.col("dist_to_tss") > 2000)).sum().alias("n_variants_gene_body"),
    ]

    # 9. aggregate
    genes = lf.group_by(gene_col).agg(agg_exprs).rename({gene_col: "gene_id"})

    # 10. enrich with metadata (lazy join)
    gtf_meta = gtf.lazy().drop(["tss", "strand", "start", "end"], strict=False).unique("gene_id")

    if "gene_id" in mane.columns:
        mane_meta = mane.lazy().select(
            [pl.col("gene_id"), pl.col("mane_transcript_id"), pl.lit(True).alias("is_mane")]
        ).drop_nulls("gene_id").unique("gene_id")
    else:
        mane_meta = mane.lazy().select(
            [pl.col("Ensembl_Gene").str.split(".").list.get(0).alias("gene_id"), pl.lit(True).alias("is_mane")]
        ).drop_nulls("gene_id").unique("gene_id")

    tpm_meta = tpm.lazy().select(["gene_id", "tpm_muscle"]).unique("gene_id")
    vgh_meta = vgh.lazy().unique("gene_id")

    enriched = (
        genes.join(gtf_meta, on="gene_id", how="left")
        .join(mane_meta, on="gene_id", how="left")
        .join(tpm_meta, on="gene_id", how="left")
        .join(vgh_meta, on="gene_id", how="left")
    )

    # 11. post-aggregation calculations (ratios & fractions)
    enriched = enriched.with_columns(
        # cv = std / mean
        cv_effect = pl.col("std_abs_effect") / pl.col("mean_abs_effect"),
        
        # fractions
        frac_high_impact_05 = pl.col("n_high_impact_gt05") / pl.col("n_variants"),
        frac_high_impact_10 = pl.col("n_high_impact_gt1") / pl.col("n_variants"),
        
        # enrichments (up vs down) - fill 0/0 nans
        enrich_up_vs_down_2kb = (pl.col("mean_abs_up2kb") / pl.col("mean_abs_down2kb")).fill_nan(0.0),
        enrich_up_vs_body = (pl.col("mean_abs_up2kb") / pl.col("mean_abs_gene_body")).fill_nan(0.0),
        
        # density (variants per kb) - genomic_length comes from gtf meta
        variants_per_kb = (pl.col("n_variants") / (pl.col("genomic_length") / 1000.0)).fill_nan(0.0)
    )

    # 12. write
    out_path.parent.mkdir(parents=True, exist_ok=True)
    enriched.sink_parquet(out_path, compression="zstd")