from __future__ import annotations

import json
from pathlib import Path

import polars as pl


def compare_vg_medians(
    gene_parquet: Path,
    observed_col: str = "vg_predicted",
    sanity_col: str = "vg_predicted_perm_sanity",
    out_path: Path | None = None,
) -> dict[str, float | int | str]:
    """Compare observed vs sanity-permuted Vg medians and optionally write a JSON report."""
    gene_parquet = gene_parquet.resolve()
    if not gene_parquet.exists():
        raise FileNotFoundError(f"Gene parquet not found: {gene_parquet}")

    lf = pl.scan_parquet(gene_parquet)
    schema_cols = set(lf.collect_schema().names())
    missing = [c for c in [observed_col, sanity_col] if c not in schema_cols]
    if missing:
        raise ValueError(
            f"Cannot run sanity check for {gene_parquet}: missing columns {missing}"
        )

    df = lf.select([observed_col, sanity_col]).collect()
    valid = df.drop_nulls([observed_col, sanity_col])

    n_total = int(df.height)
    n_valid = int(valid.height)

    if n_valid == 0:
        result: dict[str, float | int | str] = {
            "gene_parquet": str(gene_parquet),
            "observed_col": observed_col,
            "sanity_col": sanity_col,
            "n_genes_total": n_total,
            "n_genes_valid": n_valid,
            "median_observed": 0.0,
            "median_sanity": 0.0,
            "median_delta": 0.0,
            "median_ratio_obs_over_sanity": 0.0,
            "relative_delta_pct": 0.0,
        }
    else:
        median_obs = float(valid[observed_col].median())
        median_sanity = float(valid[sanity_col].median())
        median_delta = median_obs - median_sanity
        ratio = float(median_obs / median_sanity) if median_sanity != 0 else 0.0
        relative_delta_pct = float((median_delta / median_sanity) * 100.0) if median_sanity != 0 else 0.0

        result = {
            "gene_parquet": str(gene_parquet),
            "observed_col": observed_col,
            "sanity_col": sanity_col,
            "n_genes_total": n_total,
            "n_genes_valid": n_valid,
            "median_observed": median_obs,
            "median_sanity": median_sanity,
            "median_delta": median_delta,
            "median_ratio_obs_over_sanity": ratio,
            "relative_delta_pct": relative_delta_pct,
        }

    if out_path is not None:
        out_path = out_path.resolve()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(result, indent=2) + "\n")

    return result
