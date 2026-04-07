#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from alphagenome_shared import (
    DEFAULT_SEQ_LEN_1MB,
    build_layout,
    make_interval,
    read_table,
    rows_to_json,
    write_table_bundle,
)


REQUIRED_JOBS_COLUMNS = {
    "selection_id",
    "selection_group",
    "group_rank",
    "gene_id",
    "gene_symbol",
    "singleton_variant_id",
    "singleton_chrom",
    "singleton_pos",
    "singleton_ref",
    "singleton_alt",
}
REQUIRED_BACKGROUND_COLUMNS = {
    "selection_id",
    "selection_group",
    "group_rank",
    "gene_id",
    "gene_symbol",
    "singleton_variant_id",
    "singleton_chrom",
    "singleton_pos",
    "singleton_ref",
    "singleton_alt",
    "susie_variant_key",
    "susie_chrom",
    "susie_pos",
    "susie_ref",
    "susie_alt",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare scorer-ready singleton/haplotype manifests from upstream parquet inputs",
    )
    parser.add_argument("--dataset-id", help="dataset id; falls back to DATASET_ID")
    parser.add_argument("--sample-id", help="sample id; falls back to SAMPLE_ID")
    parser.add_argument("--root-dir", help="dataset root; falls back to ROOT_DIR/PDC_TMP")
    parser.add_argument(
        "--jobs-parquet",
        required=True,
        help="upstream singleton selection parquet (e.g. alphagenome_jobs.parquet)",
    )
    parser.add_argument(
        "--background-parquet",
        required=True,
        help="upstream haplotype background parquet (e.g. selected_singletons_haplotype_background.parquet)",
    )
    parser.add_argument(
        "--seq-len",
        type=int,
        default=DEFAULT_SEQ_LEN_1MB,
        help=f"AlphaGenome interval width to materialize [default: {DEFAULT_SEQ_LEN_1MB}]",
    )
    parser.add_argument(
        "--output-dir",
        help="override output directory; defaults to layout 01_inputs/singleton_haplotype",
    )
    return parser.parse_args()


def require_columns(df: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{label} missing required columns: {missing}")


def build_background_manifest(background_df: pd.DataFrame) -> pd.DataFrame:
    if background_df.empty:
        return pd.DataFrame(
            columns=[
                "selection_id",
                "selection_group",
                "group_rank",
                "gene_id",
                "gene_symbol",
                "singleton_variant_id",
                "background_variant_rank",
                "background_variant_id",
                "background_chrom",
                "background_pos",
                "background_ref",
                "background_alt",
            ]
        )

    out = background_df.copy()
    out = out[~out["matches_singleton"].fillna(False)].copy() if "matches_singleton" in out.columns else out
    out["background_variant_id"] = out["susie_variant_key"].astype(str)
    out["background_chrom"] = out["susie_chrom"].astype(str)
    out["background_pos"] = pd.to_numeric(out["susie_pos"], errors="coerce").astype("Int64")
    out["background_ref"] = out["susie_ref"].astype(str)
    out["background_alt"] = out["susie_alt"].astype(str)
    out["background_variant_rank"] = (
        out.sort_values(
            ["selection_id", "pip", "distance_to_singleton_bp", "background_variant_id"],
            ascending=[True, False, True, True],
        )
        .groupby("selection_id")
        .cumcount()
        .add(1)
    )
    keep_columns = [
        "selection_id",
        "selection_group",
        "group_rank",
        "gene_id",
        "gene_symbol",
        "singleton_variant_id",
        "singleton_chrom",
        "singleton_pos",
        "singleton_ref",
        "singleton_alt",
        "background_variant_rank",
        "background_variant_id",
        "background_chrom",
        "background_pos",
        "background_ref",
        "background_alt",
        "pip",
        "af",
        "cs_id",
        "cs_size",
        "afc",
        "afc_se",
        "distance_to_singleton_bp",
        "afc_direction",
        "matches_singleton",
    ]
    keep_columns = [column for column in keep_columns if column in out.columns]
    return out.loc[:, keep_columns].sort_values(["selection_id", "background_variant_rank"]).reset_index(drop=True)


def build_jobs_manifest(
    jobs_df: pd.DataFrame,
    background_manifest: pd.DataFrame,
    *,
    seq_len: int,
) -> pd.DataFrame:
    manifest = jobs_df.copy()
    manifest["selection_id"] = manifest["selection_id"].astype("int64")
    manifest["group_rank"] = manifest["group_rank"].astype("int64")
    manifest["variant_id"] = manifest["singleton_variant_id"].astype(str)
    manifest["chrom"] = manifest["singleton_chrom"].astype(str)
    manifest["pos"] = pd.to_numeric(manifest["singleton_pos"], errors="coerce").astype("int64")
    manifest["ref"] = manifest["singleton_ref"].astype(str)
    manifest["alt"] = manifest["singleton_alt"].astype(str)
    manifest["interval_chrom"] = manifest["chrom"]
    manifest["interval_width"] = seq_len

    interval_records = manifest.apply(
        lambda row: make_interval(seq_len, row["chrom"], row["pos"]),
        axis=1,
    )
    manifest["interval_start"] = interval_records.map(lambda interval: int(interval.start))
    manifest["interval_end"] = interval_records.map(lambda interval: int(interval.end))
    manifest["interval_str"] = interval_records.map(
        lambda interval: f"{interval.chromosome}:{int(interval.start)}-{int(interval.end)}:."
    )

    background_lookup = {}
    if not background_manifest.empty:
        for selection_id, group in background_manifest.groupby("selection_id", sort=True):
            records = group.to_dict(orient="records")
            background_lookup[int(selection_id)] = records

    manifest["background_variant_count_manifest"] = manifest["selection_id"].map(
        lambda selection_id: len(background_lookup.get(int(selection_id), []))
    )
    manifest["background_variants_json"] = manifest["selection_id"].map(
        lambda selection_id: rows_to_json(background_lookup.get(int(selection_id), []))
    )
    if "n_haplotype_background_variants" in manifest.columns:
        manifest["background_count_matches_upstream"] = (
            manifest["n_haplotype_background_variants"].fillna(0).astype(int)
            == manifest["background_variant_count_manifest"].fillna(0).astype(int)
        )

    ordered = [
        "selection_id",
        "selection_group",
        "group_rank",
        "gene_id",
        "gene_symbol",
        "variant_id",
        "singleton_variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "singleton_chrom",
        "singleton_pos",
        "singleton_ref",
        "singleton_alt",
        "singleton_af",
        "raw_score",
        "abs_raw_score",
        "interval_chrom",
        "interval_start",
        "interval_end",
        "interval_width",
        "interval_str",
        "background_variant_count_manifest",
        "background_variants_json",
    ]
    ordered += [
        column
        for column in manifest.columns
        if column not in ordered
    ]
    return manifest.loc[:, ordered].sort_values(["selection_group", "group_rank", "selection_id"]).reset_index(drop=True)


def main() -> None:
    args = parse_args()
    layout = build_layout(
        dataset_id=args.dataset_id,
        sample_id=args.sample_id,
        root_dir=args.root_dir,
    )
    output_dir = Path(args.output_dir) if args.output_dir else layout.singleton_haplotype_inputs_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    jobs_df = read_table(args.jobs_parquet)
    background_df = read_table(args.background_parquet)
    require_columns(jobs_df, REQUIRED_JOBS_COLUMNS, "jobs parquet")
    require_columns(background_df, REQUIRED_BACKGROUND_COLUMNS, "background parquet")

    background_manifest = build_background_manifest(background_df)
    jobs_manifest = build_jobs_manifest(jobs_df, background_manifest, seq_len=args.seq_len)

    jobs_parquet = output_dir / layout.singleton_haplotype_jobs_parquet.name
    jobs_tsv = output_dir / layout.singleton_haplotype_jobs_tsv.name
    background_parquet = output_dir / layout.singleton_haplotype_background_parquet.name
    background_tsv = output_dir / layout.singleton_haplotype_background_tsv.name

    write_table_bundle(jobs_manifest, parquet_path=jobs_parquet, tsv_path=jobs_tsv)
    write_table_bundle(background_manifest, parquet_path=background_parquet, tsv_path=background_tsv)

    print(f"Wrote {jobs_parquet}")
    print(f"Wrote {background_parquet}")


if __name__ == "__main__":
    main()
