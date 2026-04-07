#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from alphagenome.models import variant_scorers

from alphagenome_shared import (
    ORG,
    RNA,
    SKELETAL_MUSCLE_ONTOLOGY_CURIES,
    add_track_identity,
    apply_variants_to_sequence,
    build_layout,
    build_track_rows_from_track_data,
    filter_skeletal_muscle_rows,
    filter_target_gene,
    load_api_client,
    make_variant,
    map_friendly,
    normalize_tidy,
    open_reference_reader,
    parse_variant_rows_json,
    read_table,
    reduce_predict_sequence_delta,
    write_table_bundle,
)
from alphagenome.data import genome


REQUIRED_MANIFEST_COLUMNS = {
    "selection_id",
    "selection_group",
    "group_rank",
    "gene_id",
    "gene_symbol",
    "variant_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "interval_chrom",
    "interval_start",
    "interval_end",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score selected singleton variants against reference and haplotype backgrounds",
    )
    parser.add_argument("--dataset-id", help="dataset id; falls back to DATASET_ID")
    parser.add_argument("--sample-id", help="sample id; falls back to SAMPLE_ID")
    parser.add_argument("--root-dir", help="dataset root; falls back to ROOT_DIR/PDC_TMP")
    parser.add_argument(
        "--jobs-manifest",
        help="prepared singleton jobs manifest parquet/tsv.gz; defaults to layout singleton_haplotype_jobs.parquet",
    )
    parser.add_argument(
        "--background-manifest",
        help="optional normalized background manifest parquet/tsv.gz; used as fallback if jobs manifest lacks JSON rows",
    )
    parser.add_argument(
        "--reference-fasta",
        help="reference FASTA for building haplotype sequences; required for --mode haplotype/both",
    )
    parser.add_argument("--reference-fai", help="optional FASTA index path")
    parser.add_argument(
        "--mode",
        choices=("both", "baseline", "haplotype"),
        default="both",
        help="which scoring branch(es) to run [default: both]",
    )
    parser.add_argument(
        "--env-file",
        default=Path(__file__).resolve().with_name(".env"),
        help="dotenv file containing API_KEY_PERSONAL",
    )
    parser.add_argument(
        "--ontology-curie",
        action="append",
        dest="ontology_curies",
        help="repeat to override the default skeletal-muscle ontology terms",
    )
    parser.add_argument(
        "--selection-id",
        action="append",
        type=int,
        dest="selection_ids",
        help="optional selection_id filter; can be repeated",
    )
    parser.add_argument(
        "--max-selections",
        type=int,
        help="process only the first N manifest rows after filtering",
    )
    parser.add_argument(
        "--output-dir",
        help="override results directory; defaults to layout 03_results/singleton_haplotype",
    )
    return parser.parse_args()


def require_columns(df: pd.DataFrame, required: set[str], label: str) -> None:
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"{label} missing required columns: {missing}")


def manifest_selection_metadata(row: pd.Series) -> dict[str, object]:
    return {
        "selection_id": int(row["selection_id"]),
        "selection_group": row["selection_group"],
        "group_rank": int(row["group_rank"]),
        "gene_id": row["gene_id"],
        "gene_symbol": row["gene_symbol"],
        "variant_id": row["variant_id"],
        "singleton_variant_id": row["variant_id"],
        "chrom": row["chrom"],
        "pos": int(row["pos"]),
        "ref": row["ref"],
        "alt": row["alt"],
    }


def interval_from_row(row: pd.Series) -> genome.Interval:
    return genome.Interval(
        chromosome=str(row["interval_chrom"]),
        start=int(row["interval_start"]),
        end=int(row["interval_end"]),
    )


def fallback_background_rows(background_df: pd.DataFrame, selection_id: int) -> list[dict]:
    if background_df.empty:
        return []
    subset = background_df[background_df["selection_id"].astype(int) == int(selection_id)].copy()
    if subset.empty:
        return []
    return subset.to_dict(orient="records")


def background_variants_for_selection(
    row: pd.Series,
    background_df: pd.DataFrame,
) -> list:
    records = parse_variant_rows_json(row.get("background_variants_json"))
    if not records:
        records = fallback_background_rows(background_df, int(row["selection_id"]))
    variants = []
    for record in records:
        variants.append(
            make_variant(
                record.get("background_chrom", record.get("susie_chrom")),
                record.get("background_pos", record.get("susie_pos")),
                record.get("background_ref", record.get("susie_ref")),
                record.get("background_alt", record.get("susie_alt")),
                record.get("background_variant_id", record.get("susie_variant_key")),
            )
        )
    return variants


def baseline_rows_for_selection(
    client,
    row: pd.Series,
    interval: genome.Interval,
    singleton_variant,
) -> pd.DataFrame:
    scores = client.score_variant(
        interval,
        singleton_variant,
        variant_scorers=[variant_scorers.GeneMaskLFCScorer(requested_output=RNA)],
        organism=ORG,
    )
    tidy = variant_scorers.tidy_scores(scores, match_gene_strand=False)
    tidy = normalize_tidy(tidy)
    tidy = filter_target_gene(tidy, str(row["gene_id"]))
    tidy = filter_skeletal_muscle_rows(tidy, require_match=True)
    if tidy.empty:
        return tidy

    selection_meta = manifest_selection_metadata(row)
    for key, value in selection_meta.items():
        tidy[key] = value
    tidy = add_track_identity(tidy, default_output_type="RNA_SEQ")
    tidy["mode"] = "baseline_score_variant"
    tidy["reducer_name"] = "GeneMaskLFCScorer.raw_score"
    tidy["scalar_score"] = tidy["raw_score"]
    tidy["scorer_friendly"] = tidy["variant_scorer"].map(map_friendly)
    return tidy


def haplotype_rows_for_selection(
    client,
    reference_reader,
    row: pd.Series,
    interval: genome.Interval,
    singleton_variant,
    background_variants,
    ontology_curies: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reference_sequence = reference_reader.fetch_interval(interval)
    background_sequence = apply_variants_to_sequence(reference_sequence, interval, background_variants)
    background_plus_singleton_sequence = apply_variants_to_sequence(
        reference_sequence,
        interval,
        [*background_variants, singleton_variant],
    )

    common_kwargs = {
        "requested_outputs": [RNA],
        "ontology_terms": ontology_curies,
        "organism": ORG,
        "interval": interval,
    }
    background_output = client.predict_sequence(background_sequence, **common_kwargs)
    alternate_output = client.predict_sequence(background_plus_singleton_sequence, **common_kwargs)
    if background_output.rna_seq is None or alternate_output.rna_seq is None:
        raise ValueError(f"RNA_SEQ output missing for selection_id={row['selection_id']}")

    selection_meta = manifest_selection_metadata(row)
    selection_meta["background_variant_count_manifest"] = int(
        row.get("background_variant_count_manifest", len(background_variants))
    )
    selection_meta["interval_str"] = f"{interval.chromosome}:{int(interval.start)}-{int(interval.end)}:."

    background_rows = build_track_rows_from_track_data(
        background_output.rna_seq,
        selection_metadata=selection_meta,
        sequence_role="background_only",
        sequence_label="background_only",
        center_variant_pos=int(row["pos"]),
    )
    alternate_rows = build_track_rows_from_track_data(
        alternate_output.rna_seq,
        selection_metadata=selection_meta,
        sequence_role="background_plus_singleton",
        sequence_label="background_plus_singleton",
        center_variant_pos=int(row["pos"]),
    )
    background_rows = filter_skeletal_muscle_rows(background_rows, require_match=True)
    alternate_rows = filter_skeletal_muscle_rows(alternate_rows, require_match=True)
    comparison_rows = reduce_predict_sequence_delta(background_rows, alternate_rows)
    return background_rows, alternate_rows, comparison_rows


def build_comparison_outputs(
    baseline_df: pd.DataFrame,
    haplotype_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    baseline_scalar = pd.DataFrame()
    if not baseline_df.empty:
        baseline_scalar = baseline_df[
            [
                "selection_id",
                "selection_group",
                "group_rank",
                "gene_id",
                "gene_symbol",
                "singleton_variant_id",
                "variant_id",
                "chrom",
                "pos",
                "ref",
                "alt",
                "track_id",
                "output_type",
                "track_name",
                "track_strand",
                "assay_title",
                "ontology_curie",
                "biosample_name",
                "biosample_type",
                "gtex_tissue",
                "mode",
                "reducer_name",
                "scalar_score",
            ]
        ].copy()
        baseline_scalar = baseline_scalar.rename(columns={"scalar_score": "baseline_scalar"})

    haplotype_scalar = pd.DataFrame()
    if not haplotype_df.empty:
        haplotype_scalar = haplotype_df[
            [
                "selection_id",
                "selection_group",
                "group_rank",
                "gene_id",
                "gene_symbol",
                "singleton_variant_id",
                "variant_id",
                "chrom",
                "pos",
                "ref",
                "alt",
                "track_id",
                "output_type",
                "track_name",
                "track_strand",
                "assay_title",
                "ontology_curie",
                "biosample_name",
                "biosample_type",
                "gtex_tissue",
                "mode",
                "reducer_name",
                "scalar_score",
                "track_mean_background_only",
                "track_mean_background_plus_singleton",
                "scalar_score_center_bin_delta",
            ]
        ].copy()
        haplotype_scalar = haplotype_scalar.rename(columns={"scalar_score": "haplotype_scalar"})

    if baseline_scalar.empty and haplotype_scalar.empty:
        empty = pd.DataFrame()
        return empty, empty

    wide = baseline_scalar.merge(
        haplotype_scalar,
        on=[
            "selection_id",
            "track_id",
        ],
        how="outer",
        suffixes=("_baseline", "_haplotype"),
    )
    for column in [
        "selection_group",
        "group_rank",
        "gene_id",
        "gene_symbol",
        "singleton_variant_id",
        "variant_id",
        "chrom",
        "pos",
        "ref",
        "alt",
        "output_type",
        "track_name",
        "track_strand",
        "assay_title",
        "ontology_curie",
        "biosample_name",
        "biosample_type",
        "gtex_tissue",
    ]:
        left = f"{column}_baseline"
        right = f"{column}_haplotype"
        if left in wide.columns and right in wide.columns:
            wide[column] = wide[left].combine_first(wide[right])
            wide = wide.drop(columns=[left, right])

    wide["haplotype_minus_baseline"] = wide["haplotype_scalar"] - wide["baseline_scalar"]

    long_frames = []
    if not baseline_scalar.empty:
        frame = baseline_scalar.rename(columns={"baseline_scalar": "scalar_score"}).copy()
        long_frames.append(frame)
    if not haplotype_scalar.empty:
        frame = haplotype_scalar.rename(columns={"haplotype_scalar": "scalar_score"}).copy()
        long_frames.append(frame)
    comparison_long = pd.concat(long_frames, ignore_index=True) if long_frames else pd.DataFrame()
    if not comparison_long.empty:
        comparison_long = comparison_long.merge(
            wide[["selection_id", "track_id", "baseline_scalar", "haplotype_scalar", "haplotype_minus_baseline"]],
            on=["selection_id", "track_id"],
            how="left",
        )

    return comparison_long, wide


def main() -> None:
    args = parse_args()
    layout = build_layout(
        dataset_id=args.dataset_id,
        sample_id=args.sample_id,
        root_dir=args.root_dir,
    )
    jobs_manifest_path = Path(args.jobs_manifest) if args.jobs_manifest else layout.singleton_haplotype_jobs_parquet
    background_manifest_path = (
        Path(args.background_manifest)
        if args.background_manifest
        else layout.singleton_haplotype_background_parquet
    )
    output_dir = Path(args.output_dir) if args.output_dir else layout.workflow_results_dir("singleton_haplotype")
    output_dir.mkdir(parents=True, exist_ok=True)

    jobs_df = read_table(jobs_manifest_path)
    require_columns(jobs_df, REQUIRED_MANIFEST_COLUMNS, "jobs manifest")
    background_df = read_table(background_manifest_path) if background_manifest_path.exists() else pd.DataFrame()

    if args.selection_ids:
        keep = set(int(selection_id) for selection_id in args.selection_ids)
        jobs_df = jobs_df[jobs_df["selection_id"].astype(int).isin(keep)].copy()
    jobs_df = jobs_df.sort_values(["selection_group", "group_rank", "selection_id"]).reset_index(drop=True)
    if args.max_selections is not None:
        jobs_df = jobs_df.head(args.max_selections).copy()

    if jobs_df.empty:
        raise ValueError("no selections to score after filtering")

    ontology_curies = args.ontology_curies or list(SKELETAL_MUSCLE_ONTOLOGY_CURIES)
    if args.mode in {"both", "haplotype"} and not args.reference_fasta:
        raise ValueError("--reference-fasta is required for haplotype scoring")

    client = load_api_client(env_file=args.env_file)
    reference_reader = None
    if args.mode in {"both", "haplotype"}:
        reference_reader = open_reference_reader(args.reference_fasta, fai_path=args.reference_fai)

    baseline_frames: list[pd.DataFrame] = []
    predict_track_frames: list[pd.DataFrame] = []
    haplotype_frames: list[pd.DataFrame] = []

    for _, row in jobs_df.iterrows():
        selection_id = int(row["selection_id"])
        print(f"Scoring selection_id={selection_id} ({row['gene_symbol']})")
        interval = interval_from_row(row)
        singleton_variant = make_variant(row["chrom"], row["pos"], row["ref"], row["alt"], row["variant_id"])

        if args.mode in {"both", "baseline"}:
            baseline_rows = baseline_rows_for_selection(client, row, interval, singleton_variant)
            if baseline_rows.empty:
                print(f"  baseline branch produced no skeletal-muscle rows for selection_id={selection_id}")
            else:
                baseline_frames.append(baseline_rows)

        if args.mode in {"both", "haplotype"}:
            background_variants = background_variants_for_selection(row, background_df)
            background_rows, alternate_rows, comparison_rows = haplotype_rows_for_selection(
                client,
                reference_reader,
                row,
                interval,
                singleton_variant,
                background_variants,
                ontology_curies,
            )
            if background_rows.empty or alternate_rows.empty:
                print(f"  haplotype branch produced no skeletal-muscle rows for selection_id={selection_id}")
            else:
                predict_track_frames.extend([background_rows, alternate_rows])
                haplotype_frames.append(comparison_rows)

    baseline_df = pd.concat(baseline_frames, ignore_index=True) if baseline_frames else pd.DataFrame()
    predict_track_df = pd.concat(predict_track_frames, ignore_index=True) if predict_track_frames else pd.DataFrame()
    haplotype_df = pd.concat(haplotype_frames, ignore_index=True) if haplotype_frames else pd.DataFrame()
    comparison_long, comparison_wide = build_comparison_outputs(baseline_df, haplotype_df)

    write_table_bundle(
        baseline_df,
        parquet_path=output_dir / "baseline_score_variant_raw.parquet",
        tsv_path=output_dir / "baseline_score_variant_raw.tsv.gz",
    )
    write_table_bundle(
        predict_track_df,
        parquet_path=output_dir / "haplotype_predict_sequence_track_outputs.parquet",
        tsv_path=output_dir / "haplotype_predict_sequence_track_outputs.tsv.gz",
    )
    write_table_bundle(
        comparison_long,
        parquet_path=output_dir / "selection_track_comparison.parquet",
        tsv_path=output_dir / "selection_track_comparison.tsv.gz",
    )
    write_table_bundle(
        comparison_wide,
        parquet_path=output_dir / "selection_track_comparison_wide.parquet",
        tsv_path=output_dir / "selection_track_comparison_wide.tsv.gz",
    )

    print(f"Wrote results under {output_dir}")


if __name__ == "__main__":
    main()
