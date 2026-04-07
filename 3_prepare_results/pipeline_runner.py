from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import polars as pl

from modules.aggregator import aggregate_genes
from modules.permuted_af import add_perm_af_gene_aware, load_gene_af_pools
from modules.sanity_check import compare_vg_medians
from modules.snv_filter import filter_to_snvs
from modules.synthetic_variant_downsampler import (
    downsample_to_real_counts,
    load_real_variant_counts,
)
from modules.variant_deduplicator import deduplicate_by_gene_and_variant


def _timed(msg: str) -> float:
    print(msg, flush=True)
    return time.perf_counter()


def _elapsed(start: float) -> str:
    return f"{time.perf_counter() - start:.2f}s"


def _load_processing_columns(variant_path: Path, filter_snvs: bool) -> pl.DataFrame:
    if not variant_path.exists():
        raise FileNotFoundError(f"Input variants parquet not found: {variant_path}")

    lf = pl.scan_parquet(variant_path)
    schema_cols = set(lf.collect_schema().names())

    required = ["gene_id", "variant_id"]
    missing = [c for c in required if c not in schema_cols]
    if missing:
        raise ValueError(f"Input parquet missing required columns: {missing}")

    if filter_snvs:
        required_snv = [c for c in ["REF", "ALT"] if c not in schema_cols]
        if required_snv:
            raise ValueError(f"--filter-snvs requires columns: {required_snv}")

    return pl.read_parquet(variant_path, low_memory=True)


def _resolve_base_root(root_dir: Path | None) -> Path:
    if root_dir:
        return root_dir
    return Path(os.getenv("ROOT_DIR") or os.getenv("PDC_TMP") or "/cfs/klemming/scratch/m/mmarandi")


def _process_variants(args: argparse.Namespace, variant_path: Path, is_ism: bool) -> Path:
    t0 = _timed(f"Loading variants from {variant_path}...")
    df = _load_processing_columns(variant_path, filter_snvs=args.filter_snvs)
    print(f"  Loaded {len(df):,} rows, {len(df.columns)} columns in {_elapsed(t0)}", flush=True)

    if args.filter_snvs:
        t = _timed("Filtering to SNVs only (removing indels and MNVs)...")
        df = filter_to_snvs(df, verbose=True)
        print(f"  SNV filtering completed in {_elapsed(t)}", flush=True)

    if args.gene_list:
        t = _timed(f"Filtering variants using whitelist: {args.gene_list}")
        whitelist_df = pl.read_csv(
            args.gene_list,
            has_header=False,
            new_columns=["gene_id"],
            separator="\t",
        )
        whitelist = whitelist_df["gene_id"].to_list()
        df = df.filter(pl.col("gene_id").is_in(whitelist))
        print(f"  Rows after gene list filtering: {len(df):,} ({_elapsed(t)})", flush=True)

    t = _timed("Deduplicating by gene_id, variant_id, and track_key...")
    df = deduplicate_by_gene_and_variant(df, verbose=True)
    print(f"  Deduplication completed in {_elapsed(t)}", flush=True)

    if is_ism:
        if not args.real_reference:
            raise ValueError("Flag --real-reference is required for ISM/NULL processing")

        real_ref = args.real_reference.resolve()
        if not real_ref.exists():
            raise FileNotFoundError(f"Real reference not found: {real_ref}")

        t = _timed(f"Loading real variant counts from {real_ref}...")
        real_counts = load_real_variant_counts(real_ref, filter_snvs=args.filter_snvs)
        print(f"  Loaded real variant counts in {_elapsed(t)}", flush=True)

        t = _timed("Downsampling to match real counts...")
        df = downsample_to_real_counts(df, real_counts, seed=42, verbose=True)
        print(f"  Downsampling completed in {_elapsed(t)}", flush=True)
        suffix = "_downsampled"
    else:
        suffix = "_dedup"

    if args.filter_snvs:
        suffix += "_snv"

    permute_requested = args.permute_af or args.permute_af_sanity
    if permute_requested:
        _timed("--- Starting AF Permutation ---")
        perm_ref_source: Path | pl.DataFrame | None = None

        if args.real_reference:
            perm_ref_source = args.real_reference.resolve()
        elif not is_ism:
            print("  No external reference provided for real data; using processed variants as AF source")
            perm_ref_source = df

        if perm_ref_source is None:
            raise ValueError("Cannot permute AF: No reference provided (use --real-reference)")

        if isinstance(perm_ref_source, Path) and not perm_ref_source.exists():
            raise FileNotFoundError(f"Permutation reference not found: {perm_ref_source}")

        t = _timed(f"  Loading AF pools from {perm_ref_source}...")
        af_pools = load_gene_af_pools(perm_ref_source)
        print(f"  Loaded AF pools in {_elapsed(t)}", flush=True)

        t = _timed("  Sampling and assigning perm_AF...")
        df = add_perm_af_gene_aware(df, af_pools, seed=42, strict_no_replacement=True)
        print(f"  perm_AF assignment completed in {_elapsed(t)}", flush=True)

        if args.permute_af_sanity and not is_ism:
            suffix += "_perm_sanity"
        else:
            suffix += "_perm"

    if args.variant_out:
        processed_path = args.variant_out.resolve()
    else:
        stem = variant_path.stem
        if not stem.endswith(suffix):
            stem += suffix
        processed_path = variant_path.with_name(f"{stem}.parquet")

    t = _timed(f"Writing processed variants to {processed_path}...")
    df.write_parquet(processed_path)
    print(f"  Wrote processed variants in {_elapsed(t)}", flush=True)
    return processed_path


def main() -> None:
    parser = argparse.ArgumentParser(description="process and aggregate variant scores")
    parser.add_argument("--variants-parquet", type=Path, required=True, help="input variants parquet")
    parser.add_argument("--gene-list", type=Path, help="optional gene whitelist")
    parser.add_argument("--variant-out", type=Path, help="output variants parquet")
    parser.add_argument("--gene-out", type=Path, required=True, help="output genes parquet")
    parser.add_argument("--root-dir", type=Path, help="root directory for reference data")
    parser.add_argument("--real-reference", type=Path, help="real dataset reference for downsampling and CI")
    parser.add_argument("--deduplicate", action="store_true", help="deduplicate variants")
    parser.add_argument("--permute-af", action="store_true", help="generate perm_AF")
    parser.add_argument("--permute-af-sanity", action="store_true", help="generate perm_AF for sanity checks in real datasets")
    parser.add_argument("--sanity-report-out", type=Path, help="optional JSON output path for vg median sanity check")
    parser.add_argument("--calc-ci", action="store_true", help="calculate confidence intervals")
    parser.add_argument("--synthetic", action="store_true", help="synthetic/null mode: use only perm_AF for all Vg metrics")
    parser.add_argument("--filter-snvs", action="store_true", help="filter to keep only SNVs (removes indels and MNVs)")
    parser.add_argument(
        "--global-unique-variant-ids",
        action="store_true",
        help="apply protocol-aware unique-by-gene_id, variant_id, and track_key during aggregation",
    )

    args = parser.parse_args()

    base_root = _resolve_base_root(args.root_dir)
    variant_path = args.variants_parquet.resolve()

    is_ism = "ism" in str(variant_path).lower() or "null" in str(variant_path).lower()

    if args.deduplicate:
        variant_path = _process_variants(args, variant_path, is_ism)

    gene_path = args.gene_out.resolve()
    gene_list = args.gene_list.resolve() if args.gene_list else None

    ci_ref = None
    if args.real_reference:
        ci_ref = args.real_reference.resolve()
    elif not is_ism:
        ci_ref = variant_path

    print(f"Aggregating genes from {variant_path}...", flush=True)
    synthetic_mode = args.synthetic or is_ism
    if synthetic_mode and not args.synthetic:
        print("Detected ISM/NULL dataset; enabling synthetic aggregation mode automatically.")

    t_agg = time.perf_counter()
    aggregate_genes(
        variant_path,
        gene_path,
        base_ref=base_root,
        is_ism=is_ism,
        gene_list_path=gene_list,
        calculate_ci=args.calc_ci,
        real_reference_path=ci_ref,
        n_permutations=1000,
        is_synthetic=synthetic_mode,
        include_perm_sanity=(args.permute_af_sanity and not synthetic_mode),
        global_unique_variant_ids=args.global_unique_variant_ids,
    )
    print(f"Gene aggregation finished in {_elapsed(t_agg)}", flush=True)

    if args.permute_af_sanity and not synthetic_mode:
        report_path = args.sanity_report_out
        if report_path is None:
            report_path = gene_path.with_name(f"{gene_path.stem}_vg_median_sanity.json")

        t_sanity = _timed(f"Running Vg median sanity check on {gene_path}...")
        sanity = compare_vg_medians(gene_path, out_path=report_path)
        print(
            "  Sanity medians: "
            f"obs={sanity['median_observed']:.6g}, "
            f"perm_sanity={sanity['median_sanity']:.6g}, "
            f"ratio={sanity['median_ratio_obs_over_sanity']:.6g}, "
            f"delta={sanity['median_delta']:.6g}",
            flush=True,
        )
        print(f"  Sanity report saved to {report_path} ({_elapsed(t_sanity)})", flush=True)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()