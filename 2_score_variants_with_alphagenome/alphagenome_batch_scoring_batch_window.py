#!/usr/bin/env python3

from __future__ import annotations

import argparse
import glob
import gzip
import math
import os
import re
import shutil
import sys
import time
from pathlib import Path

import pandas as pd
from alphagenome.models import variant_scorers

REPO_ROOT = Path(__file__).resolve().parent.parent
PREPARE_RESULTS_DIR = REPO_ROOT / "3_prepare_results"
if str(PREPARE_RESULTS_DIR) not in sys.path:
    sys.path.insert(0, str(PREPARE_RESULTS_DIR))

from modules.normalisation_helper import add_protocol_track_columns

from alphagenome_shared import (
    ORG,
    RNA,
    build_layout,
    ensg_core,
    ensure_variant_ids,
    inject_into_anndata,
    load_api_client,
    load_gene_whitelist,
    make_intervals,
    make_variants,
    map_friendly,
    normalize_tidy,
    resolve_legacy_input_paths,
)

DEFAULT_BATCH = 128
DEFAULT_SEQ_LEN = 2**20
ENS_RE = re.compile(r"(ENSG\d+)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Legacy AlphaGenome batch-window scorer using score_variants()",
    )
    parser.add_argument("--dataset-id", help="dataset id; falls back to DATASET_ID")
    parser.add_argument("--sample-id", help="sample id; falls back to SAMPLE_ID")
    parser.add_argument("--root-dir", help="dataset root; falls back to ROOT_DIR/PDC_TMP")
    parser.add_argument("--variants", help="variant TSV(.gz); defaults to dataset-centric inputs")
    parser.add_argument("--targets-bed", help="TSS target BED; defaults to dataset-centric inputs")
    parser.add_argument("--gene-list", help="optional gene whitelist; defaults to dataset-centric inputs")
    parser.add_argument(
        "--env-file",
        default=Path(__file__).resolve().with_name(".env"),
        help="dotenv file containing API_KEY_PERSONAL",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=int(os.getenv("BATCH", DEFAULT_BATCH)),
        help=f"variants per batch [default: {DEFAULT_BATCH}]",
    )
    parser.add_argument(
        "--seq-len-gm",
        type=int,
        default=int(os.getenv("SEQ_LEN_GM", DEFAULT_SEQ_LEN)),
        help=f"legacy GeneMask sequence length [default: {DEFAULT_SEQ_LEN}]",
    )
    parser.add_argument(
        "--stitch",
        action="store_true",
        default=os.getenv("STITCH", "0") == "1",
        help="stitch chunk_*.tsv.gz files into a single annotated TSV",
    )
    return parser.parse_args()


def parse_existing_chunks(chunks_dir: Path) -> set[tuple[int, int]]:
    done = set()
    for path in chunks_dir.glob("chunk_*.tsv.gz"):
        start, end = path.name.replace("chunk_", "").replace(".tsv.gz", "").split("_")
        done.add((int(start), int(end)))
    return done


def env_int(name: str, default: str) -> int:
    raw = os.getenv(name, default)
    try:
        return int(raw)
    except ValueError as exc:
        raise ValueError(f"{name} must be an integer, got {raw!r}") from exc


def job_batch_span(num_batches: int, job_index: int, job_total: int) -> tuple[int, int]:
    if job_total < 1:
        raise ValueError("JOB_TOTAL must be >= 1")
    if job_index < 0 or job_index >= job_total:
        raise ValueError(f"JOB_INDEX {job_index} invalid for total {job_total}")

    base = num_batches // job_total
    remainder = num_batches % job_total
    start = job_index * base + min(job_index, remainder)
    extra = 1 if job_index < remainder else 0
    end = start + base + extra
    return start, end


def default_job_bounds(
    num_variants: int,
    batch_size: int,
    job_index: int,
    job_total: int,
) -> tuple[int, int]:
    if num_variants <= 0:
        return 0, 0
    num_batches = math.ceil(num_variants / batch_size)
    start_batch, end_batch = job_batch_span(num_batches, job_index, job_total)
    return start_batch * batch_size, min(end_batch * batch_size, num_variants)


def load_legacy_variant_table(variants_path: Path) -> pd.DataFrame:
    df = pd.read_csv(variants_path, sep="\t", compression="infer", low_memory=False)
    required = {"CHROM", "POS", "REF", "ALT", "gene_tag"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"missing columns in {variants_path}: {missing}")

    if "variant_id" not in df.columns:
        df["variant_id"] = (
            df["CHROM"].astype(str)
            + ":"
            + pd.to_numeric(df["POS"], errors="coerce").astype("Int64").astype(str)
            + ":"
            + df["REF"].astype(str)
            + ">"
            + df["ALT"].astype(str)
        )

    out = df.loc[:, ["variant_id", "gene_tag", "CHROM", "POS", "REF", "ALT"]].copy()
    out["POS"] = pd.to_numeric(out["POS"], errors="coerce")
    out = out.dropna(subset=["POS"]).copy()
    out["POS"] = out["POS"].astype(int)
    out["gene_tag_core"] = out["gene_tag"].map(ensg_core)
    if "is_anchor" not in out.columns:
        out["is_anchor"] = False
    return out


def load_targets_bed(targets_path: Path) -> pd.DataFrame:
    bed = pd.read_csv(
        targets_path,
        sep="\t",
        header=None,
        names=["chr", "start", "end", "name"],
        dtype={"chr": str, "start": int, "end": int, "name": str},
    )
    bed["chr"] = bed["chr"].astype(str).map(lambda chrom: chrom if chrom.startswith("chr") else f"chr{chrom}")
    bed["tss_pos"] = ((bed["start"].astype(int) + bed["end"].astype(int)) // 2).astype(int)
    bed["gene_tag"] = bed["name"].astype(str)
    bed["ensg_core"] = bed["name"].map(ensg_core)
    return bed


def add_anchor_rows(df: pd.DataFrame, bed: pd.DataFrame) -> pd.DataFrame:
    have_genes = set(df["gene_tag_core"].dropna().astype(str).unique())
    all_genes = set(bed["ensg_core"].dropna().astype(str).unique())
    missing_core = sorted(all_genes - have_genes)

    anchors = pd.DataFrame(columns=df.columns)
    if missing_core:
        bed_missing = bed[bed["ensg_core"].isin(missing_core)].copy()
        bed_missing["REF"] = "A"
        bed_missing["ALT"] = "C"

        clash_keys = set(
            zip(
                df["CHROM"].astype(str),
                df["POS"].astype(int),
                df["gene_tag_core"].astype(str),
            )
        )

        def maybe_flip(row: pd.Series) -> str:
            key = (str(row["chr"]), int(row["tss_pos"]), str(row["ensg_core"]))
            return "G" if key in clash_keys else "C"

        bed_missing["ALT"] = bed_missing.apply(maybe_flip, axis=1)
        anchors = pd.DataFrame(
            {
                "variant_id": (
                    bed_missing["chr"].astype(str)
                    + ":"
                    + bed_missing["tss_pos"].astype(int).astype(str)
                    + ":A>"
                    + bed_missing["ALT"].astype(str)
                ),
                "gene_tag": bed_missing["ensg_core"].astype(str),
                "gene_tag_core": bed_missing["ensg_core"].astype(str),
                "CHROM": bed_missing["chr"].astype(str),
                "POS": bed_missing["tss_pos"].astype(int),
                "REF": "A",
                "ALT": bed_missing["ALT"].astype(str),
                "is_anchor": True,
            }
        )

    out = df.copy()
    out["is_anchor"] = out["is_anchor"].fillna(False).astype(bool)
    if len(anchors):
        out = pd.concat([out, anchors.loc[:, out.columns]], ignore_index=True)
        out = out.drop_duplicates(subset=["variant_id", "gene_tag"], keep="first")
    return out


def filter_muscle_legacy(df: pd.DataFrame) -> pd.DataFrame:
    mask_muscle = pd.Series(False, index=df.index)
    if "gtex_tissue" in df.columns:
        norm = df["gtex_tissue"].astype(str).map(
            lambda value: re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")
        )
        mask_muscle = norm.str.contains("muscle", na=False) & norm.str.contains("skeletal", na=False)
    mask_uberon = pd.Series(False, index=df.index)
    if "ontology_curie" in df.columns:
        mask_uberon = df["ontology_curie"].astype(str).str.contains("UBERON:0001134", case=False, na=False)
    if mask_muscle.any() or mask_uberon.any():
        return df[mask_muscle | mask_uberon].copy()
    return df


def stitch_chunks(chunks_dir: Path, output_path: Path) -> None:
    chunk_files = sorted(glob.glob((chunks_dir / "chunk_*.tsv.gz").as_posix()))
    if not chunk_files:
        raise FileNotFoundError(f"no chunk files found under {chunks_dir}")

    tmp_final = output_path.with_suffix(output_path.suffix + ".tmp")
    with gzip.open(tmp_final, "wt") as writer:
        wrote_header = False
        for chunk_path in chunk_files:
            with gzip.open(chunk_path, "rt") as reader:
                if not wrote_header:
                    shutil.copyfileobj(reader, writer)
                    wrote_header = True
                else:
                    reader.readline()
                    shutil.copyfileobj(reader, writer)
    os.replace(tmp_final, output_path)


def main() -> None:
    args = parse_args()
    layout = build_layout(
        dataset_id=args.dataset_id,
        sample_id=args.sample_id,
        root_dir=args.root_dir,
    )
    variants_path, targets_path, gene_list_path = resolve_legacy_input_paths(
        layout,
        variants_path=args.variants,
        targets_path=args.targets_bed,
        gene_list_path=args.gene_list,
    )

    model = load_api_client(env_file=args.env_file)
    output_path = layout.results_dir / f"{layout.sample_id}.variants.annotated.tsv.gz"
    chunks_dir = layout.chunks_dir
    chunks_dir.mkdir(parents=True, exist_ok=True)

    df = load_legacy_variant_table(variants_path)
    bed = load_targets_bed(targets_path)
    df = add_anchor_rows(df, bed)
    gene_whitelist = load_gene_whitelist(gene_list_path)

    total_variants = len(df)
    num_batches = math.ceil(total_variants / args.batch_size)
    done = parse_existing_chunks(chunks_dir)

    job_total = env_int("JOB_TOTAL", "1")
    job_index = env_int("JOB_INDEX", "0")
    default_raw_start, default_raw_end = default_job_bounds(
        total_variants,
        args.batch_size,
        job_index,
        job_total,
    )
    raw_start = int(os.getenv("RAW_START", str(default_raw_start)))
    raw_end = int(os.getenv("RAW_END", str(default_raw_end)))

    slice_start = (raw_start // args.batch_size) * args.batch_size
    slice_end = ((raw_end + args.batch_size - 1) // args.batch_size) * args.batch_size
    print(
        f"job {job_index + 1}/{job_total} raw slice {raw_start}-{raw_end} "
        f"aligned to {slice_start}-{slice_end} of {total_variants} variants"
    )

    all_batch_ids = list(range(num_batches - 1, -1, -1))
    batch_ids = [
        batch_id
        for batch_id in all_batch_ids
        if (batch_id * args.batch_size) >= slice_start and (batch_id * args.batch_size) < slice_end
    ]
    if not batch_ids:
        print(f"no batches assigned for job {job_index + 1}/{job_total}, nothing to process")

    for batch_id in batch_ids:
        start = batch_id * args.batch_size
        end = min(start + args.batch_size, total_variants)
        out_chunk = chunks_dir / f"chunk_{start:07d}_{end:07d}.tsv.gz"
        lock_path = Path(f"{out_chunk}.lock")
        tmp_chunk = Path(f"{out_chunk}.tmp")

        if (start, end) in done or out_chunk.exists():
            print(f"Skipping existing chunk {start}-{end}")
            continue

        try:
            fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            print(f"Locked by another process, skip {start}-{end}")
            continue

        try:
            print(f"start batch {start}-{end} size {end - start}")
            meta = df.iloc[start:end].copy()
            intervals = make_intervals(args.seq_len_gm, meta["CHROM"].to_numpy(), meta["POS"].to_numpy())
            variants = make_variants(meta["CHROM"], meta["POS"], meta["REF"], meta["ALT"], meta["variant_id"])
            interval_to_varid = {
                f"{interval.chromosome}:{int(interval.start)}-{int(interval.end)}:.": variant_id
                for interval, variant_id in zip(intervals, meta["variant_id"])
            }

            t0 = time.time()
            scores = model.score_variants(
                intervals=intervals,
                variants=variants,
                variant_scorers=[variant_scorers.GeneMaskLFCScorer(requested_output=RNA)],
                organism=ORG,
                progress_bar=False,
            )
            print(f"score_variants done in {time.time() - t0:.1f}s")
            inject_into_anndata(scores, meta)

            tidy = variant_scorers.tidy_scores(scores, match_gene_strand=False)
            tidy = normalize_tidy(tidy)
            tidy = add_protocol_track_columns(tidy)
            ensure_variant_ids(tidy, interval_to_varid)

            anchor_map = dict(zip(meta["variant_id"], meta["is_anchor"]))
            tidy["is_anchor"] = tidy["variant_id"].map(anchor_map).fillna(False)
            tidy["seq_len"] = args.seq_len_gm
            tidy["scorer_friendly"] = tidy["variant_scorer"].map(map_friendly)

            out_df = tidy
            if gene_whitelist:
                out_df = out_df.copy()
                out_df["gene_id_core"] = out_df["gene_id"].map(ensg_core)
                out_df = out_df[out_df["gene_id_core"].isin(gene_whitelist)].copy()
                out_df = out_df.drop(columns=["gene_id_core"], errors="ignore")
                if out_df.empty:
                    print(f"no whitelist genes in chunk {start}-{end}, skip write")
                    continue

            out_df = filter_muscle_legacy(out_df)
            out_df.to_csv(tmp_chunk, sep="\t", index=False, compression="gzip")
            os.replace(tmp_chunk, out_chunk)
            print(f"wrote {out_chunk.name}, progress {end}/{total_variants}")
        finally:
            try:
                os.close(fd)
            except Exception:
                pass
            try:
                os.remove(lock_path)
            except FileNotFoundError:
                pass

    if args.stitch:
        stitch_chunks(chunks_dir, output_path)
        print(f"Done. Wrote {output_path}")


if __name__ == "__main__":
    main()
