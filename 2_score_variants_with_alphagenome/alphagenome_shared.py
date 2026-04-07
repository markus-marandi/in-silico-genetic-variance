from __future__ import annotations

import json
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from dotenv import load_dotenv


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from helpers.path_manager import DEFAULT_ROOT, ProjectLayout  # noqa: E402

from alphagenome.data import genome  # noqa: E402
from alphagenome.data import track_data as track_data_lib  # noqa: E402
from alphagenome.models import dna_client, variant_scorers  # noqa: E402
from alphagenome.models import dna_client as _dc  # noqa: E402


DEFAULT_SEQ_LEN_1MB = 2**20
RNA = dna_client.OutputType.RNA_SEQ
ORG = _dc.Organism.HOMO_SAPIENS
DEFAULT_API_KEY_ENV = "API_KEY_PERSONAL"
SKELETAL_MUSCLE_ONTOLOGY_CURIES = (
    "UBERON:0011907",
    "UBERON:0001134",
    "CL:0000515",
    "CL:0000594",
)
TRACK_ID_COLUMNS = (
    "output_type",
    "track_name",
    "track_strand",
    "assay_title",
    "ontology_curie",
    "biosample_name",
    "biosample_type",
    "gtex_tissue",
)
RENAME_FRIENDLY = {
    "GeneMaskLFCScorer": "gene_exonmask_delta_log2",
    "CenterMaskScorer(width=10001).DIFF_MEAN": "center10.001kb_diff_mean",
    "CenterMaskScorer(width=100001).DIFF_MEAN": "center100.001kb_diff_mean",
    "CenterMaskScorer(width=10001).L2_DIFF_LOG1P": "center10.001kb_l2_log1p",
    "CenterMaskScorer(width=100001).L2_DIFF_LOG1P": "center100.001kb_l2_log1p",
}
ENS_RE = re.compile(r"(ENSG\d+)")


def build_layout(
    *,
    dataset_id: str | None = None,
    sample_id: str | None = None,
    root_dir: str | Path | None = None,
) -> ProjectLayout:
    dataset_id = dataset_id or os.getenv("DATASET_ID")
    sample_id = sample_id or os.getenv("SAMPLE_ID")
    if not dataset_id or not sample_id:
        raise ValueError("dataset_id/sample_id or DATASET_ID/SAMPLE_ID are required")

    resolved_root = (
        Path(root_dir)
        if root_dir is not None
        else Path(os.getenv("ROOT_DIR", os.getenv("PDC_TMP", str(DEFAULT_ROOT))))
    )
    layout = ProjectLayout(
        dataset_id=dataset_id,
        sample_id=sample_id,
        root_dir=resolved_root,
    )
    layout.make_dirs()
    return layout


def resolve_legacy_input_paths(
    layout: ProjectLayout,
    *,
    variants_path: str | Path | None = None,
    targets_path: str | Path | None = None,
    gene_list_path: str | Path | None = None,
) -> tuple[Path, Path, Path | None]:
    variants = Path(variants_path) if variants_path is not None else layout.resolve_input_variants()
    targets = Path(targets_path) if targets_path is not None else layout.resolve_input_targets()
    if gene_list_path is None:
        gene_list = layout.resolve_input_gene_list(required=False)
    else:
        gene_list = Path(gene_list_path)
    return variants, targets, gene_list


def load_api_client(
    *,
    env_file: str | Path | None = None,
    api_key_env: str = DEFAULT_API_KEY_ENV,
):
    if env_file is not None:
        load_dotenv(env_file)
    else:
        load_dotenv()
    api_key = os.getenv(api_key_env)
    if not api_key:
        raise ValueError(f"set {api_key_env} in your environment or .env file")
    return dna_client.create(api_key=api_key)


def read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    suffixes = path.suffixes
    if suffixes[-1:] == [".parquet"] or suffixes[-2:] == [".parquet", ".gz"]:
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def write_table_bundle(
    df: pd.DataFrame,
    *,
    parquet_path: str | Path | None = None,
    tsv_path: str | Path | None = None,
) -> None:
    if parquet_path is not None:
        parquet_path = Path(parquet_path)
        parquet_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(parquet_path, index=False)
    if tsv_path is not None:
        tsv_path = Path(tsv_path)
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(tsv_path, sep="\t", index=False, compression="gzip")


def _interval_to_str(x: Any) -> str:
    if isinstance(x, str):
        return x
    if hasattr(x, "chromosome") and hasattr(x, "start") and hasattr(x, "end"):
        return f"{x.chromosome}:{int(x.start)}-{int(x.end)}:."
    return str(x)


def map_friendly(name: str) -> str:
    text = str(name)
    for needle, replacement in RENAME_FRIENDLY.items():
        if needle in text:
            return replacement
    return text


def ensg_core(value: Any) -> str | None:
    if value is None or pd.isna(value):
        return None
    match = ENS_RE.search(str(value))
    return match.group(1) if match else None


def normalize_label(value: Any) -> str:
    value = "" if value is None else str(value)
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def normalize_tidy(tidy: pd.DataFrame) -> pd.DataFrame:
    tidy = tidy.copy()
    if "variant_id" in tidy.columns:
        tidy["variant_id"] = tidy["variant_id"].map(
            lambda value: value.name if hasattr(value, "name") else str(value)
        )
    elif "variant" in tidy.columns:
        tidy["variant_id"] = tidy["variant"].map(
            lambda value: value.name if hasattr(value, "name") else str(value)
        )
    if "scored_interval" in tidy.columns:
        tidy["scored_interval_str"] = tidy["scored_interval"].map(_interval_to_str)
    else:
        tidy["scored_interval_str"] = pd.NA
    return tidy


def load_gene_whitelist(path: str | Path | None) -> set[str]:
    if not path:
        return set()
    genes = pd.read_csv(path, header=None, names=["gene_id"], dtype=str)
    genes["gene_id"] = genes["gene_id"].astype(str).str.strip()
    genes = genes[genes["gene_id"] != ""]
    genes["gene_id_core"] = genes["gene_id"].map(ensg_core)
    return set(genes["gene_id_core"].dropna())


def inject_into_anndata(scores, meta_df: pd.DataFrame) -> None:
    for per_variant, (_, row) in zip(scores, meta_df.iterrows()):
        if not isinstance(per_variant, (list, tuple)):
            per_variant = [per_variant]
        for anndata_obj in per_variant:
            anndata_obj.obs["variant_id"] = str(row.variant_id)
            anndata_obj.obs["CHROM"] = str(row.CHROM)
            anndata_obj.obs["POS"] = int(row.POS)
            anndata_obj.obs["REF"] = str(row.REF)
            anndata_obj.obs["ALT"] = str(row.ALT)
            anndata_obj.obs["gene_tag"] = str(row.gene_tag) if pd.notna(row.gene_tag) else ""


def ensure_variant_ids(df: pd.DataFrame, interval_to_varid: dict[str, str]) -> None:
    if "variant_id" not in df.columns:
        df["variant_id"] = pd.NA
    mask = df["variant_id"].isna() | (df["variant_id"].astype(str).str.strip() == "")
    if mask.any():
        df.loc[mask, "variant_id"] = df.loc[mask, "scored_interval_str"].map(interval_to_varid)


def make_interval(seq_len: int, chrom: Any, pos: Any) -> genome.Interval:
    half = seq_len // 2
    start = max(int(pos) - half, 1)
    end = start + seq_len
    return genome.Interval(chromosome=str(chrom), start=int(start), end=int(end))


def make_intervals(seq_len: int, chroms, poses) -> list[genome.Interval]:
    return [make_interval(seq_len, chrom, pos) for chrom, pos in zip(chroms, poses)]


def make_variant(
    chrom: Any,
    pos: Any,
    ref: Any,
    alt: Any,
    name: Any,
    *,
    info: dict[str, Any] | None = None,
) -> genome.Variant:
    return genome.Variant(
        chromosome=str(chrom),
        position=int(pos),
        reference_bases=str(ref),
        alternate_bases=str(alt),
        name=str(name),
        info=info or {},
    )


def make_variants(chroms, poses, refs, alts, names) -> list[genome.Variant]:
    return [
        make_variant(chrom, pos, ref, alt, name)
        for chrom, pos, ref, alt, name in zip(chroms, poses, refs, alts, names)
    ]


def normalize_track_columns(df: pd.DataFrame, *, default_output_type: str | None = None) -> pd.DataFrame:
    df = df.copy()
    rename_map = {}
    if "name" in df.columns and "track_name" not in df.columns:
        rename_map["name"] = "track_name"
    if "strand" in df.columns and "track_strand" not in df.columns:
        rename_map["strand"] = "track_strand"
    if "Assay title" in df.columns and "assay_title" not in df.columns:
        rename_map["Assay title"] = "assay_title"
    if rename_map:
        df = df.rename(columns=rename_map)

    for column in TRACK_ID_COLUMNS:
        if column not in df.columns:
            if column == "output_type" and default_output_type is not None:
                df[column] = default_output_type
            else:
                df[column] = pd.NA

    return df


def add_track_identity(df: pd.DataFrame, *, default_output_type: str | None = None) -> pd.DataFrame:
    df = normalize_track_columns(df, default_output_type=default_output_type)
    parts = []
    for column in TRACK_ID_COLUMNS:
        parts.append(df[column].astype("string").fillna(""))
    df["track_id"] = parts[0]
    for series in parts[1:]:
        df["track_id"] = df["track_id"] + "||" + series
    return df


def filter_skeletal_muscle_rows(df: pd.DataFrame, *, require_match: bool = False) -> pd.DataFrame:
    df = normalize_track_columns(df)
    ontology = df["ontology_curie"].astype("string").fillna("")
    gtex = df["gtex_tissue"].astype("string").fillna("").map(normalize_label)
    biosample = df["biosample_name"].astype("string").fillna("").map(normalize_label)

    mask = ontology.isin(SKELETAL_MUSCLE_ONTOLOGY_CURIES)
    mask = mask | ((gtex.str.contains("muscle", na=False)) & (gtex.str.contains("skeletal", na=False)))
    mask = mask | ((biosample.str.contains("skeletal", na=False)) & (biosample.str.contains("muscle", na=False)))

    if require_match or mask.any():
        return df.loc[mask].copy()
    return df


def filter_target_gene(df: pd.DataFrame, gene_id: str | None) -> pd.DataFrame:
    if gene_id is None or "gene_id" not in df.columns:
        return df
    gene_core = ensg_core(gene_id)
    if gene_core is None:
        return df
    out = df.copy()
    out["gene_id_core"] = out["gene_id"].map(ensg_core)
    out = out[out["gene_id_core"] == gene_core].copy()
    return out.drop(columns=["gene_id_core"], errors="ignore")


@dataclass(frozen=True)
class FastaIndexRecord:
    length: int
    offset: int
    line_bases: int
    line_width: int


class IndexedFastaReader:
    def __init__(self, fasta_path: str | Path, fai_path: str | Path | None = None):
        self.fasta_path = Path(fasta_path)
        self.fai_path = Path(fai_path) if fai_path is not None else Path(f"{self.fasta_path}.fai")
        if not self.fasta_path.exists():
            raise FileNotFoundError(f"reference FASTA not found: {self.fasta_path}")
        if not self.fai_path.exists():
            raise FileNotFoundError(f"reference FASTA index not found: {self.fai_path}")
        if self.fasta_path.suffix == ".gz":
            raise ValueError(
                "compressed FASTA requires pysam/pyfaidx or an uncompressed FASTA + .fai; "
                f"got {self.fasta_path}"
            )
        self.index = self._load_index()

    def _load_index(self) -> dict[str, FastaIndexRecord]:
        entries: dict[str, FastaIndexRecord] = {}
        with self.fai_path.open() as handle:
            for line in handle:
                chrom, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")[:5]
                entries[chrom] = FastaIndexRecord(
                    length=int(length),
                    offset=int(offset),
                    line_bases=int(line_bases),
                    line_width=int(line_width),
                )
        return entries

    def _resolve_chrom(self, chrom: str) -> str:
        if chrom in self.index:
            return chrom
        if chrom.startswith("chr") and chrom[3:] in self.index:
            return chrom[3:]
        if f"chr{chrom}" in self.index:
            return f"chr{chrom}"
        raise KeyError(f"chromosome {chrom!r} not found in FASTA index")

    def fetch(self, chrom: str, start: int, end: int) -> str:
        resolved_chrom = self._resolve_chrom(str(chrom))
        record = self.index[resolved_chrom]
        if start < 0 or end < start or end > record.length:
            raise ValueError(
                f"invalid interval {resolved_chrom}:{start}-{end} for contig length {record.length}"
            )

        sequence_chunks: list[str] = []
        current = start
        with self.fasta_path.open("rb") as handle:
            while current < end:
                line_index = current // record.line_bases
                line_offset = current % record.line_bases
                take = min(record.line_bases - line_offset, end - current)
                byte_offset = record.offset + (line_index * record.line_width) + line_offset
                handle.seek(byte_offset)
                sequence_chunks.append(handle.read(take).decode("ascii"))
                current += take
        return "".join(sequence_chunks).upper()

    def fetch_interval(self, interval: genome.Interval) -> str:
        return self.fetch(interval.chromosome, int(interval.start), int(interval.end))


def open_reference_reader(
    fasta_path: str | Path,
    *,
    fai_path: str | Path | None = None,
):
    try:
        import pysam  # type: ignore
    except ModuleNotFoundError:
        return IndexedFastaReader(fasta_path, fai_path=fai_path)

    class PysamReferenceReader:
        def __init__(self, path: str | Path):
            self._path = Path(path)
            self._fasta = pysam.FastaFile(self._path.as_posix())

        def fetch_interval(self, interval: genome.Interval) -> str:
            chrom = interval.chromosome
            if chrom not in self._fasta.references:
                if chrom.startswith("chr") and chrom[3:] in self._fasta.references:
                    chrom = chrom[3:]
                elif f"chr{chrom}" in self._fasta.references:
                    chrom = f"chr{chrom}"
            return self._fasta.fetch(chrom, int(interval.start), int(interval.end)).upper()

    return PysamReferenceReader(fasta_path)


def variant_span0(variant: genome.Variant) -> tuple[int, int]:
    start = int(variant.position) - 1
    end = start + len(str(variant.reference_bases))
    return start, end


def variants_overlap(left: genome.Variant, right: genome.Variant) -> bool:
    left_start, left_end = variant_span0(left)
    right_start, right_end = variant_span0(right)
    return left_start < right_end and right_start < left_end


def apply_variants_to_sequence(
    reference_sequence: str,
    interval: genome.Interval,
    variants: Sequence[genome.Variant],
) -> str:
    if not variants:
        return reference_sequence.upper()

    overlap_order = sorted(
        variants,
        key=lambda variant: (int(variant.position), len(str(variant.reference_bases))),
    )
    for first, second in zip(overlap_order, overlap_order[1:]):
        if variants_overlap(first, second):
            raise ValueError(f"overlapping variants are not supported: {first.name} and {second.name}")

    sequence = reference_sequence.upper()
    sorted_variants = sorted(
        variants,
        key=lambda variant: (int(variant.position), len(str(variant.reference_bases))),
        reverse=True,
    )

    for variant in sorted_variants:
        rel_start = (int(variant.position) - 1) - int(interval.start)
        rel_end = rel_start + len(str(variant.reference_bases))
        if rel_start < 0 or rel_end > len(sequence):
            raise ValueError(
                f"variant {variant.name} falls outside interval {_interval_to_str(interval)}"
            )
        observed_ref = sequence[rel_start:rel_end].upper()
        expected_ref = str(variant.reference_bases).upper()
        if observed_ref != expected_ref:
            raise ValueError(
                f"reference mismatch for {variant.name}: expected {expected_ref}, observed {observed_ref}"
            )
        sequence = sequence[:rel_start] + str(variant.alternate_bases).upper() + sequence[rel_end:]
    return sequence


def _json_default(value: Any) -> Any:
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if pd.isna(value):
        return None
    return value


def rows_to_json(rows: Sequence[Mapping[str, Any]]) -> str:
    return json.dumps([{key: _json_default(value) for key, value in row.items()} for row in rows])


def parse_variant_rows_json(raw: str | None) -> list[dict[str, Any]]:
    if raw is None or (isinstance(raw, float) and pd.isna(raw)) or raw == "":
        return []
    parsed = json.loads(raw)
    if not isinstance(parsed, list):
        raise ValueError("background_variants_json must encode a list of variant rows")
    return [dict(item) for item in parsed]


def build_track_rows_from_track_data(
    track_data: track_data_lib.TrackData,
    *,
    selection_metadata: Mapping[str, Any],
    sequence_role: str,
    sequence_label: str,
    center_variant_pos: int | None = None,
    output_type: str = "RNA_SEQ",
) -> pd.DataFrame:
    values = np.asarray(track_data.values)
    if values.ndim == 1:
        values = values[None, :]

    metadata = normalize_track_columns(track_data.metadata, default_output_type=output_type)
    rows: list[dict[str, Any]] = []
    for track_idx, (_, meta_row) in enumerate(metadata.iterrows()):
        per_track = np.asarray(values[..., track_idx])
        flat = per_track.astype(float).reshape(-1)
        row = dict(selection_metadata)
        row.update(
            {
                "sequence_role": sequence_role,
                "sequence_label": sequence_label,
                "output_type": output_type,
                "track_name": meta_row.get("track_name", pd.NA),
                "track_strand": meta_row.get("track_strand", pd.NA),
                "assay_title": meta_row.get("assay_title", pd.NA),
                "ontology_curie": meta_row.get("ontology_curie", pd.NA),
                "biosample_name": meta_row.get("biosample_name", pd.NA),
                "biosample_type": meta_row.get("biosample_type", pd.NA),
                "gtex_tissue": meta_row.get("gtex_tissue", pd.NA),
                "track_resolution_bp": int(track_data.resolution),
                "track_shape_json": json.dumps(list(per_track.shape)),
                "track_values_json": json.dumps(per_track.tolist()),
                "track_mean": float(np.nanmean(flat)),
                "track_sum": float(np.nansum(flat)),
                "track_min": float(np.nanmin(flat)),
                "track_max": float(np.nanmax(flat)),
                "track_num_values": int(flat.size),
            }
        )

        center_bin_idx = pd.NA
        center_bin_value = pd.NA
        if (
            center_variant_pos is not None
            and track_data.interval is not None
            and per_track.ndim >= 1
            and track_data.resolution > 0
        ):
            relative_bp = (int(center_variant_pos) - 1) - int(track_data.interval.start)
            if relative_bp >= 0:
                center_bin_idx = int(relative_bp // int(track_data.resolution))
                if center_bin_idx < per_track.shape[0]:
                    center_bin_value = float(np.asarray(per_track)[center_bin_idx].mean())

        row["center_bin_index"] = center_bin_idx
        row["center_bin_value"] = center_bin_value
        rows.append(row)

    out = pd.DataFrame(rows)
    return add_track_identity(out, default_output_type=output_type)


def reduce_predict_sequence_delta(
    background_rows: pd.DataFrame,
    alternate_rows: pd.DataFrame,
) -> pd.DataFrame:
    bg = add_track_identity(background_rows, default_output_type="RNA_SEQ").copy()
    alt = add_track_identity(alternate_rows, default_output_type="RNA_SEQ").copy()

    join_cols = ["selection_id", "track_id"]
    keep_cols = [
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
    ]
    merged = bg[join_cols + keep_cols + ["track_mean", "track_sum", "center_bin_value"]].merge(
        alt[join_cols + ["track_mean", "track_sum", "center_bin_value"]],
        on=join_cols,
        how="inner",
        suffixes=("_background_only", "_background_plus_singleton"),
    )
    merged["mode"] = "haplotype_predict_sequence"
    merged["reducer_name"] = "mean_alt_minus_background_full_track"
    merged["scalar_score"] = (
        merged["track_mean_background_plus_singleton"] - merged["track_mean_background_only"]
    )
    merged["scalar_score_center_bin_delta"] = (
        merged["center_bin_value_background_plus_singleton"] - merged["center_bin_value_background_only"]
    )
    return merged
