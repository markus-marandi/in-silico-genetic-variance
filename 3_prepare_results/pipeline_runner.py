from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

# --- FIX 1: Add current directory to sys.path to allow absolute imports ---
# This fixes "ImportError: attempted relative import with no known parent package"
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from modules.stitcher import stitch_variants, load_gene_list
from modules.normalizer import normalize_and_backfill
from modules.annotator import annotate_af, annotate_gnomad
from modules.aggregator import aggregate_genes
from modules.external_data_loader import ExternalDataLoader

DATE_FMT = "%Y%m%d"

def _load_project_layout():
    # Adjusted path logic to work reliably
    helper_path = HERE.parent / "helpers" / "path_manager.py"
    
    # Check if file exists to give a better error
    if not helper_path.exists():
        raise FileNotFoundError(f"Path manager helper not found at: {helper_path}")

    spec = importlib.util.spec_from_file_location("path_manager", helper_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"unable to load path manager from {helper_path}")
    
    module = importlib.util.module_from_spec(spec)
    # FIX 3: Register in sys.modules to prevent dataclass errors later
    sys.modules["path_manager"] = module
    
    spec.loader.exec_module(module)  # type: ignore[arg-type]
    return module.ProjectLayout

@dataclass(frozen=True)
class PipelineSpec:
    dataset_id: str
    sample_id: str
    root_dir: Path
    chunks_dir: Path
    variant_output: Path
    gene_output: Path
    variants_af_path: Path
    gnomad_path: Path | None
    is_ism: bool

    @classmethod
    def from_args(
        cls,
        dataset_id: str | None,
        sample_id: str | None,
        root_dir: str | None,
        chunk_dir: str | None,
        variant_out: str | None,
        gene_out: str | None,
        variants_af: str | None,
        gnomad_af: str | None,
    ) -> "PipelineSpec":
        layout_cls = _load_project_layout()
        base_root = Path(root_dir) if root_dir else Path(
            os.getenv("ROOT_DIR") or os.getenv("PDC_TMP") or "/cfs/klemming/scratch/m/mmarandi"
        )

        if dataset_id and sample_id:
            layout = layout_cls(dataset_id=dataset_id, sample_id=sample_id, root_dir=base_root)
            layout.make_dirs()
            resolved_chunk_dir = Path(chunk_dir) if chunk_dir else layout.chunks_dir
        else:
            if not chunk_dir:
                raise ValueError("chunks_dir is required when dataset_id/sample_id are not provided")
            resolved_chunk_dir = Path(chunk_dir).resolve()
            sample_id, dataset_id = cls._infer_ids_from_chunks(resolved_chunk_dir)
            layout = layout_cls(dataset_id=dataset_id, sample_id=sample_id, root_dir=base_root)
            layout.make_dirs()

        tag = datetime.now().strftime(DATE_FMT)
        variant_path = Path(variant_out).resolve() if variant_out else layout.results_dir / f"{layout.sample_id}_variants_{tag}.parquet"
        gene_path = Path(gene_out).resolve() if gene_out else layout.results_dir / f"{layout.sample_id}_genes_{tag}.parquet"
        variants_af_path = Path(variants_af).resolve() if variants_af else layout.inputs_dir / "variants.tsv.gz"

        return cls(
            dataset_id=layout.dataset_id,
            sample_id=layout.sample_id,
            root_dir=layout.root_dir,
            chunks_dir=resolved_chunk_dir,
            variant_output=variant_path,
            gene_output=gene_path,
            variants_af_path=variants_af_path,
            gnomad_path=Path(gnomad_af) if gnomad_af else None,
            is_ism="ism" in layout.sample_id.lower(),
        )

    @staticmethod
    def _infer_ids_from_chunks(chunks_dir: Path) -> tuple[str, str]:
        parts = chunks_dir.resolve().parts
        if len(parts) < 3:
            raise ValueError(f"cannot infer dataset/sample from {chunks_dir}")
        sample_id = parts[-2]
        dataset_id = parts[-3]
        return sample_id, dataset_id


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Stitch, annotate, and summarize chunked scores.")
    parser.add_argument("--chunks-dir", type=Path, help="directory containing chunk_*.tsv.gz")
    parser.add_argument("--dataset-id", type=str, help="dataset id (preferred)")
    parser.add_argument("--sample-id", type=str, help="sample id (preferred)")
    parser.add_argument("--root-dir", type=Path, help="override root dir (defaults to $PDC_TMP)")
    parser.add_argument("--variant-out", type=Path, help="override variant parquet output path")
    parser.add_argument("--gene-out", type=Path, help="override gene parquet output path")
    parser.add_argument("--gene-list", type=Path, help="optional gene whitelist (one id per line)")
    parser.add_argument("--gnomad-af", type=Path, help="optional gnomAD parquet with CHROM,POS,REF,ALT,AF")
    parser.add_argument("--variants-af", type=Path, help="initial variants tsv/tsv.gz with AF (defaults to 01_inputs/variants.tsv.gz)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    spec = PipelineSpec.from_args(
        dataset_id=args.dataset_id or os.getenv("DATASET_ID"),
        sample_id=args.sample_id or os.getenv("SAMPLE_ID"),
        root_dir=str(args.root_dir) if args.root_dir else os.getenv("ROOT_DIR") or os.getenv("PDC_TMP"),
        chunk_dir=str(args.chunks_dir) if args.chunks_dir else None,
        variant_out=str(args.variant_out) if args.variant_out else None,
        gene_out=str(args.gene_out) if args.gene_out else None,
        variants_af=str(args.variants_af) if args.variants_af else None,
        gnomad_af=str(args.gnomad_af) if args.gnomad_af else None,
    )

    gene_list = args.gene_list
    if gene_list is None:
        default_gene = spec.root_dir / "experiments" / spec.dataset_id / spec.sample_id / "01_inputs" / "gene_list.tsv"
        if default_gene.exists():
            gene_list = default_gene

    print(f"Stitching variants for {spec.sample_id}...")
    lf = stitch_variants(spec.chunks_dir, gene_list_path=gene_list)
    lf = normalize_and_backfill(lf)

    print("Annotating AF...")
    lf, has_af = annotate_af(lf, spec.variants_af_path)
    
    if spec.gnomad_path:
        print("Annotating gnomAD...")
        lf = annotate_gnomad(lf, spec.gnomad_path)

    print(f"Writing variants to {spec.variant_output}...")
    # Ensure parent dir exists
    spec.variant_output.parent.mkdir(parents=True, exist_ok=True)
    lf.sink_parquet(spec.variant_output, compression="zstd")

    is_ism = spec.is_ism or not has_af
    print(f"Aggregating genes (ISM mode={is_ism})...")
    aggregate_genes(spec.variant_output, spec.gene_output, base_ref=spec.root_dir, is_ism=is_ism)
    print("Done.")

if __name__ == "__main__":
    main()