from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import polars as pl

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from modules.stitcher import stitch_variants, load_gene_list
from modules.normalizer import normalize_and_backfill
from modules.annotator import annotate_af, annotate_gnomad
from modules.aggregator import aggregate_genes
from modules.external_data_loader import ExternalDataLoader
from modules.variant_deduplicator import deduplicate_by_gene_and_variant
from modules.synthetic_variant_downsampler import load_real_variant_counts, downsample_to_real_counts
from modules.permuted_af import load_gene_af_pools, add_perm_af_gene_aware

DATE_FMT = "%Y%m%d"

def _load_project_layout():
    helper_path = HERE.parent / "helpers" / "path_manager.py"
    
    if not helper_path.exists():
        raise FileNotFoundError(f"Path manager helper not found at: {helper_path}")

    spec = importlib.util.spec_from_file_location("path_manager", helper_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"unable to load path manager from {helper_path}")
    
    module = importlib.util.module_from_spec(spec)
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
    real_reference: Path | None

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
        real_ref: str | None,
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
        is_ism = ("ism" in layout.sample_id.lower()) or ("null" in layout.sample_id.lower())
        
        # all datasets get _dedup, ism/null also get _downsampled
        if not variant_out:
            if is_ism:
                variant_path = layout.results_dir / f"{layout.sample_id}_variants_dedup_downsampled_{tag}.parquet"
            else:
                variant_path = layout.results_dir / f"{layout.sample_id}_variants_dedup_{tag}.parquet"
        else:
            variant_path = Path(variant_out).resolve()
        
        gene_path = Path(gene_out).resolve() if gene_out else layout.results_dir / f"{layout.sample_id}_genes_{tag}.parquet"
        variants_af_path = Path(variants_af).resolve() if variants_af else layout.inputs_dir / "variants.tsv.gz"
        real_reference = Path(real_ref).resolve() if real_ref else None

        return cls(
            dataset_id=layout.dataset_id,
            sample_id=layout.sample_id,
            root_dir=layout.root_dir,
            chunks_dir=resolved_chunk_dir,
            variant_output=variant_path,
            gene_output=gene_path,
            variants_af_path=variants_af_path,
            gnomad_path=Path(gnomad_af) if gnomad_af else None,
            is_ism=is_ism,
            real_reference=real_reference,
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
    parser.add_argument("--variants-parquet", type=Path, help="existing variants parquet to aggregate directly")
    parser.add_argument("--real-reference", type=Path, help="real dataset reference for downsampling (gene-level or variant-level parquet, required for ISM/NULL)")
    parser.add_argument("--deduplicate", action="store_true", help="deduplicate existing variants parquet (and downsample if ISM/NULL)")
    parser.add_argument("--permute-af", action="store_true", help="Generate perm_AF from real reference and compute vg_predicted_perm")
    parser.add_argument("--calc-ci", action="store_true", help="Calculate 5th/95th CI via Monte Carlo")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.variants_parquet:
        base_root = Path(args.root_dir) if args.root_dir else Path(
            os.getenv("ROOT_DIR") or os.getenv("PDC_TMP") or "/cfs/klemming/scratch/m/mmarandi"
        )
        variant_path = Path(args.variants_parquet).resolve()
        
        # 1. Infer dataset type
        is_ism = "ism" in str(variant_path).lower() or "null" in str(variant_path).lower()
        
        # 2. Handle Processing
        if args.deduplicate:
            print(f"Loading variants from {variant_path}...")
            df = pl.read_parquet(variant_path)

            # [OPTIMIZATION] Filter FIRST, before expensive deduplication
            if args.gene_list:
                print(f"Filtering variants using whitelist: {args.gene_list}")
                whitelist_df = pl.read_csv(
                    args.gene_list, 
                    has_header=False, 
                    new_columns=["gene_id"], 
                    separator="\t"
                )
                whitelist = whitelist_df["gene_id"].to_list()
                df = df.filter(pl.col("gene_id").is_in(whitelist))
                print(f"  Rows after filtering: {len(df)}")
            
            # Step A: Rigid Deduplication
            print("Deduplicating by gene_id and variant_id...")
            df = deduplicate_by_gene_and_variant(df, verbose=True)
            
            # Step B: Downsampling (ISM/NULL only)
            if is_ism:
                if not args.real_reference:
                    raise ValueError("Flag --real-reference is required for ISM/NULL processing.")
                
                real_ref = Path(args.real_reference).resolve()
                if not real_ref.exists():
                    raise FileNotFoundError(f"Real reference not found: {real_ref}")
                
                print(f"Loading real variant counts from {real_ref}...")
                real_counts = load_real_variant_counts(real_ref)
                
                print("Downsampling to match real counts...")
                df = downsample_to_real_counts(df, real_counts, seed=42, verbose=True)
                
                suffix = "_downsampled"
            else:
                suffix = "_dedup"

            # Step C: AF Permutation (For both ISM and Sanity Check)
            if args.permute_af:
                print("--- Starting AF Permutation ---")
                perm_ref_path = None
                
                if args.real_reference:
                    # Use explicit reference if provided (Standard for ISM)
                    perm_ref_path = Path(args.real_reference).resolve()
                elif not is_ism: 
                    # Sanity Check Logic: Use input variants as source if no ref provided
                    print("  No external reference provided for real data; using input variants as AF source.")
                    perm_ref_path = variant_path
                
                if not perm_ref_path:
                    raise ValueError("Cannot permute AF: No reference provided (use --real-reference).")

                if not perm_ref_path.exists():
                    raise FileNotFoundError(f"Permutation reference not found: {perm_ref_path}")

                print(f"  Loading AF pools from {perm_ref_path}...")
                af_pools = load_gene_af_pools(perm_ref_path)

                print("  Sampling and assigning perm_AF...")
                df = add_perm_af_gene_aware(df, af_pools, seed=42)
                
                suffix += "_perm"

            # Determine output path safely
            if args.variant_out:
                processed_path = Path(args.variant_out).resolve()
            else:
                # Auto-naming: check if suffix exists to avoid duplication
                stem = variant_path.stem
                if not stem.endswith(suffix):
                    stem += suffix
                processed_path = variant_path.with_name(f"{stem}.parquet")

            print(f"Writing processed variants to {processed_path}...")
            df.write_parquet(processed_path)
            
            # Update variant_path so aggregation uses the NEW file
            variant_path = processed_path
        
        # 3. Aggregate Genes
        gene_path = Path(args.gene_out).resolve() if args.gene_out else variant_path.with_name(f"{variant_path.stem}_genes.parquet")
        gene_list = Path(args.gene_list).resolve() if args.gene_list else None

        # Determine reference for CI
        ci_ref = None
        if args.real_reference:
            ci_ref = Path(args.real_reference).resolve()
        elif not is_ism:
            ci_ref = variant_path # Self-reference for sanity check

        print(f"Aggregating genes from {variant_path}...")
        aggregate_genes(
            variant_path,
            gene_path,
            base_ref=base_root,
            is_ism=is_ism,
            gene_list_path=gene_list,
            calculate_ci=args.calc_ci,
            real_reference_path=ci_ref, # [FIXED] Pass the calculated ci_ref, not args.real_reference
            n_permutations=1000,
        )
        print("Done.")
        return

    # ... (rest of the file is fine) ...
    spec = PipelineSpec.from_args(
        dataset_id=args.dataset_id or os.getenv("DATASET_ID"),
        sample_id=args.sample_id or os.getenv("SAMPLE_ID"),
        root_dir=str(args.root_dir) if args.root_dir else os.getenv("ROOT_DIR") or os.getenv("PDC_TMP"),
        chunk_dir=str(args.chunks_dir) if args.chunks_dir else None,
        variant_out=str(args.variant_out) if args.variant_out else None,
        gene_out=str(args.gene_out) if args.gene_out else None,
        variants_af=str(args.variants_af) if args.variants_af else None,
        gnomad_af=str(args.gnomad_af) if args.gnomad_af else None,
        real_ref=str(args.real_reference) if args.real_reference else None,
    )
    
    if spec.is_ism and not spec.real_reference:
        raise ValueError(
            f"ISM/NULL dataset '{spec.sample_id}' requires --real-reference for downsampling"
        )
    
    if spec.real_reference and not spec.real_reference.exists():
        raise FileNotFoundError(f"real reference not found: {spec.real_reference}")

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

    is_ism = spec.is_ism or not has_af
    
    print("Collecting and deduplicating...")
    df = lf.collect()
    
    print("  deduplicating by gene_id and variant_id...")
    df = deduplicate_by_gene_and_variant(df, verbose=True)
    
    if is_ism:
        print(f"  loading real variant counts from {spec.real_reference}...")
        real_counts = load_real_variant_counts(spec.real_reference)
        
        print("  downsampling to match real counts...")
        df = downsample_to_real_counts(df, real_counts, seed=42, verbose=True)
    
    print(f"Writing variants to {spec.variant_output}...")
    spec.variant_output.parent.mkdir(parents=True, exist_ok=True)
    df.write_parquet(spec.variant_output, compression="zstd")

    print(f"Aggregating genes (ISM mode={is_ism})...")
    aggregate_genes(
        spec.variant_output,
        spec.gene_output,
        base_ref=spec.root_dir,
        is_ism=is_ism,
        gene_list_path=gene_list,
    )
    print("Done.")

if __name__ == "__main__":
    main()