from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = Path(
    os.getenv(
        "PDC_TMP",
        REPO_ROOT / "2_score_variants_with_alphagenome" / "tmp",
    )
)


@dataclass(frozen=True)
class ProjectLayout:
    """derive dataset-centric paths for a single experiment run.

    args:
        dataset_id (str): logical dataset grouping, e.g. dataset5_null.
        sample_id (str): run-specific sample name, e.g. background_ISM.
        root_dir (Path): base scratch root; defaults to $PDC_TMP or shared scratch.

    returns:
        ProjectLayout: immutable view of the experiment tree.
    """

    dataset_id: str
    sample_id: str
    root_dir: Path = DEFAULT_ROOT

    @property
    def experiment_root(self) -> Path:
        return self.root_dir / "experiments" / self.dataset_id / self.sample_id

    @property
    def inputs_dir(self) -> Path:
        return self.experiment_root / "01_inputs"

    @property
    def chunks_dir(self) -> Path:
        return self.experiment_root / "02_chunks"

    @property
    def results_dir(self) -> Path:
        return self.experiment_root / "03_results"

    @property
    def input_variants(self) -> Path:
        return self.inputs_dir / "variants.tsv.gz"

    @property
    def input_targets(self) -> Path:
        return self.inputs_dir / "targets.bed"

    @property
    def input_gene_list(self) -> Path:
        return self.inputs_dir / "gene_list.tsv"

    @property
    def run_info(self) -> Path:
        return self.results_dir / "run_info.json"

    def chunk_path(self, start: int, end: int) -> Path:
        return self.chunks_dir / f"chunk_{start:07d}_{end:07d}.tsv.gz"

    def variant_parquet(self, tag: str = "annotated") -> Path:
        return self.results_dir / f"{self.sample_id}.variants.{tag}.parquet"

    def gene_parquet(self, tag: str = "summary") -> Path:
        return self.results_dir / f"{self.sample_id}.genes.{tag}.parquet"

    @property
    def singleton_haplotype_inputs_dir(self) -> Path:
        return self.inputs_dir / "singleton_haplotype"

    @property
    def singleton_haplotype_jobs_parquet(self) -> Path:
        return self.singleton_haplotype_inputs_dir / "singleton_haplotype_jobs.parquet"

    @property
    def singleton_haplotype_jobs_tsv(self) -> Path:
        return self.singleton_haplotype_inputs_dir / "singleton_haplotype_jobs.tsv.gz"

    @property
    def singleton_haplotype_background_parquet(self) -> Path:
        return self.singleton_haplotype_inputs_dir / "singleton_haplotype_background_variants.parquet"

    @property
    def singleton_haplotype_background_tsv(self) -> Path:
        return self.singleton_haplotype_inputs_dir / "singleton_haplotype_background_variants.tsv.gz"

    @property
    def is_ism_run(self) -> bool:
        return ("ism" in self.sample_id.lower()) or ("null" in self.sample_id.lower())

    def workflow_results_dir(self, workflow: str) -> Path:
        return self.results_dir / workflow

    def make_dirs(self) -> None:
        """ensure experiment directories exist."""
        for path in (
            self.inputs_dir,
            self.chunks_dir,
            self.results_dir,
            self.singleton_haplotype_inputs_dir,
        ):
            path.mkdir(parents=True, exist_ok=True)

    def _resolve_existing_path(
        self,
        preferred: Path,
        patterns: tuple[str, ...],
        *,
        description: str,
        required: bool = True,
    ) -> Path | None:
        if preferred.exists():
            return preferred

        matches: list[Path] = []
        for pattern in patterns:
            matches.extend(sorted(self.inputs_dir.glob(pattern)))

        unique_matches: list[Path] = []
        seen: set[Path] = set()
        for match in matches:
            if match in seen or not match.exists():
                continue
            unique_matches.append(match)
            seen.add(match)

        if len(unique_matches) == 1:
            return unique_matches[0]
        if len(unique_matches) > 1:
            msg = ", ".join(path.name for path in unique_matches)
            raise FileExistsError(
                f"multiple candidate {description} files found under {self.inputs_dir}: {msg}"
            )
        if required:
            raise FileNotFoundError(
                f"missing {description} under {self.inputs_dir}; checked {preferred.name}"
            )
        return None

    def resolve_input_variants(self) -> Path:
        return self._resolve_existing_path(
            self.input_variants,
            ("variants.tsv.gz", "*_variants.tsv.gz", "*.variants.tsv.gz", "*_variants.tsv"),
            description="variant table",
        )

    def resolve_input_targets(self) -> Path:
        return self._resolve_existing_path(
            self.input_targets,
            ("targets.bed", "*_tss*.bed", "*.bed"),
            description="target BED",
        )

    def resolve_input_gene_list(self, required: bool = False) -> Path | None:
        return self._resolve_existing_path(
            self.input_gene_list,
            ("gene_list.tsv", "ensg_list.txt", "*gene*list*.tsv", "*gene*subset*.tsv"),
            description="gene list",
            required=required,
        )

    @classmethod
    def from_env(cls) -> "ProjectLayout":
        """build layout from env vars DATASET_ID, SAMPLE_ID, ROOT_DIR/PDC_TMP."""
        dataset_id = os.getenv("DATASET_ID")
        sample_id = os.getenv("SAMPLE_ID")
        if not dataset_id or not sample_id:
            raise ValueError("DATASET_ID and SAMPLE_ID must be set")
        root = Path(os.getenv("ROOT_DIR", os.getenv("PDC_TMP", str(DEFAULT_ROOT))))
        return cls(dataset_id=dataset_id, sample_id=sample_id, root_dir=root)


def migrate_chunks(old_chunk_dir: Path, layout: ProjectLayout, move: bool = False) -> list[Path]:
    """relocate legacy chunk_*.tsv.gz files into the dataset-centric tree.

    args:
        old_chunk_dir (Path): source directory containing chunk_*.tsv.gz.
        layout (ProjectLayout): destination layout.
        move (bool): move instead of copy when true.

    returns:
        list[Path]: paths written under the new layout.
    """
    if not old_chunk_dir.exists():
        raise FileNotFoundError(f"source chunk dir not found: {old_chunk_dir}")

    layout.make_dirs()
    written: list[Path] = []
    for src in sorted(old_chunk_dir.glob("chunk_*.tsv.gz")):
        dest = layout.chunks_dir / src.name
        if dest.exists():
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        if move:
            shutil.move(src.as_posix(), dest.as_posix())
        else:
            shutil.copy2(src.as_posix(), dest.as_posix())
        written.append(dest)
    return written


def copy_inputs(src_variants: Path, src_targets: Path, layout: ProjectLayout) -> None:
    """copy prepared inputs into the layout inputs_dir."""
    layout.make_dirs()
    for source, dest in (
        (src_variants, layout.input_variants),
        (src_targets, layout.input_targets),
    ):
        if not source.exists():
            raise FileNotFoundError(f"missing input: {source}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source.as_posix(), dest.as_posix())


@dataclass(frozen=True)
class SyntheticVariantPaths:
    """manage paths for synthetic variant generation workflow.
    
    args:
        dataset_name (str): name of the dataset, e.g. ClinGen or Background.
        base_dir (Path): base directory for intermediate files.
        output_dir (Path): directory for final outputs.
    
    returns:
        SyntheticVariantPaths: immutable view of paths.
    """
    
    dataset_name: str
    base_dir: Path
    output_dir: Path
    
    @property
    def raw_variants_tsv(self) -> Path:
        return self.base_dir / f"{self.dataset_name}_variants.tsv"
    
    @property
    def keyed_ht(self) -> Path:
        return self.base_dir / f"{self.dataset_name}_variants.keyed.ht"
    
    @property
    def bed_file(self) -> Path:
        return self.base_dir / f"{self.dataset_name}_gene_set±10kb.tab.bed"
    
    @property
    def candidates_ht(self) -> Path:
        return self.output_dir / f"{self.dataset_name}.matched_synthetic.candidates.ht"
    
    @property
    def sampled_ht(self) -> Path:
        return self.output_dir / f"{self.dataset_name}.matched_synthetic.sampled.ht"
    
    @property
    def sampled_tsv(self) -> Path:
        return self.output_dir / f"{self.dataset_name}.matched_synthetic.sampled.tsv.bgz"
    
    @property
    def observed_dist_plot(self) -> Path:
        return self.output_dir / f"{self.dataset_name}.observed_dist.png"
    
    @property
    def sampled_dist_plot(self) -> Path:
        return self.output_dir / f"{self.dataset_name}.sampled_dist.png"
