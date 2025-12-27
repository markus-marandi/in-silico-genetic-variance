from __future__ import annotations

import os
import shutil
from dataclasses import dataclass
from pathlib import Path
import sys


DEFAULT_ROOT = Path(os.getenv("PDC_TMP", "/cfs/klemming/scratch/m/mmarandi"))


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
    def run_info(self) -> Path:
        return self.results_dir / "run_info.json"

    def chunk_path(self, start: int, end: int) -> Path:
        return self.chunks_dir / f"chunk_{start:07d}_{end:07d}.tsv.gz"

    def variant_parquet(self, tag: str = "annotated") -> Path:
        return self.results_dir / f"{self.sample_id}.variants.{tag}.parquet"

    def gene_parquet(self, tag: str = "summary") -> Path:
        return self.results_dir / f"{self.sample_id}.genes.{tag}.parquet"

    @property
    def is_ism_run(self) -> bool:
        return "ism" in self.sample_id.lower()

    def make_dirs(self) -> None:
        """ensure experiment directories exist."""
        for path in (self.inputs_dir, self.chunks_dir, self.results_dir):
            path.mkdir(parents=True, exist_ok=True)

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
