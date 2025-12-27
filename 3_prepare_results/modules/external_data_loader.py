from __future__ import annotations

from pathlib import Path

import pandas as pd
import polars as pl

from modules.normalisation_helper import (
    load_mane_select,
    parse_gtf_gene_with_lengths,
    parse_gtf_utr5_lengths,
    load_tpm_muscle,
    load_vgh_metrics,
)


class ExternalDataLoader:
    """bridge: load metadata via pandas helpers and convert to polars."""

    def __init__(self, base: Path):
        self.base = base
        self.mane_path = base / "initial_data_external/MANE.GRCh38.v1.4.summary.txt"
        self.gtf_path = base / "initial_data_external/gencode.v49.annotation.gtf"
        self.tpm_path = base / "initial_data_external/GTEx_Analysis_2022-06-06_v10_RNASeQCv2.4.2_gene_median_tpm.gct"
        self.vgh_path = base / "initial_data_external/gene_metrics_vgh_202407.tsv.gz"
        self.aneva_path = base / "initial_data_external/Vg_GTEx_v7.tsv.gz"

    @staticmethod
    def _to_polars(df: pd.DataFrame) -> pl.DataFrame:
        pldf = pl.from_pandas(df)
        if "gene_id" in pldf.columns:
            pldf = pldf.with_columns(pl.col("gene_id").cast(pl.Utf8))
        return pldf

    def load_mane(self) -> pl.DataFrame:
        if not self.mane_path.exists():
            raise FileNotFoundError(f"MANE summary not found: {self.mane_path}")
        df = load_mane_select(self.mane_path)
        return self._to_polars(df)

    def load_gtf_genes(self) -> pl.DataFrame:
        if not self.gtf_path.exists():
            raise FileNotFoundError(f"GTF not found: {self.gtf_path}")
        df_genes = parse_gtf_gene_with_lengths(self.gtf_path)
        try:
            df_mane = load_mane_select(self.mane_path)
            df_utr5 = parse_gtf_utr5_lengths(self.gtf_path, df_mane)
            df_genes = df_genes.merge(df_utr5, on="gene_id", how="left")
        except Exception as exc:  # pragma: no cover - warning only
            print(f"warning: could not calculate UTR5 lengths: {exc}")
        return self._to_polars(df_genes)

    def load_tpm(self) -> pl.DataFrame:
        if not self.tpm_path.exists():
            raise FileNotFoundError(f"TPM GCT not found: {self.tpm_path}")
        df = load_tpm_muscle(self.tpm_path)
        return self._to_polars(df)

    def load_vgh(self) -> pl.DataFrame:
        if not self.vgh_path.exists():
            raise FileNotFoundError(f"VGH metrics not found: {self.vgh_path}")
        df = load_vgh_metrics(self.vgh_path)
        return self._to_polars(df)

    def load_aneva(self) -> pl.DataFrame:
        """load Mohammadi et al. genetic variance reference; expects TSV/CSV/Parquet."""
        if self.aneva_path.suffix in {".tsv", ".gz", ".bgz"} and self.aneva_path.exists():
            df = pl.read_csv(self.aneva_path, separator="\t", null_values=["", "NA", "None"])
            return df
        if self.aneva_path.suffix == ".parquet" and self.aneva_path.exists():
            return pl.read_parquet(self.aneva_path)
        # fallback: attempt pyreadr for .rda/.rds
        if self.aneva_path.suffix in {".rda", ".rds"}:
            try:
                import pyreadr  # type: ignore
            except Exception as exc:  # pragma: no cover
                raise ImportError(
                    f"pyreadr not available to load {self.aneva_path}; convert to tsv/parquet."
                ) from exc
            res = pyreadr.read_r(self.aneva_path)  # type: ignore
            if not res:
                raise ValueError(f"no tables found in {self.aneva_path}")
            df_pd = next(iter(res.values()))
            df = pl.from_pandas(df_pd)
            return df
        raise FileNotFoundError(f"aneva reference not found or unsupported format: {self.aneva_path}")
