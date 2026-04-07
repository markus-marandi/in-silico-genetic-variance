from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_PLOTTING = True
except ModuleNotFoundError:
    matplotlib = None
    plt = None
    sns = None
    HAS_PLOTTING = False


REPO_ROOT = Path(__file__).resolve().parent.parent
PATH_MANAGER = REPO_ROOT / "helpers/path_manager.py"


def _load_layout():
    if not PATH_MANAGER.exists():
        raise FileNotFoundError(f"path manager not found: {PATH_MANAGER}")

    spec = importlib.util.spec_from_file_location("ag_path_manager", PATH_MANAGER)
    if spec is None or spec.loader is None:
        raise ImportError(f"unable to load path manager from {PATH_MANAGER}")

    module = importlib.util.module_from_spec(spec)
    sys.modules["ag_path_manager"] = module
    spec.loader.exec_module(module)
    return module.ProjectLayout.from_env()


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot AlphaGenome batch stability outputs.")
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help="Batch stability result directory. Defaults to the dataset/sample layout path.",
    )
    parser.add_argument(
        "--top-k",
        type=int,
        default=6,
        help="Number of highest-variability variants to show in the trajectory figure.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory where summaries and plots should be written. Defaults to results-dir.",
    )
    return parser.parse_args()


def _resolve_results_dir(cli_results_dir: Path | None) -> Path:
    if cli_results_dir is not None:
        return cli_results_dir
    layout = _load_layout()
    return layout.results_dir / "batch_stability"


def _parse_variant_id(variant_id: str) -> tuple[str, int, str]:
    chrom, pos, change = str(variant_id).split(":", 2)
    return chrom, int(pos), change


def _short_variant_label(variant_id: str, max_change_len: int = 18) -> str:
    chrom, pos, change = _parse_variant_id(variant_id)
    if len(change) > max_change_len:
        change = f"{change[:max_change_len - 1]}..."
    return f"{chrom}:{pos} {change}"


def _load_iteration_frame(iterations_dir: Path) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    pattern = re.compile(r"iteration_(\d+)\.tsv\.gz$")

    for path in sorted(iterations_dir.glob("iteration_*.tsv.gz")):
        match = pattern.match(path.name)
        if match is None:
            continue
        iteration = int(match.group(1))
        frame = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
        frame["iteration"] = iteration
        frames.append(frame)

    if not frames:
        raise FileNotFoundError(f"no iteration files found in {iterations_dir}")

    all_iterations = pd.concat(frames, ignore_index=True)
    required = {"iteration", "variant_id", "gene_id", "raw_score"}
    missing = required - set(all_iterations.columns)
    if missing:
        raise ValueError(f"iteration outputs are missing required columns: {sorted(missing)}")

    all_iterations["variant_id"] = all_iterations["variant_id"].astype("string").str.strip()
    all_iterations["gene_id"] = all_iterations["gene_id"].astype("string").str.strip()
    all_iterations["raw_score"] = pd.to_numeric(all_iterations["raw_score"], errors="coerce")
    all_iterations = all_iterations.dropna(subset=["variant_id", "gene_id", "raw_score"]).copy()
    if all_iterations.empty:
        raise ValueError("iteration outputs contain no usable rows after cleaning")

    assay_col = "Assay title" if "Assay title" in all_iterations.columns else (
        "Assay_title" if "Assay_title" in all_iterations.columns else None
    )
    assay_source = all_iterations[assay_col] if assay_col else pd.Series(pd.NA, index=all_iterations.index)
    track_source = all_iterations["track_name"] if "track_name" in all_iterations.columns else pd.Series(pd.NA, index=all_iterations.index)
    protocol_source = assay_source.fillna(track_source)
    protocol_norm = (
        protocol_source.astype("string")
        .fillna("")
        .str.lower()
        .str.replace(r"[^a-z0-9]+", "_", regex=True)
        .str.strip("_")
    )

    all_iterations["protocol_group"] = np.where(
        protocol_norm.str.contains("polya") & protocol_norm.str.contains("rna_seq"),
        "polyA_plus_rna_seq",
        np.where(
            protocol_norm.str.contains("total") & protocol_norm.str.contains("rna_seq"),
            "total_rna_seq",
            "other",
        ),
    )

    return all_iterations


def _summarize_protocol_competition(
    agg: pd.DataFrame,
    iterations: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    pair_summary_rows: list[dict[str, object]] = []
    pair_protocol_rows: list[dict[str, object]] = []

    iter_protocol_counts = (
        iterations.groupby(["iteration", "protocol_group"])
        .size()
        .rename("winner_count")
        .reset_index()
        .pivot(index="iteration", columns="protocol_group", values="winner_count")
        .fillna(0)
        .astype(int)
        .reset_index()
    )

    for (variant_id, gene_id), subset in iterations.groupby(["variant_id", "gene_id"], sort=False):
        subset = subset.sort_values("iteration").copy()
        protocol_stats = (
            subset.groupby("protocol_group")["raw_score"]
            .agg(["count", "mean", "std"])
            .reset_index()
            .rename(columns={"count": "n_iterations", "mean": "mean_score", "std": "std_score"})
        )
        protocol_stats["std_score"] = protocol_stats["std_score"].fillna(0.0)

        overall_mean = float(subset["raw_score"].mean())
        ss_total = float(((subset["raw_score"] - overall_mean) ** 2).sum())
        ss_between = float(
            (protocol_stats["n_iterations"] * (protocol_stats["mean_score"] - overall_mean) ** 2).sum()
        )
        if subset["protocol_group"].nunique() <= 1 or ss_total <= 1e-15:
            protocol_eta_sq = 0.0
        else:
            protocol_eta_sq = ss_between / ss_total

        protocol_switches = int((subset["protocol_group"] != subset["protocol_group"].shift()).sum() - 1)
        protocol_switch_rate = protocol_switches / max(len(subset) - 1, 1)
        within_protocol_weighted_std = float(
            np.average(protocol_stats["std_score"], weights=protocol_stats["n_iterations"])
        )

        dominant_protocol = (
            protocol_stats.sort_values(["n_iterations", "mean_score"], ascending=[False, False])
            .iloc[0]["protocol_group"]
        )
        protocol_counts = {
            row["protocol_group"]: int(row["n_iterations"])
            for _, row in protocol_stats.iterrows()
        }

        poly_mean = protocol_stats.loc[
            protocol_stats["protocol_group"] == "polyA_plus_rna_seq", "mean_score"
        ]
        total_mean = protocol_stats.loc[
            protocol_stats["protocol_group"] == "total_rna_seq", "mean_score"
        ]
        protocol_gap = (
            float(abs(poly_mean.iloc[0] - total_mean.iloc[0]))
            if not poly_mean.empty and not total_mean.empty
            else np.nan
        )

        pair_summary_rows.append(
            {
                "variant_id": variant_id,
                "gene_id": gene_id,
                "n_iterations": int(len(subset)),
                "n_protocols_observed": int(subset["protocol_group"].nunique()),
                "dominant_protocol": dominant_protocol,
                "protocol_switches": protocol_switches,
                "protocol_switch_rate": float(protocol_switch_rate),
                "protocol_eta_sq": float(protocol_eta_sq),
                "within_protocol_weighted_std": float(within_protocol_weighted_std),
                "protocol_gap_polyA_vs_total": protocol_gap,
                "polyA_wins": protocol_counts.get("polyA_plus_rna_seq", 0),
                "total_rna_wins": protocol_counts.get("total_rna_seq", 0),
                "other_wins": protocol_counts.get("other", 0),
            }
        )

        for _, row in protocol_stats.iterrows():
            pair_protocol_rows.append(
                {
                    "variant_id": variant_id,
                    "gene_id": gene_id,
                    "protocol_group": row["protocol_group"],
                    "n_iterations": int(row["n_iterations"]),
                    "mean_score": float(row["mean_score"]),
                    "std_score": float(row["std_score"]),
                }
            )

    pair_summary = pd.DataFrame(pair_summary_rows)
    pair_protocol_summary = pd.DataFrame(pair_protocol_rows)
    pair_summary = pair_summary.merge(
        agg[["variant_id", "gene_id", "mean_score", "std_score"]],
        on=["variant_id", "gene_id"],
        how="left",
    ).rename(
        columns={
            "mean_score": "overall_mean_score",
            "std_score": "overall_std_score",
        }
    )

    total_pairs = len(pair_summary)
    mixed_pairs = int((pair_summary["n_protocols_observed"] > 1).sum())
    always_poly = int((pair_summary["polyA_wins"] == pair_summary["n_iterations"]).sum())
    always_total = int((pair_summary["total_rna_wins"] == pair_summary["n_iterations"]).sum())
    top_variable = pair_summary.nlargest(min(20, total_pairs), "overall_std_score")
    top_variable_protocol_driven = int((top_variable["protocol_eta_sq"] >= 0.8).sum())

    stable_single = pair_summary[pair_summary["n_protocols_observed"] == 1]["overall_std_score"]
    mixed_std = pair_summary[pair_summary["n_protocols_observed"] > 1]["overall_std_score"]

    overview_lines = [
        "AlphaGenome Batch Stability: protocol-aware interpretation",
        f"Total variant/gene pairs: {total_pairs}",
        f"Pairs with >1 winning protocol across iterations: {mixed_pairs} ({mixed_pairs / max(total_pairs, 1) * 100:.1f}%)",
        f"Pairs always won by polyA: {always_poly}",
        f"Pairs always won by total RNA: {always_total}",
        f"Median overall std, single-protocol winners: {stable_single.median():.6f}" if not stable_single.empty else "Median overall std, single-protocol winners: NA",
        f"Median overall std, mixed-protocol winners: {mixed_std.median():.6f}" if not mixed_std.empty else "Median overall std, mixed-protocol winners: NA",
        f"Top variable pairs with protocol_eta_sq >= 0.8: {top_variable_protocol_driven}/{len(top_variable)}",
        (
            "Interpretation: this rerun is fully deterministic. There is no protocol switching, no track switching, and no residual score variability across the 50 repeated runs."
            if mixed_pairs == 0 and float(pair_summary["overall_std_score"].max()) == 0.0
            else "Interpretation: when a pair has high overall variability and high protocol_eta_sq, the instability is largely explained by winner switching between protocols under max-|score| aggregation rather than broad per-run noise."
        ),
    ]

    return pair_summary, pair_protocol_summary, iter_protocol_counts, overview_lines


def _load_protocol_iteration_frame(path: Path) -> pd.DataFrame:
    protocol_iterations = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    required = {"iteration", "variant_id", "gene_id", "protocol_group", "raw_score"}
    missing = required - set(protocol_iterations.columns)
    if missing:
        raise ValueError(f"protocol iteration outputs are missing required columns: {sorted(missing)}")

    protocol_iterations["variant_id"] = protocol_iterations["variant_id"].astype("string").str.strip()
    protocol_iterations["gene_id"] = protocol_iterations["gene_id"].astype("string").str.strip()
    protocol_iterations["protocol_group"] = protocol_iterations["protocol_group"].astype("string").str.strip()
    protocol_iterations["raw_score"] = pd.to_numeric(protocol_iterations["raw_score"], errors="coerce")
    protocol_iterations = protocol_iterations.dropna(
        subset=["iteration", "variant_id", "gene_id", "protocol_group", "raw_score"]
    ).copy()
    if protocol_iterations.empty:
        raise ValueError("protocol iteration outputs contain no usable rows after cleaning")
    return protocol_iterations


def _load_protocol_aggregate_frame(path: Path) -> pd.DataFrame:
    agg_protocol = pd.read_csv(path, sep="\t", compression="gzip", low_memory=False)
    required = {"variant_id", "gene_id", "protocol_group", "mean_score", "std_score"}
    missing = required - set(agg_protocol.columns)
    if missing:
        raise ValueError(f"protocol aggregated statistics are missing required columns: {sorted(missing)}")

    agg_protocol["variant_id"] = agg_protocol["variant_id"].astype("string").str.strip()
    agg_protocol["gene_id"] = agg_protocol["gene_id"].astype("string").str.strip()
    agg_protocol["protocol_group"] = agg_protocol["protocol_group"].astype("string").str.strip()
    agg_protocol["mean_score"] = pd.to_numeric(agg_protocol["mean_score"], errors="coerce")
    agg_protocol["std_score"] = pd.to_numeric(agg_protocol["std_score"], errors="coerce").fillna(0.0)
    agg_protocol["abs_mean_score"] = agg_protocol["mean_score"].abs()
    return agg_protocol.dropna(subset=["variant_id", "gene_id", "protocol_group", "mean_score", "std_score"]).copy()


def _summarize_protocol_separated(
    agg_protocol: pd.DataFrame,
    protocol_iterations: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    counts = (
        agg_protocol.groupby("protocol_group")
        .agg(
            n_series=("variant_id", "size"),
            median_std=("std_score", "median"),
            max_std=("std_score", "max"),
            median_abs_mean=("abs_mean_score", "median"),
        )
        .reset_index()
    )

    lines = [
        "AlphaGenome Batch Stability: protocol-separated interpretation",
        f"Total protocol-specific series: {len(agg_protocol)}",
    ]

    for _, row in counts.iterrows():
        lines.append(
            f"{row['protocol_group']}: n={int(row['n_series'])}, median std={row['median_std']:.6f}, "
            f"max std={row['max_std']:.6f}, median |mean|={row['median_abs_mean']:.6f}"
        )

    if "track_switches" in agg_protocol.columns:
        for protocol, sub in agg_protocol.groupby("protocol_group"):
            lines.append(
                f"{protocol}: series with >0 within-protocol track switches = "
                f"{int((sub['track_switches'] > 0).sum())}/{len(sub)}"
            )

    if float(agg_protocol["std_score"].max()) == 0.0:
        lines.append(
            "Interpretation: after separating polyA and total RNA, this batch is fully deterministic within each protocol."
        )
    else:
        lines.append(
            "Interpretation: any remaining variability here reflects within-protocol instability, not polyA-versus-total-RNA mixing."
        )

    protocol_counts_by_iteration = (
        protocol_iterations.groupby(["iteration", "protocol_group"])
        .size()
        .rename("series_count")
        .reset_index()
        .pivot(index="iteration", columns="protocol_group", values="series_count")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    return protocol_counts_by_iteration, lines


def _plot_overview(agg: pd.DataFrame, plot_dir: Path) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    sns.histplot(agg["mean_score"], bins=24, kde=True, ax=axes[0, 0], color="#1f77b4")
    axes[0, 0].axvline(0.0, color="black", linestyle="--", linewidth=1)
    axes[0, 0].set_title("Distribution of Mean Scores")
    axes[0, 0].set_xlabel("Mean score across 50 iterations")

    sns.histplot(agg["std_score"], bins=24, kde=True, ax=axes[0, 1], color="#d62728")
    axes[0, 1].set_title("Distribution of Stability Std")
    axes[0, 1].set_xlabel("Standard deviation across 50 iterations")

    sns.scatterplot(
        data=agg,
        x="mean_score",
        y="std_score",
        hue="abs_mean_score",
        palette="viridis",
        s=70,
        ax=axes[1, 0],
    )
    axes[1, 0].axvline(0.0, color="black", linestyle="--", linewidth=1)
    axes[1, 0].set_title("Effect Size vs Variability")
    axes[1, 0].set_xlabel("Mean score")
    axes[1, 0].set_ylabel("Std score")
    axes[1, 0].legend(title="|mean|", loc="upper right", frameon=True)

    sns.boxplot(x=agg["std_score"], ax=axes[1, 1], color="#ff9896")
    sns.stripplot(x=agg["std_score"], ax=axes[1, 1], color="#8c564b", size=4, alpha=0.55)
    axes[1, 1].set_title("Spread of Variant Stability")
    axes[1, 1].set_xlabel("Std score")

    fig.suptitle("AlphaGenome Batch Stability Overview", fontsize=16)
    fig.tight_layout()

    out_path = plot_dir / "01_stability_overview.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_locus_map(agg: pd.DataFrame, plot_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(14, 6))
    std_max = float(agg["std_score"].max())
    size_scale = 1.0 if std_max <= 0 else std_max

    scatter = ax.scatter(
        agg["position"],
        agg["mean_score"],
        c=agg["std_score"],
        s=80 + 220 * (agg["std_score"] / size_scale),
        cmap="magma",
        alpha=0.9,
        edgecolors="white",
        linewidths=0.5,
    )
    ax.axhline(0.0, color="black", linestyle="--", linewidth=1)
    ax.set_title("Variant Position vs Mean Score")
    ax.set_xlabel("Genomic position")
    ax.set_ylabel("Mean score across 50 iterations")

    top_labels = agg.nlargest(8, "std_score")
    for _, row in top_labels.iterrows():
        ax.annotate(
            _short_variant_label(row["variant_id"], max_change_len=10),
            (row["position"], row["mean_score"]),
            textcoords="offset points",
            xytext=(5, 5),
            fontsize=8,
        )

    cbar = fig.colorbar(scatter, ax=ax)
    cbar.set_label("Std score")

    fig.tight_layout()
    out_path = plot_dir / "02_locus_score_map.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_top_variable(agg: pd.DataFrame, plot_dir: Path, top_n: int = 15) -> Path:
    top = agg.nlargest(top_n, "std_score").copy()
    top = top.sort_values("std_score", ascending=True)
    top["label"] = top["variant_id"].map(lambda value: _short_variant_label(value, max_change_len=10))

    fig, ax = plt.subplots(figsize=(12, 8))
    bars = ax.barh(top["label"], top["std_score"], color="#ff7f0e", alpha=0.9)
    ax.set_title(f"Top {top_n} Most Variable Variants")
    ax.set_xlabel("Std score across 50 iterations")
    ax.set_ylabel("Variant")

    for bar, (_, row) in zip(bars, top.iterrows()):
        ax.text(
            bar.get_width() + 0.0001,
            bar.get_y() + bar.get_height() / 2,
            f"mean={row['mean_score']:.4f}",
            va="center",
            fontsize=8,
        )

    fig.tight_layout()
    out_path = plot_dir / "03_top_variable_variants.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_trajectories(
    agg: pd.DataFrame,
    iterations: pd.DataFrame,
    plot_dir: Path,
    top_k: int,
) -> Path:
    top = agg.nlargest(top_k, "std_score")[["variant_id", "gene_id", "mean_score", "std_score"]].copy()
    top_ids = set(zip(top["variant_id"], top["gene_id"]))
    traj = iterations[iterations[["variant_id", "gene_id"]].apply(tuple, axis=1).isin(top_ids)].copy()
    traj = traj.merge(top, on=["variant_id", "gene_id"], how="left", suffixes=("", "_agg"))

    n_cols = 2
    n_rows = (len(top) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(14, 4.2 * n_rows), sharex=True)
    axes_flat = list(axes.flat) if hasattr(axes, "flat") else [axes]

    for ax, (_, row) in zip(axes_flat, top.iterrows()):
        subset = traj[(traj["variant_id"] == row["variant_id"]) & (traj["gene_id"] == row["gene_id"])].copy()
        subset = subset.sort_values("iteration")
        ax.plot(subset["iteration"], subset["raw_score"], color="#7f7f7f", linewidth=1.3, alpha=0.7)
        if "protocol_group" in subset.columns:
            palette = {
                "polyA_plus_rna_seq": "#1f77b4",
                "total_rna_seq": "#d62728",
                "other": "#2ca02c",
            }
            for protocol, protocol_df in subset.groupby("protocol_group", sort=False):
                ax.scatter(
                    protocol_df["iteration"],
                    protocol_df["raw_score"],
                    s=28,
                    color=palette.get(protocol, "#2ca02c"),
                    label=protocol,
                    zorder=3,
                )
            ax.legend(loc="best", fontsize=8, frameon=True)
        else:
            ax.scatter(subset["iteration"], subset["raw_score"], s=28, color="#1f77b4", zorder=3)
        ax.axhline(row["mean_score"], color="#d62728", linestyle="--", linewidth=1)
        ax.set_title(
            f"{_short_variant_label(row['variant_id'], max_change_len=12)}\n"
            f"mean={row['mean_score']:.4f}, std={row['std_score']:.4f}",
            fontsize=10,
        )
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Raw score")

    for ax in axes_flat[len(top):]:
        ax.axis("off")

    fig.suptitle("Iteration-by-Iteration Trajectories for Most Variable Variants", fontsize=16)
    fig.tight_layout()
    out_path = plot_dir / "04_top_variable_trajectories.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_protocol_vs_variability(pair_summary: pd.DataFrame, plot_dir: Path) -> Path:
    fig, ax = plt.subplots(figsize=(12, 7))
    plot_df = pair_summary.copy()
    plot_df["protocol_mode"] = np.where(
        plot_df["n_protocols_observed"] > 1,
        "mixed winners",
        "single winner",
    )

    sns.scatterplot(
        data=plot_df,
        x="protocol_eta_sq",
        y="overall_std_score",
        hue="protocol_mode",
        style="dominant_protocol",
        s=90,
        ax=ax,
    )
    ax.set_title("How Much Variability Is Explained by Winner Protocol")
    ax.set_xlabel("Protocol eta^2 (between-protocol share of score variance)")
    ax.set_ylabel("Overall std score across iterations")
    ax.set_xlim(-0.02, 1.02)
    ax.legend(loc="best", frameon=True)

    fig.tight_layout()
    out_path = plot_dir / "05_protocol_vs_variability.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_protocol_counts(iter_protocol_counts: pd.DataFrame, plot_dir: Path) -> Path:
    plot_df = iter_protocol_counts.melt(
        id_vars="iteration",
        var_name="protocol_group",
        value_name="winner_count",
    )

    fig, ax = plt.subplots(figsize=(12, 6))
    sns.lineplot(
        data=plot_df,
        x="iteration",
        y="winner_count",
        hue="protocol_group",
        marker="o",
        linewidth=2,
        ax=ax,
    )
    ax.set_title("Winning Protocol Counts Across Iterations")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Number of winner rows")
    ax.legend(loc="best", frameon=True, title="Protocol")

    fig.tight_layout()
    out_path = plot_dir / "06_protocol_winner_counts.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_protocol_separated_overview(agg_protocol: pd.DataFrame, plot_dir: Path) -> Path:
    fig, axes = plt.subplots(2, 2, figsize=(16, 11))

    sns.histplot(
        data=agg_protocol,
        x="mean_score",
        hue="protocol_group",
        bins=24,
        kde=True,
        ax=axes[0, 0],
        element="step",
        common_norm=False,
    )
    axes[0, 0].axvline(0.0, color="black", linestyle="--", linewidth=1)
    axes[0, 0].set_title("Mean Scores by Protocol")
    axes[0, 0].set_xlabel("Mean score across 50 iterations")

    sns.histplot(
        data=agg_protocol,
        x="std_score",
        hue="protocol_group",
        bins=24,
        kde=False,
        ax=axes[0, 1],
        element="step",
        common_norm=False,
    )
    axes[0, 1].set_title("Stability Std by Protocol")
    axes[0, 1].set_xlabel("Standard deviation across 50 iterations")

    sns.scatterplot(
        data=agg_protocol,
        x="mean_score",
        y="std_score",
        hue="protocol_group",
        size="abs_mean_score",
        sizes=(40, 180),
        ax=axes[1, 0],
    )
    axes[1, 0].axvline(0.0, color="black", linestyle="--", linewidth=1)
    axes[1, 0].set_title("Effect Size vs Variability by Protocol")
    axes[1, 0].set_xlabel("Mean score")
    axes[1, 0].set_ylabel("Std score")

    sns.boxplot(data=agg_protocol, x="protocol_group", y="std_score", ax=axes[1, 1], color="#f4b6c2")
    sns.stripplot(
        data=agg_protocol,
        x="protocol_group",
        y="std_score",
        ax=axes[1, 1],
        color="#8c564b",
        size=4,
        alpha=0.55,
    )
    axes[1, 1].set_title("Spread of Stability Within Each Protocol")
    axes[1, 1].set_xlabel("Protocol")
    axes[1, 1].set_ylabel("Std score")
    axes[1, 1].tick_params(axis="x", rotation=12)

    fig.suptitle("AlphaGenome Batch Stability Overview by Protocol", fontsize=16)
    fig.tight_layout()

    out_path = plot_dir / "01_protocol_stability_overview.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_protocol_separated_trajectories(
    agg_protocol: pd.DataFrame,
    protocol_iterations: pd.DataFrame,
    plot_dir: Path,
    top_k: int,
) -> Path:
    protocols = list(agg_protocol["protocol_group"].dropna().unique())
    per_protocol = max(1, top_k // max(len(protocols), 1))
    top_frames: list[pd.DataFrame] = []
    for protocol in protocols:
        sub = agg_protocol[agg_protocol["protocol_group"] == protocol].copy()
        if sub.empty:
            continue
        top_frames.append(sub.nlargest(per_protocol, "std_score"))

    top = pd.concat(top_frames, ignore_index=True) if top_frames else agg_protocol.nlargest(top_k, "std_score").copy()
    if len(top) < top_k:
        extras = agg_protocol.loc[
            ~agg_protocol.set_index(["variant_id", "gene_id", "protocol_group"]).index.isin(
                top.set_index(["variant_id", "gene_id", "protocol_group"]).index
            )
        ].nlargest(top_k - len(top), "std_score")
        top = pd.concat([top, extras], ignore_index=True)

    keys = set(zip(top["variant_id"], top["gene_id"], top["protocol_group"]))
    traj = protocol_iterations[
        protocol_iterations[["variant_id", "gene_id", "protocol_group"]].apply(tuple, axis=1).isin(keys)
    ].copy()
    traj = traj.merge(
        top[["variant_id", "gene_id", "protocol_group", "mean_score", "std_score"]],
        on=["variant_id", "gene_id", "protocol_group"],
        how="left",
    )

    n_cols = 2
    n_rows = (len(top) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4.5 * n_rows), sharex=True)
    axes_flat = list(axes.flat) if hasattr(axes, "flat") else [axes]

    for ax, (_, row) in zip(axes_flat, top.iterrows()):
        subset = traj[
            (traj["variant_id"] == row["variant_id"])
            & (traj["gene_id"] == row["gene_id"])
            & (traj["protocol_group"] == row["protocol_group"])
        ].copy().sort_values("iteration")
        ax.plot(subset["iteration"], subset["raw_score"], color="#7f7f7f", linewidth=1.2, alpha=0.7)
        color_col = "track_key" if "track_key" in subset.columns else "protocol_group"
        if color_col in subset.columns:
            for label, label_df in subset.groupby(color_col, sort=False):
                ax.scatter(label_df["iteration"], label_df["raw_score"], s=28, label=str(label), zorder=3)
            if subset[color_col].nunique() <= 4:
                ax.legend(loc="best", fontsize=7, frameon=True)
        else:
            ax.scatter(subset["iteration"], subset["raw_score"], s=28, zorder=3)

        track_switch_note = ""
        if "track_switches" in row.index and pd.notna(row.get("track_switches")):
            track_switch_note = f", switches={int(row['track_switches'])}"
        ax.axhline(row["mean_score"], color="#d62728", linestyle="--", linewidth=1)
        ax.set_title(
            f"{row['protocol_group']}: {_short_variant_label(row['variant_id'], max_change_len=12)}\n"
            f"mean={row['mean_score']:.4f}, std={row['std_score']:.4f}{track_switch_note}",
            fontsize=10,
        )
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Raw score")

    for ax in axes_flat[len(top):]:
        ax.axis("off")

    fig.suptitle("Top Variable Trajectories Within Each Protocol", fontsize=16)
    fig.tight_layout()
    out_path = plot_dir / "02_protocol_top_variable_trajectories.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_within_protocol_track_switches(agg_protocol: pd.DataFrame, plot_dir: Path) -> Path:
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))

    if "n_distinct_selected_tracks" in agg_protocol.columns:
        sns.boxplot(
            data=agg_protocol,
            x="protocol_group",
            y="n_distinct_selected_tracks",
            ax=axes[0],
            color="#a1d99b",
        )
        sns.stripplot(
            data=agg_protocol,
            x="protocol_group",
            y="n_distinct_selected_tracks",
            ax=axes[0],
            color="#2ca25f",
            size=4,
            alpha=0.6,
        )
        axes[0].set_title("Distinct Selected Tracks per Series")
        axes[0].set_xlabel("Protocol")
        axes[0].set_ylabel("Count")
        axes[0].tick_params(axis="x", rotation=12)
    else:
        axes[0].text(0.5, 0.5, "Track identity not available", ha="center", va="center")
        axes[0].axis("off")

    if "track_switches" in agg_protocol.columns:
        sns.boxplot(
            data=agg_protocol,
            x="protocol_group",
            y="track_switches",
            ax=axes[1],
            color="#9ecae1",
        )
        sns.stripplot(
            data=agg_protocol,
            x="protocol_group",
            y="track_switches",
            ax=axes[1],
            color="#3182bd",
            size=4,
            alpha=0.6,
        )
        axes[1].set_title("Within-Protocol Track Switches Across Iterations")
        axes[1].set_xlabel("Protocol")
        axes[1].set_ylabel("Switch count")
        axes[1].tick_params(axis="x", rotation=12)
    else:
        axes[1].text(0.5, 0.5, "Track switch summary not available", ha="center", va="center")
        axes[1].axis("off")

    fig.suptitle("Within-Protocol Track Competition Summary", fontsize=16)
    fig.tight_layout()
    out_path = plot_dir / "03_within_protocol_track_switches.png"
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main() -> None:
    args = _parse_args()
    results_dir = _resolve_results_dir(args.results_dir)
    output_dir = args.output_dir.resolve() if args.output_dir is not None else results_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    agg_path = results_dir / "aggregated_statistics.tsv.gz"
    agg_protocol_path = results_dir / "aggregated_statistics_by_protocol.tsv.gz"
    iterations_dir = results_dir / "iterations"
    protocol_iterations_path = results_dir / "protocol_iteration_scores.tsv.gz"
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    if agg_protocol_path.exists() and protocol_iterations_path.exists():
        agg_protocol = _load_protocol_aggregate_frame(agg_protocol_path)
        protocol_iterations = _load_protocol_iteration_frame(protocol_iterations_path)
        iter_protocol_counts, overview_lines = _summarize_protocol_separated(agg_protocol, protocol_iterations)

        overview_path = output_dir / "protocol_separated_overview.txt"
        iter_protocol_path = output_dir / "protocol_iteration_counts.tsv.gz"
        overview_path.write_text("\n".join(overview_lines) + "\n")
        iter_protocol_counts.to_csv(iter_protocol_path, sep="\t", index=False, compression="gzip")

        print("Protocol-separated interpretation:")
        for line in overview_lines:
            print(f"  {line}")

        print("Saved protocol summaries:")
        for path in [iter_protocol_path, overview_path]:
            print(path)

        if not HAS_PLOTTING:
            print("Plotting skipped because matplotlib/seaborn are not installed in this environment.")
            return

        sns.set_theme(style="whitegrid", context="talk")
        outputs = [
            _plot_protocol_separated_overview(agg_protocol, plot_dir),
            _plot_protocol_separated_trajectories(agg_protocol, protocol_iterations, plot_dir, top_k=args.top_k),
            _plot_within_protocol_track_switches(agg_protocol, plot_dir),
            _plot_protocol_counts(iter_protocol_counts, plot_dir),
        ]
    else:
        agg = pd.read_csv(agg_path, sep="\t", compression="gzip", low_memory=False)
        required = {"variant_id", "gene_id", "mean_score", "std_score"}
        missing = required - set(agg.columns)
        if missing:
            raise ValueError(f"aggregated statistics are missing required columns: {sorted(missing)}")

        parsed = agg["variant_id"].map(_parse_variant_id)
        agg["chromosome"] = parsed.map(lambda item: item[0])
        agg["position"] = parsed.map(lambda item: item[1])
        agg["change"] = parsed.map(lambda item: item[2])
        agg["abs_mean_score"] = agg["mean_score"].abs()
        agg = agg.sort_values("position").reset_index(drop=True)

        iterations = _load_iteration_frame(iterations_dir)
        pair_summary, pair_protocol_summary, iter_protocol_counts, overview_lines = _summarize_protocol_competition(agg, iterations)

        pair_summary_path = output_dir / "protocol_competition_summary.tsv.gz"
        pair_protocol_path = output_dir / "protocol_pair_statistics.tsv.gz"
        iter_protocol_path = output_dir / "iteration_protocol_winner_counts.tsv.gz"
        overview_path = output_dir / "protocol_competition_overview.txt"

        pair_summary.to_csv(pair_summary_path, sep="\t", index=False, compression="gzip")
        pair_protocol_summary.to_csv(pair_protocol_path, sep="\t", index=False, compression="gzip")
        iter_protocol_counts.to_csv(iter_protocol_path, sep="\t", index=False, compression="gzip")
        overview_path.write_text("\n".join(overview_lines) + "\n")

        print("Protocol-aware interpretation (from globally selected rows only):")
        print("  Note: protocol-separated iteration files were not found, so this fallback still reflects cross-protocol max-abs selection.")
        for line in overview_lines:
            print(f"  {line}")

        print("Saved protocol summaries:")
        for path in [pair_summary_path, pair_protocol_path, iter_protocol_path, overview_path]:
            print(path)

        if not HAS_PLOTTING:
            print("Plotting skipped because matplotlib/seaborn are not installed in this environment.")
            return

        sns.set_theme(style="whitegrid", context="talk")
        outputs = [
            _plot_overview(agg, plot_dir),
            _plot_locus_map(agg, plot_dir),
            _plot_top_variable(agg, plot_dir),
            _plot_trajectories(agg, iterations, plot_dir, top_k=args.top_k),
            _plot_protocol_vs_variability(pair_summary, plot_dir),
            _plot_protocol_counts(iter_protocol_counts, plot_dir),
        ]

    print("Saved plots:")
    for path in outputs:
        print(path)


if __name__ == "__main__":
    main()
