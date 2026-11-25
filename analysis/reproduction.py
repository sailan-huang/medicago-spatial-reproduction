from __future__ import annotations

import hashlib
import json
import logging
import platform
import subprocess
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.stats import kruskal, mannwhitneyu
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.neighbors import NearestNeighbors

LOGGER = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]
PAPER_DPI = 300

TISSUE_COLORS = {
    "Cortex_1": "#4C78A8",
    "Cortex_2": "#9ECAE9",
    "Vascular bundle": "#8C8C8C",
    "Meristematic zone": "#E45756",
    "Proximal infection zone": "#F2A93B",
    "Distal infection zone": "#B279A2",
    "Interzone": "#72B7B2",
    "Nitrogen fixation zone": "#54A24B",
    "Senescence zone": "#B07AA1",
}
TISSUE_ABBREVIATIONS = {
    "Cortex_1": "C1",
    "Cortex_2": "C2",
    "Vascular bundle": "VB",
    "Meristematic zone": "MZ",
    "Proximal infection zone": "PIZ",
    "Distal infection zone": "DIZ",
    "Interzone": "IZ",
    "Nitrogen fixation zone": "NFZ",
    "Senescence zone": "SZ",
}


@dataclass(frozen=True)
class Paths:
    root: Path = ROOT

    @property
    def config(self) -> Path:
        return self.root / "config" / "analysis.yaml"

    @property
    def figure_dir(self) -> Path:
        return self.root / "results" / "figures"

    @property
    def table_dir(self) -> Path:
        return self.root / "results" / "tables"


def configure_logging() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)


def configure_publication_style() -> None:
    sns.set_theme(
        style="ticks",
        context="paper",
        font="DejaVu Sans",
        rc={
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.title_fontsize": 9,
            "figure.dpi": 120,
            "savefig.dpi": PAPER_DPI,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        },
    )


def save_figure(fig: plt.Figure, png_path: Path) -> None:
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=PAPER_DPI, bbox_inches="tight", facecolor="white")
    fig.savefig(
        png_path.with_suffix(".pdf"),
        bbox_inches="tight",
        facecolor="white",
        metadata={"CreationDate": None, "ModDate": None},
    )
    fig.savefig(
        png_path.with_suffix(".svg"),
        bbox_inches="tight",
        facecolor="white",
        metadata={"Date": None},
    )


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        va="top",
        ha="left",
        clip_on=False,
    )


def add_scale_bar(
    ax: plt.Axes,
    spatial: np.ndarray,
    config: dict[str, Any],
) -> None:
    publication = config["analysis"]["publication"]
    scale_um = float(publication["scale_bar_um"])
    um_per_unit = float(publication["spatial_um_per_coordinate"])
    scale_units = scale_um / um_per_unit
    x_min, x_max = np.nanmin(spatial[:, 0]), np.nanmax(spatial[:, 0])
    y_min, y_max = np.nanmin(spatial[:, 1]), np.nanmax(spatial[:, 1])
    x0 = x_min + 0.03 * (x_max - x_min)
    y0 = y_max - 0.035 * (y_max - y_min)
    ax.plot(
        [x0, x0 + scale_units],
        [y0, y0],
        color="white",
        linewidth=5,
        solid_capstyle="butt",
        zorder=10,
    )
    ax.plot(
        [x0, x0 + scale_units],
        [y0, y0],
        color="black",
        linewidth=2,
        solid_capstyle="butt",
        zorder=11,
    )
    ax.text(
        x0 + scale_units / 2,
        y0 - 0.012 * (y_max - y_min),
        f"{scale_um:g} µm",
        ha="center",
        va="bottom",
        fontsize=8,
        color="black",
        zorder=12,
    )


def style_spatial_axes(ax: plt.Axes, title: str) -> None:
    ax.set_title(title)
    ax.set_xlabel("x (bin coordinate)")
    ax.set_ylabel("y (bin coordinate)")
    ax.set_aspect("equal")
    ax.grid(False)
    ax.tick_params(length=2)


def safe_git_revision() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def git_worktree_dirty() -> bool | None:
    try:
        status = subprocess.check_output(
            ["git", "status", "--porcelain"],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        )
        return bool(status.strip())
    except (OSError, subprocess.CalledProcessError):
        return None


def load_config(path: Path | None = None) -> dict[str, Any]:
    import yaml

    path = path or Paths().config
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_inputs(config: dict[str, Any]) -> tuple[Path, Path]:
    paths = Paths()
    sample = paths.root / "data" / "raw" / config["data"]["sample_filename"]
    annotation = (
        paths.root / "data" / "external" / config["data"]["annotation_filename"]
    )
    for label, path, expected in (
        ("sample", sample, config["data"]["sample_sha256"]),
        ("annotation", annotation, config["data"]["annotation_sha256"]),
    ):
        if not path.exists():
            raise FileNotFoundError(
                f"Missing {label}: {path}. Run `python scripts/download_data.py`."
            )
        observed = sha256(path)
        if observed.lower() != expected.lower():
            raise ValueError(
                f"SHA-256 mismatch for {label}: expected {expected}, observed {observed}"
            )
    return sample, annotation


def resolve_supplementary(config: dict[str, Any]) -> Path:
    path = Paths().root / "data" / "external" / config["data"]["supplementary_filename"]
    expected = config["data"]["supplementary_sha256"]
    if not path.exists():
        raise FileNotFoundError(
            f"Missing supplementary table: {path}. "
            "Run `python scripts/download_data.py`."
        )
    observed = sha256(path)
    if observed.lower() != expected.lower():
        raise ValueError(
            "SHA-256 mismatch for supplementary table: "
            f"expected {expected}, observed {observed}"
        )
    return path


def classify_gene(gene_id: str) -> str:
    if gene_id.startswith("MtrunA17"):
        return "Medicago"
    if gene_id.startswith("SM-"):
        return "Sinorhizobium"
    return "unknown"


def load_and_validate(sample_path: Path) -> ad.AnnData:
    adata = sc.read_h5ad(sample_path)
    required_obs = {"organ", "phase", "group", "array_row", "array_col"}
    missing_obs = required_obs.difference(adata.obs.columns)
    if missing_obs:
        raise ValueError(f"Missing required observation fields: {sorted(missing_obs)}")
    if "counts" not in adata.layers:
        raise ValueError("The input must contain a raw-count layer named 'counts'.")
    if "spatial" not in adata.obsm:
        raise ValueError(
            "The input must contain spatial coordinates in obsm['spatial']."
        )
    if not adata.var_names.is_unique:
        raise ValueError("Gene identifiers are not unique.")
    spatial = np.asarray(adata.obsm["spatial"])
    if spatial.shape != (adata.n_obs, 2):
        raise ValueError(
            f"Spatial coordinates must have shape ({adata.n_obs}, 2), got {spatial.shape}."
        )

    # The public file stores counts in both .X and layers['counts'].
    adata.X = adata.layers["counts"].copy()
    if not sparse.issparse(adata.X):
        adata.X = sparse.csr_matrix(adata.X)
    if adata.X.data.size and adata.X.data.min() < 0:
        raise ValueError("The raw-count layer contains negative values.")
    if adata.X.data.size and not np.allclose(adata.X.data, np.round(adata.X.data)):
        raise ValueError("The raw-count layer contains non-integer values.")
    adata.var["species"] = [classify_gene(g) for g in adata.var_names]
    if (adata.var["species"] == "unknown").any():
        examples = adata.var_names[adata.var["species"] == "unknown"][:10].tolist()
        raise ValueError(f"Unclassified gene identifiers found: {examples}")
    if not np.isfinite(spatial).all():
        raise ValueError("Spatial coordinates contain non-finite values.")
    if adata.obs["organ"].isna().any():
        raise ValueError("Tissue annotations contain missing values.")
    recomputed_counts = np.asarray(adata.X.sum(axis=1)).ravel()
    recomputed_features = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    if "nCount_Spatial" in adata.obs and not np.allclose(
        recomputed_counts, adata.obs["nCount_Spatial"].to_numpy()
    ):
        raise ValueError("Stored nCount_Spatial values do not match the count layer.")
    if "nFeature_Spatial" in adata.obs and not np.allclose(
        recomputed_features, adata.obs["nFeature_Spatial"].to_numpy()
    ):
        raise ValueError("Stored nFeature_Spatial values do not match the count layer.")
    if "array_row" in adata.obs and not np.allclose(
        spatial[:, 0], adata.obs["array_row"].to_numpy()
    ):
        raise ValueError("Spatial coordinate column 1 does not match array_row.")
    if "array_col" in adata.obs and not np.allclose(
        spatial[:, 1], adata.obs["array_col"].to_numpy()
    ):
        raise ValueError("Spatial coordinate column 2 does not match array_col.")
    return adata


def add_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    plant_mask = (adata.var["species"] == "Medicago").to_numpy()
    rhizobia_mask = (adata.var["species"] == "Sinorhizobium").to_numpy()
    total = np.asarray(adata.X.sum(axis=1)).ravel()
    plant = np.asarray(adata.X[:, plant_mask].sum(axis=1)).ravel()
    rhizobia = np.asarray(adata.X[:, rhizobia_mask].sum(axis=1)).ravel()
    detected = np.asarray((adata.X > 0).sum(axis=1)).ravel()

    adata.obs["total_counts_recomputed"] = total
    adata.obs["detected_genes_recomputed"] = detected
    adata.obs["plant_counts"] = plant
    adata.obs["rhizobia_counts"] = rhizobia
    adata.obs["rhizobia_fraction"] = np.divide(
        rhizobia,
        total,
        out=np.full(total.shape, np.nan, dtype=float),
        where=total > 0,
    )
    adata.obs["log10_total_counts_recomputed"] = np.log10(total + 1)
    return adata


def load_marker_mapping(annotation_path: Path) -> pd.DataFrame:
    mapping = pd.read_excel(annotation_path, header=2)
    mapping.columns = [
        "mt5_locus",
        "acronym",
        "mt5_r19_new",
        "mt5_r18_new",
        "new_locus_tag",
        "description",
        "publications",
        "mt5_r16",
        "mt4",
        "mt20120830",
        "affy",
        "nimblegen",
    ]
    mapping = mapping[["mt5_locus", "acronym", "mt4", "description"]].copy()
    mapping = mapping.dropna(subset=["mt5_locus", "acronym"])
    mapping["gene_id"] = (
        mapping["mt5_locus"].astype(str).str.replace("_", "", regex=False).str.strip()
    )
    mapping["acronym"] = mapping["acronym"].astype(str)
    return mapping


def load_supplementary_marker_evidence(path: Path) -> pd.DataFrame:
    evidence = pd.read_excel(
        path,
        sheet_name="S1-5",
        skiprows=4,
        header=None,
        usecols=[0, 1, 2, 3, 22],
    )
    evidence.columns = [
        "cluster",
        "annotation_source",
        "gene_id",
        "gene_symbol",
        "functional_annotation",
    ]
    evidence = evidence.dropna(subset=["cluster", "gene_id"]).copy()
    evidence["gene_id_norm"] = (
        evidence["gene_id"].astype(str).str.replace("_", "", regex=False).str.strip()
    )
    evidence["cluster"] = evidence["cluster"].astype(str).str.strip()
    evidence["gene_symbol"] = evidence["gene_symbol"].astype(str).str.strip()
    return evidence


def load_sample_reference_metadata(path: Path, sample_id: str) -> pd.DataFrame:
    metadata = pd.read_excel(path, sheet_name="S1-3", header=2)
    metadata = metadata[
        metadata["Sample/lasso ID"].astype(str).str.strip() == sample_id
    ].copy()
    if len(metadata) != 1:
        raise ValueError(
            f"Expected one Table S1-3 row for {sample_id}, found {len(metadata)}."
        )
    return metadata


def resolve_markers(
    marker_config: dict[str, list[str]],
    preferred_loci: dict[str, str],
    mapping: pd.DataFrame,
    available_genes: set[str],
    supplementary_evidence: pd.DataFrame | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for zone, acronyms in marker_config.items():
        for acronym in acronyms:
            preferred_locus = str(preferred_loci.get(acronym, "")).replace("_", "")
            preferred_rows = mapping[mapping["gene_id"] == preferred_locus].copy()
            if not preferred_rows.empty:
                best = preferred_rows.iloc[0]
                rows.append(
                    _marker_row(
                        zone=zone,
                        acronym=acronym,
                        gene_id=str(best["gene_id"]),
                        matched_acronym=str(best["acronym"]),
                        description=str(best["description"]),
                        available=str(best["gene_id"]) in available_genes,
                        supplementary_evidence=supplementary_evidence,
                    )
                )
                continue
            candidates = mapping[
                mapping["acronym"]
                .str.split("|", regex=False)
                .apply(lambda values, target=acronym: target in values)
            ].copy()
            if candidates.empty:
                rows.append(
                    _marker_row(
                        zone=zone,
                        acronym=acronym,
                        gene_id="",
                        matched_acronym="",
                        description="",
                        available=False,
                        supplementary_evidence=supplementary_evidence,
                    )
                )
                continue
            # Resolve each alias to the first Mt5.0 locus present in this sample.
            candidates["available"] = candidates["gene_id"].isin(available_genes)
            candidates["exact_alias"] = candidates["acronym"].astype(str) == acronym
            candidates = candidates.sort_values(
                ["available", "exact_alias"],
                ascending=[False, False],
            )
            best = candidates.iloc[0]
            rows.append(
                _marker_row(
                    zone=zone,
                    acronym=acronym,
                    gene_id=str(best["gene_id"]),
                    matched_acronym=str(best["acronym"]),
                    description=str(best["description"]),
                    available=bool(best["available"]),
                    supplementary_evidence=supplementary_evidence,
                )
            )
    return pd.DataFrame(rows)


def _marker_row(
    *,
    zone: str,
    acronym: str,
    gene_id: str,
    matched_acronym: str,
    description: str,
    available: bool,
    supplementary_evidence: pd.DataFrame | None,
) -> dict[str, str]:
    expected_cluster = {
        "meristematic_zone": "Meristematic zone",
        "proximal_infection_zone": "Proximal infection zone",
        "distal_infection_zone": "Distal infection zone",
        "interzone": "Interzone",
        "nitrogen_fixation_zone": "Nitrogen fixation zone",
        "vascular_bundle": "Vascular bundle",
    }[zone]
    evidence_row = pd.DataFrame()
    if supplementary_evidence is not None:
        evidence_row = supplementary_evidence[
            (supplementary_evidence["gene_id_norm"] == gene_id)
            & (supplementary_evidence["cluster"] == expected_cluster)
        ]
    evidence_match = bool(
        not evidence_row.empty and acronym in set(evidence_row["gene_symbol"])
    )
    return {
        "zone": zone,
        "requested_acronym": acronym,
        "resolved_gene_id": gene_id,
        "matched_acronym": matched_acronym,
        "description": description,
        "available_in_sample": str(available),
        "supplementary_evidence_match": str(evidence_match),
        "supplementary_annotation": (
            str(evidence_row.iloc[0]["functional_annotation"])
            if not evidence_row.empty
            else ""
        ),
    }


def plot_qc_overview(adata: ad.AnnData, path: Path) -> None:
    spatial = np.asarray(adata.obsm["spatial"])
    config = load_config()
    fig, axes = plt.subplots(2, 3, figsize=(15, 9), constrained_layout=True)

    sns.histplot(
        np.log10(adata.obs["total_counts_recomputed"] + 1),
        bins=35,
        ax=axes[0, 0],
        color="#35618f",
    )
    axes[0, 0].set(
        xlabel="log10(total counts + 1)",
        ylabel="Spatial bins",
    )

    sns.histplot(
        adata.obs["detected_genes_recomputed"],
        bins=35,
        ax=axes[0, 1],
        color="#c87941",
    )
    axes[0, 1].set(
        xlabel="Detected genes",
        ylabel="Spatial bins",
    )

    plant_rhizo = axes[0, 2].scatter(
        np.log10(adata.obs["plant_counts"] + 1),
        np.log10(adata.obs["rhizobia_counts"] + 1),
        c=adata.obs["rhizobia_fraction"],
        cmap="magma",
        s=16,
        alpha=0.8,
        vmin=0,
        vmax=np.nanquantile(adata.obs["rhizobia_fraction"], 0.99),
    )
    axes[0, 2].set(
        xlabel="log10(plant counts + 1)",
        ylabel="log10(rhizobial counts + 1)",
    )
    fig.colorbar(plant_rhizo, ax=axes[0, 2], shrink=0.8, label="Rhizobia / total")

    total_points = axes[1, 0].scatter(
        spatial[:, 0],
        spatial[:, 1],
        c=adata.obs["log10_total_counts_recomputed"],
        cmap="viridis",
        s=14,
        marker="s",
        linewidths=0,
        rasterized=True,
    )
    axes[1, 0].invert_yaxis()
    style_spatial_axes(axes[1, 0], "Sequencing depth")
    add_scale_bar(axes[1, 0], spatial, config)
    fig.colorbar(
        total_points, ax=axes[1, 0], shrink=0.8, label="log10(total counts + 1)"
    )

    rhizobia_points = axes[1, 1].scatter(
        spatial[:, 0],
        spatial[:, 1],
        c=adata.obs["rhizobia_fraction"],
        cmap="magma",
        s=14,
        marker="s",
        linewidths=0,
        rasterized=True,
        vmin=0,
        vmax=np.nanquantile(adata.obs["rhizobia_fraction"], 0.99),
    )
    axes[1, 1].invert_yaxis()
    style_spatial_axes(axes[1, 1], "Rhizobial count fraction")
    add_scale_bar(axes[1, 1], spatial, config)
    fig.colorbar(rhizobia_points, ax=axes[1, 1], shrink=0.8, label="Rhizobia / total")

    sns.histplot(
        adata.obs["rhizobia_fraction"],
        bins=35,
        ax=axes[1, 2],
        color="#8b4b8b",
    )
    axes[1, 2].set(
        xlabel="Rhizobia / total counts",
        ylabel="Spatial bins",
        xlim=(0, np.nanquantile(adata.obs["rhizobia_fraction"], 0.99)),
    )

    for label, ax in zip("ABCDEF", axes.flat):
        add_panel_label(ax, label)
    fig.suptitle("Quality control of the P3 nodule spatial transcriptome", fontsize=12)
    save_figure(fig, path)
    plt.close(fig)


def plot_annotation_and_depth(adata: ad.AnnData, path: Path) -> None:
    spatial = np.asarray(adata.obsm["spatial"])
    config = load_config()
    labels = list(adata.obs["organ"].cat.categories)
    palette = {label: TISSUE_COLORS.get(label, "#777777") for label in labels}

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    for label, frame in adata.obs.groupby("organ", observed=True):
        points = spatial[frame.index.map(adata.obs.index.get_loc)]
        axes[0].scatter(
            points[:, 0],
            points[:, 1],
            s=15,
            marker="s",
            linewidths=0,
            label=label,
            color=palette[label],
            rasterized=True,
        )
    axes[0].invert_yaxis()
    style_spatial_axes(axes[0], "Author tissue annotation")
    add_scale_bar(axes[0], spatial, config)
    axes[0].legend(
        title="Tissue zone",
        bbox_to_anchor=(1.01, 1.0),
        loc="upper left",
        frameon=False,
        fontsize=7,
    )

    order = [label for label in TISSUE_COLORS if label in labels]
    sns.boxplot(
        data=adata.obs,
        x="organ",
        y="log10_total_counts_recomputed",
        hue="organ",
        order=order,
        ax=axes[1],
        palette=palette,
        showfliers=False,
        width=0.65,
        linewidth=0.8,
        legend=False,
    )
    sns.stripplot(
        data=adata.obs,
        x="organ",
        y="log10_total_counts_recomputed",
        order=order,
        ax=axes[1],
        color="black",
        size=1.2,
        alpha=0.28,
        jitter=0.2,
        rasterized=True,
    )
    axes[1].tick_params(axis="x", rotation=45)
    axes[1].set(
        xlabel="",
        ylabel="log10(total counts + 1)",
    )
    axes[1].grid(axis="y", color="#d9d9d9", linewidth=0.6)
    axes[1].set_axisbelow(True)
    add_panel_label(axes[0], "A")
    add_panel_label(axes[1], "B")
    fig.suptitle("Tissue annotation and sequencing depth", fontsize=12)
    save_figure(fig, path)
    plt.close(fig)


def plot_normalization_contrast(
    plant: ad.AnnData,
    gene_id: str,
    symbol: str,
    path: Path,
) -> None:
    gene_id = gene_id.replace("_", "")
    if gene_id not in plant.var_names:
        raise ValueError(f"Normalization example gene is absent: {gene_id}")
    config = load_config()
    spatial = np.asarray(plant.obsm["spatial"])
    raw = np.asarray(plant[:, gene_id].layers["counts"].toarray()).ravel()
    normalized = np.asarray(plant[:, gene_id].X.toarray()).ravel()
    raw_upper = np.quantile(raw[raw > 0], 0.99) if np.any(raw > 0) else 1
    norm_upper = (
        np.quantile(normalized[normalized > 0], 0.99) if np.any(normalized > 0) else 1
    )

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    for ax, values, upper, title, cmap in (
        (axes[0], raw, raw_upper, f"{symbol}: raw counts", "viridis"),
        (axes[1], normalized, norm_upper, f"{symbol}: log-normalized", "viridis"),
    ):
        points = ax.scatter(
            spatial[:, 0],
            spatial[:, 1],
            c=values,
            cmap=cmap,
            s=15,
            marker="s",
            linewidths=0,
            vmin=0,
            vmax=upper,
            rasterized=True,
        )
        ax.invert_yaxis()
        style_spatial_axes(ax, title)
        add_scale_bar(ax, spatial, config)
        fig.colorbar(
            points,
            ax=ax,
            shrink=0.75,
            label="Raw counts" if "raw" in title else "log-normalized expression",
        )

    axes[2].scatter(
        np.log10(raw + 1),
        normalized,
        s=12,
        color="#4C78A8",
        alpha=0.7,
        rasterized=True,
    )
    axes[2].set(
        xlabel="log10(raw counts + 1)",
        ylabel="log-normalized expression",
        title=f"{symbol}: per-bin relationship",
    )
    axes[2].grid(color="#e5e5e5", linewidth=0.6)
    axes[2].set_axisbelow(True)
    for label, ax in zip("ABC", axes):
        add_panel_label(ax, label)
    fig.suptitle("Raw counts compared with total-count normalization", fontsize=12)
    save_figure(fig, path)
    plt.close(fig)


def prepare_plant_analysis(
    adata: ad.AnnData, config: dict[str, Any]
) -> tuple[ad.AnnData, np.ndarray]:
    analysis_config = config["analysis"]
    plant_mask = (adata.var["species"] == "Medicago").to_numpy()
    plant = adata[:, plant_mask].copy()
    plant.obs["author_tissue"] = adata.obs["organ"].astype(str)
    plant.obs["author_cluster"] = adata.obs["seurat_clusters"].astype(str)
    plant.layers["counts"] = plant.X.copy()
    sc.pp.filter_cells(plant, min_counts=int(analysis_config["min_plant_counts"]))
    sc.pp.filter_genes(plant, min_cells=int(analysis_config["min_plant_cells"]))
    plant.layers["counts"] = sparse.csr_matrix(plant.layers["counts"])
    sc.pp.normalize_total(plant, target_sum=float(analysis_config["target_sum"]))
    sc.pp.log1p(plant)

    sc.pp.highly_variable_genes(
        plant,
        n_top_genes=int(analysis_config["n_top_genes"]),
        flavor="seurat",
    )
    sc.tl.pca(
        plant,
        n_comps=int(analysis_config["n_pcs"]),
        random_state=int(analysis_config["random_seed"]),
    )
    sc.pp.neighbors(
        plant,
        n_neighbors=int(analysis_config["n_neighbors"]),
        n_pcs=int(analysis_config["n_pcs"]),
        random_state=int(analysis_config["random_seed"]),
    )
    sc.tl.umap(
        plant,
        random_state=int(analysis_config["random_seed"]),
    )
    for key, resolution in analysis_config["resolutions"].items():
        sc.tl.leiden(
            plant,
            resolution=float(resolution),
            key_added=f"leiden_{key}",
            flavor="igraph",
            n_iterations=2,
            directed=False,
            random_state=int(analysis_config["random_seed"]),
        )
    return plant, plant_mask


def rank_region_markers(
    plant: ad.AnnData,
    mapping: pd.DataFrame,
    path: Path,
) -> pd.DataFrame:
    """Rank descriptive region markers after clustering, using author labels only for grouping."""
    sc.tl.rank_genes_groups(
        plant,
        groupby="author_tissue",
        method="wilcoxon",
        n_genes=25,
        use_raw=False,
    )
    symbol_lookup = mapping.drop_duplicates("gene_id").set_index("gene_id")
    frames: list[pd.DataFrame] = []
    for region in plant.obs["author_tissue"].unique():
        frame = sc.get.rank_genes_groups_df(plant, group=str(region)).copy()
        frame.insert(0, "author_region", str(region))
        frame["gene_symbol"] = frame["names"].map(symbol_lookup["acronym"])
        frame["functional_annotation"] = frame["names"].map(
            symbol_lookup["description"]
        )
        frames.append(frame)
    result = pd.concat(frames, ignore_index=True)
    result.to_csv(path, index=False, encoding="utf-8-sig")
    return result


def transfer_analysis_results(adata: ad.AnnData, plant: ad.AnnData) -> ad.AnnData:
    for column in ("leiden_coarse", "leiden_fine", "author_tissue", "author_cluster"):
        adata.obs[column] = pd.Series(index=adata.obs_names, dtype="object")
        adata.obs.loc[plant.obs_names, column] = plant.obs[column].astype(str).values
    for basis_name, basis in (
        ("umap", plant.obsm["X_umap"]),
        ("pca", plant.obsm["X_pca"]),
    ):
        result = np.full((adata.n_obs, basis.shape[1]), np.nan)
        mask = adata.obs_names.isin(plant.obs_names)
        result[mask, :] = basis
        adata.obsm[f"{basis_name}_independent"] = result
    adata.uns["analysis_parameters"] = json.loads(
        json.dumps(plant.uns.get("pca", {}), default=str)
    )
    return adata


def spatial_neighbor_purity(
    labels: pd.Series,
    spatial: np.ndarray,
    n_neighbors: int = 6,
) -> float:
    labels = labels.astype(str).reset_index(drop=True)
    neighbor_model = NearestNeighbors(n_neighbors=min(n_neighbors + 1, len(labels)))
    indices = neighbor_model.fit(spatial).kneighbors(return_distance=False)[:, 1:]
    purity = np.mean(
        [
            np.mean(labels.iloc[row_indices].to_numpy() == labels.iloc[index])
            for index, row_indices in enumerate(indices)
        ]
    )
    return float(purity)


def plot_clustering(
    plant: ad.AnnData,
    config: dict[str, Any],
    figure_path: Path,
    metrics_path: Path,
) -> dict[str, float]:
    spatial = np.asarray(plant.obsm["spatial"])
    primary = config["analysis"]["primary_resolution"]
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    grid = fig.add_gridspec(2, 3)
    axes = [fig.add_subplot(grid[row, col]) for row in range(2) for col in range(3)]

    sc.pl.umap(
        plant,
        color="leiden_coarse",
        ax=axes[0],
        show=False,
        title="Independent UMAP: resolution 0.3",
        frameon=True,
    )
    sc.pl.umap(
        plant,
        color="leiden_fine",
        ax=axes[1],
        show=False,
        title="Independent UMAP: resolution 0.8",
        frameon=True,
    )

    for ax, key, title in (
        (axes[2], "author_tissue", "Author tissue annotation"),
        (axes[3], "leiden_coarse", "Own clusters: resolution 0.3"),
        (axes[4], f"leiden_{primary}", "Own clusters: resolution 0.8"),
    ):
        categories = plant.obs[key].astype("category").cat.categories
        palette = (
            [TISSUE_COLORS.get(str(category), "#777777") for category in categories]
            if key == "author_tissue"
            else sns.color_palette("tab20", n_colors=len(categories))
        )
        for color, category in zip(palette, categories):
            mask = plant.obs[key].astype(str).to_numpy() == str(category)
            ax.scatter(
                spatial[mask, 0],
                spatial[mask, 1],
                color=color,
                label=str(category),
                marker="s",
                linewidths=0,
                rasterized=True,
                s=14,
            )
        ax.invert_yaxis()
        ax.set_aspect("equal")
        style_spatial_axes(ax, title)
        add_scale_bar(ax, spatial, config)

    # Use the spare panel for the author-zone legend; cluster legends are attached to UMAPs.
    handles, labels = axes[2].get_legend_handles_labels()
    axes[5].axis("off")
    axes[5].legend(
        handles,
        labels,
        title="Author tissue zone",
        loc="center",
        frameon=False,
        ncol=1,
        fontsize=8,
    )
    for label, ax in zip("ABCDEF", axes):
        add_panel_label(ax, label)
    fig.suptitle("Independent plant-only clustering", fontsize=12)
    save_figure(fig, figure_path)
    plt.close(fig)

    metrics: dict[str, float] = {}
    for key in ("leiden_coarse", "leiden_fine"):
        metrics[f"{key}_ari_vs_author_tissue"] = float(
            adjusted_rand_score(plant.obs["author_tissue"], plant.obs[key])
        )
        metrics[f"{key}_nmi_vs_author_tissue"] = float(
            normalized_mutual_info_score(plant.obs["author_tissue"], plant.obs[key])
        )
        metrics[f"{key}_n_clusters"] = float(plant.obs[key].nunique())
        metrics[f"{key}_silhouette_pca20"] = float(
            silhouette_score(plant.obsm["X_pca"][:, :20], plant.obs[key])
        )
        metrics[f"{key}_spatial_neighbor_purity"] = spatial_neighbor_purity(
            plant.obs[key],
            spatial,
        )
    metrics["author_tissue_spatial_neighbor_purity"] = spatial_neighbor_purity(
        plant.obs["author_tissue"],
        spatial,
    )
    metrics["plant_bins_analyzed"] = float(plant.n_obs)
    metrics_path.write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return metrics


def plot_cluster_confusion(
    plant: ad.AnnData,
    config: dict[str, Any],
    path: Path,
) -> None:
    primary = config["analysis"]["primary_resolution"]
    tissue_order = [
        label for label in TISSUE_COLORS if label in set(plant.obs["author_tissue"])
    ]
    cross_tab = pd.crosstab(
        pd.Categorical(
            plant.obs["author_tissue"],
            categories=tissue_order,
            ordered=True,
        ),
        plant.obs[f"leiden_{primary}"],
    )
    row_normalized = cross_tab.div(cross_tab.sum(axis=1).replace(0, np.nan), axis=0)

    fig, ax = plt.subplots(figsize=(9, 5.5), constrained_layout=True)
    sns.heatmap(
        row_normalized,
        annot=cross_tab,
        fmt="d",
        cmap="Blues",
        vmin=0,
        vmax=1,
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": "Row-normalized fraction", "shrink": 0.8},
        ax=ax,
    )
    ax.set(
        xlabel=f"Independent Leiden cluster (resolution {config['analysis']['resolutions'][primary]})",
        ylabel="Author tissue annotation",
        title="Cluster composition by author tissue label",
    )
    add_panel_label(ax, "A")
    save_figure(fig, path)
    plt.close(fig)


def plot_marker_dot_matrix(
    means: pd.DataFrame,
    fractions: pd.DataFrame,
    genes: list[str],
    group_order: list[str],
    title: str,
    panel_label: str,
    path: Path,
) -> None:
    means = means.loc[group_order, genes]
    fractions = fractions.loc[group_order, genes]
    scaled = means.copy()
    for gene in genes:
        values = means[gene]
        span = values.max() - values.min()
        scaled[gene] = 0.5 if span == 0 else (values - values.min()) / span

    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    x_coordinates, y_coordinates = np.meshgrid(
        np.arange(len(genes)),
        np.arange(len(group_order)),
    )
    matrix = ax.scatter(
        x_coordinates.ravel(),
        y_coordinates.ravel(),
        c=scaled.to_numpy().ravel(),
        s=8 + 150 * fractions.to_numpy().ravel(),
        cmap="Reds",
        vmin=0,
        vmax=1,
        edgecolors="#555555",
        linewidths=0.25,
    )
    ax.set_xticks(np.arange(len(genes)), genes, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(group_order)), group_order)
    ax.set_xlim(-0.6, len(genes) - 0.4)
    ax.set_ylim(len(group_order) - 0.4, -0.6)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title)
    ax.grid(False)
    ax.tick_params(length=2)
    colorbar = fig.colorbar(matrix, ax=ax, shrink=0.75, pad=0.02)
    colorbar.set_label("Scaled mean expression")
    size_handles = [
        ax.scatter([], [], s=8 + 150 * fraction, c="#999999", edgecolors="#555555")
        for fraction in (0.25, 0.5, 1.0)
    ]
    ax.legend(
        size_handles,
        ["25", "50", "100"],
        title="Detected bins (%)",
        loc="center left",
        bbox_to_anchor=(1.14, 0.5),
        frameon=False,
    )
    add_panel_label(ax, panel_label)
    save_figure(fig, path)
    plt.close(fig)


def plot_marker_genes(
    plant: ad.AnnData,
    marker_table: pd.DataFrame,
    config: dict[str, Any],
    dotplot_path: Path,
    spatial_path: Path,
    table_path: Path,
) -> pd.DataFrame:
    available = marker_table.loc[
        marker_table["available_in_sample"].astype(str).str.lower() == "true"
    ].copy()
    # Family-annotation rows can contain several aliases mapping to the same locus.
    available = available.drop_duplicates("resolved_gene_id", keep="first")
    genes = available["resolved_gene_id"].tolist()
    acronyms = available["requested_acronym"].tolist()
    display = {gene: acronym for gene, acronym in zip(genes, acronyms)}

    # Rename only the selected markers for readable plot labels.
    plot_data = plant[:, genes].copy()
    plot_data.var["symbol"] = [display[g] for g in plot_data.var_names]
    plot_data.var_names = plot_data.var["symbol"].astype(str).tolist()

    primary = config["analysis"]["primary_resolution"]
    expression = pd.DataFrame(
        plot_data.X.toarray() if sparse.issparse(plot_data.X) else plot_data.X,
        index=plot_data.obs_names,
        columns=plot_data.var_names,
    )
    expression["author_tissue"] = plant.obs["author_tissue"].values
    means = expression.groupby("author_tissue", observed=True).mean()
    detections = expression.drop(columns="author_tissue").gt(0)
    detections["author_tissue"] = plant.obs["author_tissue"].values
    fractions = detections.groupby("author_tissue", observed=True).mean()

    tissue_order = [
        label
        for label in TISSUE_COLORS
        if label in set(plant.obs["author_tissue"].astype(str))
    ]
    cluster_order = sorted(
        plant.obs[f"leiden_{primary}"].astype(str).unique(),
        key=lambda value: int(value) if value.isdigit() else value,
    )
    genes_display = plot_data.var_names.tolist()
    plot_marker_dot_matrix(
        means,
        fractions,
        genes_display,
        tissue_order,
        "Published markers by author tissue annotation",
        "A",
        dotplot_path,
    )
    cluster_expression = expression.drop(columns="author_tissue").copy()
    cluster_expression["cluster"] = plant.obs[f"leiden_{primary}"].astype(str).values
    cluster_means = cluster_expression.groupby("cluster", observed=True).mean()
    cluster_detections = expression.drop(columns="author_tissue").gt(0)
    cluster_detections["cluster"] = plant.obs[f"leiden_{primary}"].astype(str).values
    cluster_fractions = cluster_detections.groupby("cluster", observed=True).mean()
    plot_marker_dot_matrix(
        cluster_means,
        cluster_fractions,
        genes_display,
        cluster_order,
        f"Published markers by independent cluster ({primary})",
        "B",
        dotplot_path.with_name(dotplot_path.stem + "_by_cluster" + dotplot_path.suffix),
    )

    spatial = np.asarray(plant.obsm["spatial"])
    selected = plot_data.var_names.tolist()
    ncols = 3
    nrows = int(np.ceil(len(selected) / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(15, 5 * nrows),
        constrained_layout=True,
    )
    axes = np.asarray(axes).reshape(-1)
    for ax, gene in zip(axes, selected):
        values = np.asarray(plot_data[:, gene].X.toarray()).ravel()
        cutoff = np.quantile(values[values > 0], 0.99) if np.any(values > 0) else 1
        points = ax.scatter(
            spatial[:, 0],
            spatial[:, 1],
            c=values,
            cmap="viridis",
            s=15,
            marker="s",
            linewidths=0,
            vmin=0,
            vmax=cutoff,
            rasterized=True,
        )
        ax.invert_yaxis()
        style_spatial_axes(ax, gene)
        add_scale_bar(ax, spatial, config)
        fig.colorbar(points, ax=ax, shrink=0.75, label="log-normalized expression")
    for ax in axes[len(selected) :]:
        ax.axis("off")

    for label, ax in zip("ABCDEFGHI", axes):
        add_panel_label(ax, label)
    fig.suptitle("Spatial expression of published region markers", fontsize=12)
    save_figure(fig, spatial_path)
    plt.close(fig)

    rows: list[dict[str, Any]] = []
    for marker in available.itertuples(index=False):
        gene_name = marker.requested_acronym
        if gene_name not in means.columns:
            continue
        expected_region = {
            "meristematic_zone": "Meristematic zone",
            "proximal_infection_zone": "Proximal infection zone",
            "distal_infection_zone": "Distal infection zone",
            "interzone": "Interzone",
            "nitrogen_fixation_zone": "Nitrogen fixation zone",
            "vascular_bundle": "Vascular bundle",
        }[marker.zone]
        target_mean = float(means.loc[expected_region, gene_name])
        other_mean = float(means.drop(index=expected_region)[gene_name].mean())
        rows.append(
            {
                "zone": marker.zone,
                "acronym": gene_name,
                "gene_id": marker.resolved_gene_id,
                "best_author_region": means[gene_name].idxmax(),
                "mean_expression_best_region": float(means[gene_name].max()),
                "expected_region_mean_expression": target_mean,
                "other_regions_mean_expression": other_mean,
                "expected_region_specificity_log2fc": float(
                    np.log2((target_mean + 1e-6) / (other_mean + 1e-6))
                ),
                "detection_fraction_best_region": float(
                    fractions[gene_name].loc[means[gene_name].idxmax()]
                ),
                "detection_fraction_overall": float(detections[gene_name].mean()),
                "supplementary_evidence_match": marker.supplementary_evidence_match,
            }
        )
    marker_summary = pd.DataFrame(rows)
    marker_summary["expected_region"] = marker_summary["zone"].map(
        {
            "meristematic_zone": "Meristematic zone",
            "proximal_infection_zone": "Proximal infection zone",
            "distal_infection_zone": "Distal infection zone",
            "interzone": "Interzone",
            "nitrogen_fixation_zone": "Nitrogen fixation zone",
            "vascular_bundle": "Vascular bundle",
        }
    )
    marker_summary["matches_expected_region"] = (
        marker_summary["best_author_region"] == marker_summary["expected_region"]
    )
    marker_summary.to_csv(table_path, index=False, encoding="utf-8-sig")
    return marker_summary


def plot_rhizobia_analysis(
    adata: ad.AnnData,
    figure_path: Path,
    table_path: Path,
    statistics_path: Path,
) -> pd.DataFrame:
    spatial = np.asarray(adata.obsm["spatial"])
    config = load_config()
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)

    depth = axes[0, 0].scatter(
        spatial[:, 0],
        spatial[:, 1],
        c=adata.obs["log10_total_counts_recomputed"],
        cmap="viridis",
        s=15,
        marker="s",
        linewidths=0,
        rasterized=True,
    )
    axes[0, 0].invert_yaxis()
    style_spatial_axes(axes[0, 0], "Total counts")
    add_scale_bar(axes[0, 0], spatial, config)
    fig.colorbar(depth, ax=axes[0, 0], shrink=0.75, label="log10(total counts + 1)")

    fraction = axes[0, 1].scatter(
        spatial[:, 0],
        spatial[:, 1],
        c=adata.obs["rhizobia_fraction"],
        cmap="magma",
        s=15,
        marker="s",
        linewidths=0,
        rasterized=True,
        vmin=0,
        vmax=np.nanquantile(adata.obs["rhizobia_fraction"], 0.99),
    )
    axes[0, 1].invert_yaxis()
    style_spatial_axes(axes[0, 1], "Rhizobial count fraction")
    add_scale_bar(axes[0, 1], spatial, config)
    fig.colorbar(fraction, ax=axes[0, 1], shrink=0.75, label="Rhizobia / total")

    order = [
        label for label in TISSUE_COLORS if label in set(adata.obs["organ"].astype(str))
    ]
    sns.violinplot(
        data=adata.obs,
        x="organ",
        y="rhizobia_fraction",
        order=order,
        hue="organ",
        palette={label: TISSUE_COLORS.get(label, "#777777") for label in order},
        cut=0,
        inner="quartile",
        ax=axes[1, 0],
        legend=False,
    )
    axes[1, 0].tick_params(axis="x", rotation=45)
    axes[1, 0].set(
        xlabel="",
        ylabel="Rhizobia / total counts",
        ylim=(0, np.nanquantile(adata.obs["rhizobia_fraction"], 0.99)),
    )

    sns.scatterplot(
        data=adata.obs,
        x="log10_total_counts_recomputed",
        y="rhizobia_fraction",
        hue="organ",
        palette={
            label: TISSUE_COLORS.get(label, "#777777")
            for label in set(adata.obs["organ"].astype(str))
        },
        s=28,
        alpha=0.65,
        ax=axes[1, 1],
        legend=False,
        rasterized=True,
    )
    axes[1, 1].set(
        xlabel="log10(total counts + 1)",
        ylabel="Rhizobia / total counts",
        ylim=(0, np.nanquantile(adata.obs["rhizobia_fraction"], 0.99)),
    )
    for label, ax in zip("ABCD", axes.flat):
        add_panel_label(ax, label)
    fig.suptitle("Host–rhizobium co-detection", fontsize=12)
    save_figure(fig, figure_path)
    plt.close(fig)

    summary = (
        adata.obs.groupby("organ", observed=True)
        .agg(
            n_bins=("rhizobia_fraction", "size"),
            median_rhizobia_fraction=("rhizobia_fraction", "median"),
            q25_rhizobia_fraction=("rhizobia_fraction", lambda x: x.quantile(0.25)),
            q75_rhizobia_fraction=("rhizobia_fraction", lambda x: x.quantile(0.75)),
            median_total_counts=("total_counts_recomputed", "median"),
            median_rhizobia_counts=("rhizobia_counts", "median"),
        )
        .reset_index()
    )
    statistics_rows: list[dict[str, Any]] = []
    if "Nitrogen fixation zone" in set(adata.obs["organ"]):
        nitrogen = adata.obs.loc[
            adata.obs["organ"] == "Nitrogen fixation zone", "rhizobia_fraction"
        ].dropna()
        other = adata.obs.loc[
            adata.obs["organ"] != "Nitrogen fixation zone", "rhizobia_fraction"
        ].dropna()
        mw = mannwhitneyu(nitrogen, other, alternative="two-sided")
        statistics_rows.append(
            {
                "test": "Mann-Whitney U",
                "comparison": "Nitrogen fixation zone vs all other regions",
                "n_group_1": len(nitrogen),
                "n_group_2": len(other),
                "statistic": float(mw.statistic),
                "p_value": float(mw.pvalue),
                "interpretation": (
                    "descriptive; spatial bins are not independent biological replicates"
                ),
            }
        )
    groups = [
        group["rhizobia_fraction"].dropna().to_numpy()
        for _, group in adata.obs.groupby("organ", observed=True)
    ]
    if len(groups) > 2:
        kw = kruskal(*groups)
        statistics_rows.append(
            {
                "test": "Kruskal-Wallis",
                "comparison": "all author tissue regions",
                "n_group_1": len(groups),
                "n_group_2": int(sum(len(group) for group in groups)),
                "statistic": float(kw.statistic),
                "p_value": float(kw.pvalue),
                "interpretation": (
                    "descriptive; spatial bins are not independent biological replicates"
                ),
            }
        )
    pd.DataFrame(statistics_rows).to_csv(
        statistics_path,
        index=False,
        encoding="utf-8-sig",
    )
    summary.to_csv(table_path, index=False, encoding="utf-8-sig")
    return summary


def write_run_metadata(
    adata: ad.AnnData,
    plant: ad.AnnData,
    metrics: dict[str, float],
    config: dict[str, Any],
    sample_reference: pd.DataFrame,
    path: Path,
) -> None:
    mode_counts = {
        "X_equals_counts": bool((adata.X != adata.layers["counts"]).nnz == 0)
    }
    metadata = {
        "generated_at_utc": datetime.now(UTC).isoformat(),
        "git_revision": safe_git_revision(),
        "git_worktree_dirty_at_run": git_worktree_dirty(),
        "paper": {
            "title": (
                "Spatiotemporal co-transcriptomics uncovers coordinated "
                "rhizobia-legume stress programs supporting nitrogen fixation in Medicago"
            ),
            "doi": "10.1016/j.cell.2026.07.033",
            "journal": "Cell",
            "year": 2026,
        },
        "sample": config["data"]["sample_filename"],
        "sample_sha256": config["data"]["sample_sha256"],
        "annotation_sha256": config["data"]["annotation_sha256"],
        "supplementary_sha256": config["data"]["supplementary_sha256"],
        "author_reported_sample_metadata": sample_reference.iloc[0].to_dict(),
        "public_object_bins_minus_table_s13": int(
            adata.n_obs - sample_reference.iloc[0]["Cell Number"]
        ),
        "analysis_parameters": config["analysis"],
        "input_dimensions": {
            "bins": int(adata.n_obs),
            "genes": int(adata.n_vars),
            "plant_genes": int((adata.var["species"] == "Medicago").sum()),
            "rhizobia_genes": int((adata.var["species"] == "Sinorhizobium").sum()),
            "nonzeros": int(adata.X.nnz),
            "total_counts": float(adata.X.sum()),
        },
        "integrity_checks": {
            **mode_counts,
            "unique_gene_ids": bool(adata.var_names.is_unique),
            "finite_coordinates": bool(np.isfinite(adata.obsm["spatial"]).all()),
            "missing_tissue_labels": int(adata.obs["organ"].isna().sum()),
            "plant_analysis_bins": int(plant.n_obs),
            "plant_analysis_genes": int(plant.n_vars),
            "n_count_matches_reconstructed": bool(
                np.allclose(
                    adata.obs["nCount_Spatial"],
                    adata.obs["total_counts_recomputed"],
                )
            ),
            "n_feature_matches_reconstructed": bool(
                np.allclose(
                    adata.obs["nFeature_Spatial"],
                    adata.obs["detected_genes_recomputed"],
                )
            ),
        },
        "software": {
            "python": platform.python_version(),
            "scanpy": version("scanpy"),
            "anndata": version("anndata"),
            "pandas": str(pd.__version__),
            "numpy": str(np.__version__),
            "scipy": version("scipy"),
            "scikit_learn": version("scikit-learn"),
            "matplotlib": version("matplotlib"),
            "seaborn": version("seaborn"),
        },
        "clustering_metrics": metrics,
    }
    path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )


def run_analysis(config_path: Path | None = None) -> dict[str, Any]:
    configure_logging()
    configure_publication_style()
    paths = Paths()
    config = load_config(config_path)
    paths.figure_dir.mkdir(parents=True, exist_ok=True)
    paths.table_dir.mkdir(parents=True, exist_ok=True)
    sample_path, annotation_path = resolve_inputs(config)
    supplementary_path = resolve_supplementary(config)
    sample_reference = load_sample_reference_metadata(
        supplementary_path,
        config["data"]["sample_filename"].removesuffix(".h5ad"),
    )
    sample_reference.to_csv(
        paths.table_dir / "sample_reference.csv",
        index=False,
        encoding="utf-8-sig",
    )

    LOGGER.info("Loading and validating public spatial data")
    adata = load_and_validate(sample_path)
    adata = add_qc_metrics(adata)
    LOGGER.info("Writing QC figures")
    plot_qc_overview(adata, paths.figure_dir / "01_qc_overview.png")
    plot_annotation_and_depth(adata, paths.figure_dir / "02_annotation_and_depth.png")

    LOGGER.info("Preparing plant-only independent clustering")
    plant, _ = prepare_plant_analysis(adata, config)
    adata = transfer_analysis_results(adata, plant)
    normalization_gene = config["analysis"]["normalization_example"]
    plot_normalization_contrast(
        plant,
        normalization_gene["gene_id"],
        normalization_gene["symbol"],
        paths.figure_dir / "01b_normalization_contrast.png",
    )
    LOGGER.info("Writing clustering comparison")
    metrics = plot_clustering(
        plant,
        config,
        paths.figure_dir / "03_independent_clustering.png",
        paths.table_dir / "clustering_metrics.json",
    )
    plot_cluster_confusion(
        plant,
        config,
        paths.figure_dir / "03b_cluster_annotation_confusion.png",
    )
    rank_region_markers(
        plant,
        load_marker_mapping(annotation_path),
        paths.table_dir / "region_marker_discovery.csv",
    )

    LOGGER.info("Resolving official Mt5.0 marker IDs")
    mapping = load_marker_mapping(annotation_path)
    supplementary_evidence = load_supplementary_marker_evidence(supplementary_path)
    marker_table = resolve_markers(
        config["marker_genes"],
        config["marker_preferred_loci"],
        mapping,
        set(plant.var_names),
        supplementary_evidence,
    )
    marker_table.to_csv(
        paths.table_dir / "marker_resolution.csv",
        index=False,
        encoding="utf-8-sig",
    )
    marker_summary = plot_marker_genes(
        plant,
        marker_table,
        config,
        paths.figure_dir / "04_marker_dotplots.png",
        paths.figure_dir / "05_marker_spatial.png",
        paths.table_dir / "marker_specificity.csv",
    )

    LOGGER.info("Analyzing host-rhizobium co-detection")
    rhizobia_summary = plot_rhizobia_analysis(
        adata,
        paths.figure_dir / "06_rhizobia_fraction.png",
        paths.table_dir / "rhizobia_by_region.csv",
        paths.table_dir / "rhizobia_statistics.csv",
    )

    adata.obs.to_csv(
        paths.table_dir / "bin_level_results.csv",
        encoding="utf-8-sig",
    )
    write_run_metadata(
        adata,
        plant,
        metrics,
        config,
        sample_reference,
        paths.table_dir / "run_metadata.json",
    )

    summary = {
        "clustering": metrics,
        "marker_rows": len(marker_summary),
        "rhizobia_region_rows": len(rhizobia_summary),
    }
    LOGGER.info("Analysis complete: %s", summary)
    return summary


if __name__ == "__main__":
    run_analysis()
