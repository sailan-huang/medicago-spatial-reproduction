from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass
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
from scipy.stats import mannwhitneyu
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

LOGGER = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[1]


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
    path = (
        Paths().root
        / "data"
        / "external"
        / config["data"]["supplementary_filename"]
    )
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
        raise ValueError("The input must contain spatial coordinates in obsm['spatial'].")
    if not adata.var_names.is_unique:
        raise ValueError("Gene identifiers are not unique.")

    # The public file stores counts in both .X and layers['counts'].
    adata.X = adata.layers["counts"].copy()
    adata.var["species"] = [classify_gene(g) for g in adata.var_names]
    if (adata.var["species"] == "unknown").any():
        examples = adata.var_names[adata.var["species"] == "unknown"][:10].tolist()
        raise ValueError(f"Unclassified gene identifiers found: {examples}")
    if not np.isfinite(np.asarray(adata.obsm["spatial"])).all():
        raise ValueError("Spatial coordinates contain non-finite values.")
    if adata.obs["organ"].isna().any():
        raise ValueError("Tissue annotations contain missing values.")
    if not isinstance(adata.X, sparse.spmatrix):
        adata.X = sparse.csr_matrix(adata.X)
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
        mapping["mt5_locus"]
        .astype(str)
        .str.replace("_", "", regex=False)
        .str.strip()
    )
    mapping["acronym"] = mapping["acronym"].astype(str)
    return mapping


def resolve_markers(
    marker_config: dict[str, list[str]],
    preferred_loci: dict[str, str],
    mapping: pd.DataFrame,
    available_genes: set[str],
) -> pd.DataFrame:
    rows: list[dict[str, str]] = []
    for zone, acronyms in marker_config.items():
        for acronym in acronyms:
            preferred_locus = str(preferred_loci.get(acronym, "")).replace("_", "")
            preferred_rows = mapping[mapping["gene_id"] == preferred_locus].copy()
            if not preferred_rows.empty:
                best = preferred_rows.iloc[0]
                rows.append(
                    {
                        "zone": zone,
                        "requested_acronym": acronym,
                        "resolved_gene_id": str(best["gene_id"]),
                        "matched_acronym": str(best["acronym"]),
                        "description": str(best["description"]),
                        "available_in_sample": str(
                            str(best["gene_id"]) in available_genes
                        ),
                    }
                )
                continue
            candidates = mapping[
                mapping["acronym"]
                .str.split("|", regex=False)
                .apply(lambda values: acronym in values)
            ].copy()
            if candidates.empty:
                rows.append(
                    {
                        "zone": zone,
                        "requested_acronym": acronym,
                        "resolved_gene_id": "",
                        "matched_acronym": "",
                        "description": "",
                        "available_in_sample": "False",
                    }
                )
                continue
            # Resolve each alias to the first Mt5.0 locus present in this sample.
            candidates["available"] = candidates["gene_id"].isin(available_genes)
            candidates["exact_alias"] = candidates["acronym"].apply(
                lambda value: str(value).split("|").count(acronym) == 1
                and len(str(value).split("|")) == 1
            )
            candidates = candidates.sort_values(
                ["available", "exact_alias"],
                ascending=[False, False],
            )
            best = candidates.iloc[0]
            rows.append(
                {
                    "zone": zone,
                    "requested_acronym": acronym,
                    "resolved_gene_id": str(best["gene_id"]),
                    "matched_acronym": str(best["acronym"]),
                    "description": str(best["description"]),
                    "available_in_sample": str(bool(best["available"])),
                }
            )
    return pd.DataFrame(rows)


def plot_qc_overview(adata: ad.AnnData, path: Path) -> None:
    spatial = np.asarray(adata.obsm["spatial"])
    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)

    sns.histplot(
        np.log10(adata.obs["total_counts_recomputed"] + 1),
        bins=35,
        ax=axes[0, 0],
        color="#35618f",
    )
    axes[0, 0].set(
        title="Total counts per spatial bin",
        xlabel="log10(total counts + 1)",
        ylabel="Bins",
    )

    sns.histplot(
        adata.obs["detected_genes_recomputed"],
        bins=35,
        ax=axes[0, 1],
        color="#c87941",
    )
    axes[0, 1].set(
        title="Detected genes per spatial bin",
        xlabel="Detected genes",
        ylabel="Bins",
    )

    total_points = axes[1, 0].scatter(
        spatial[:, 0],
        spatial[:, 1],
        c=adata.obs["log10_total_counts_recomputed"],
        cmap="viridis",
        s=18,
    )
    axes[1, 0].invert_yaxis()
    axes[1, 0].set_aspect("equal")
    axes[1, 0].set(title="Spatial distribution of sequencing depth", xlabel="x", ylabel="y")
    fig.colorbar(total_points, ax=axes[1, 0], shrink=0.8, label="log10(total counts + 1)")

    rhizobia_points = axes[1, 1].scatter(
        spatial[:, 0],
        spatial[:, 1],
        c=adata.obs["rhizobia_fraction"],
        cmap="magma",
        s=18,
        vmin=0,
        vmax=np.nanquantile(adata.obs["rhizobia_fraction"], 0.99),
    )
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_aspect("equal")
    axes[1, 1].set(title="Rhizobial count fraction", xlabel="x", ylabel="y")
    fig.colorbar(rhizobia_points, ax=axes[1, 1], shrink=0.8, label="Rhizobia / total")

    fig.suptitle("P3 nodule sample: quality-control overview", fontweight="bold")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_annotation_and_depth(adata: ad.AnnData, path: Path) -> None:
    spatial = np.asarray(adata.obsm["spatial"])
    labels = list(adata.obs["organ"].cat.categories)
    palette = dict(zip(labels, sns.color_palette("tab10", n_colors=len(labels))))

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)
    for label, frame in adata.obs.groupby("organ", observed=True):
        points = spatial[frame.index.map(adata.obs.index.get_loc)]
        axes[0].scatter(points[:, 0], points[:, 1], s=22, label=label, color=palette[label])
    axes[0].invert_yaxis()
    axes[0].set_aspect("equal")
    axes[0].set(title="Author tissue annotation", xlabel="x", ylabel="y")
    axes[0].legend(
        title="Tissue zone",
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        frameon=False,
        fontsize=8,
    )

    sns.boxplot(
        data=adata.obs,
        x="organ",
        y="log10_total_counts_recomputed",
        hue="organ",
        ax=axes[1],
        palette=palette,
        showfliers=False,
        legend=False,
    )
    sns.stripplot(
        data=adata.obs,
        x="organ",
        y="log10_total_counts_recomputed",
        ax=axes[1],
        color="black",
        size=1.5,
        alpha=0.35,
    )
    axes[1].tick_params(axis="x", rotation=55)
    axes[1].set(
        title="Sequencing depth by author annotation",
        xlabel="",
        ylabel="log10(total counts + 1)",
    )
    fig.suptitle("Tissue labels and depth are inspected together", fontweight="bold")
    fig.savefig(path, dpi=220)
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
    for basis_name, basis in (("umap", plant.obsm["X_umap"]), ("pca", plant.obsm["X_pca"])):
        result = np.full((adata.n_obs, basis.shape[1]), np.nan)
        mask = adata.obs_names.isin(plant.obs_names)
        result[mask, :] = basis
        adata.obsm[f"{basis_name}_independent"] = result
    adata.uns["analysis_parameters"] = json.loads(json.dumps(plant.uns.get("pca", {}), default=str))
    return adata


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
        palette = sns.color_palette("tab20", n_colors=len(categories))
        for color, category in zip(palette, categories):
            mask = plant.obs[key].astype(str).to_numpy() == str(category)
            ax.scatter(
                spatial[mask, 0],
                spatial[mask, 1],
                s=20,
                color=color,
                label=str(category),
            )
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set(title=title, xlabel="x", ylabel="y")

    # Compact legend for the higher-cardinality clustering result.
    handles, labels = axes[4].get_legend_handles_labels()
    axes[5].axis("off")
    axes[5].legend(
        handles,
        labels,
        title="Resolution 0.8 cluster",
        loc="center",
        frameon=False,
        ncol=2,
        fontsize=9,
    )
    fig.suptitle(
        "Independent plant-only clustering compared with author annotation",
        fontweight="bold",
    )
    fig.savefig(figure_path, dpi=220)
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
    metrics["plant_bins_analyzed"] = float(plant.n_obs)
    metrics_path.write_text(
        json.dumps(metrics, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return metrics


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
    dot_by_tissue = sc.pl.dotplot(
        plot_data,
        var_names=plot_data.var_names.tolist(),
        groupby="author_tissue",
        standard_scale="var",
        show=False,
        title="Markers by author tissue annotation",
        return_fig=True,
    )
    dot_by_cluster = sc.pl.dotplot(
        plot_data,
        var_names=plot_data.var_names.tolist(),
        groupby=f"leiden_{primary}",
        standard_scale="var",
        show=False,
        title=f"Markers by independent cluster ({primary})",
        return_fig=True,
    )
    dot_by_tissue.make_figure()
    dot_by_tissue.fig.set_size_inches(9, 7)
    dot_by_tissue.fig.savefig(dotplot_path, dpi=240, bbox_inches="tight")
    dot_by_cluster.make_figure()
    dot_by_cluster.fig.set_size_inches(9, 7)
    dot_by_cluster.fig.savefig(
        dotplot_path.with_name(dotplot_path.stem + "_by_cluster" + dotplot_path.suffix),
        dpi=240,
        bbox_inches="tight",
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
            s=22,
            vmin=0,
            vmax=cutoff,
        )
        ax.invert_yaxis()
        ax.set_aspect("equal")
        ax.set(title=f"Normalized expression: {gene}", xlabel="x", ylabel="y")
        fig.colorbar(points, ax=ax, shrink=0.75, label="log-normalized expression")
    for ax in axes[len(selected) :]:
        ax.axis("off")

    fig.suptitle(
        "Official Medicago Mt5.0 marker mapping and spatial expression",
        fontweight="bold",
    )
    fig.savefig(spatial_path, dpi=220)
    plt.close(fig)

    # Descriptive specificity: mean expression per author region and its best region.
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

    rows: list[dict[str, Any]] = []
    for marker in available.itertuples(index=False):
        gene_name = marker.requested_acronym
        if gene_name not in means.columns:
            continue
        rows.append(
            {
                "zone": marker.zone,
                "acronym": gene_name,
                "gene_id": marker.resolved_gene_id,
                "best_author_region": means[gene_name].idxmax(),
                "mean_expression_best_region": float(means[gene_name].max()),
                "detection_fraction_best_region": float(
                    fractions[gene_name].loc[means[gene_name].idxmax()]
                ),
                "detection_fraction_overall": float(
                    fractions[gene_name].mean()
                ),
            }
        )
    marker_summary = pd.DataFrame(rows)
    expected_regions = {
        "meristematic_zone": "Meristematic zone",
        "proximal_infection_zone": "Proximal infection zone",
        "distal_infection_zone": "Distal infection zone",
        "interzone": "Interzone",
        "nitrogen_fixation_zone": "Nitrogen fixation zone",
        "vascular_bundle": "Vascular bundle",
    }
    marker_summary["expected_region"] = marker_summary["zone"].map(expected_regions)
    marker_summary["matches_expected_region"] = (
        marker_summary["best_author_region"]
        == marker_summary["expected_region"]
    )
    marker_summary.to_csv(table_path, index=False, encoding="utf-8-sig")
    return marker_summary


def plot_rhizobia_analysis(
    adata: ad.AnnData,
    figure_path: Path,
    table_path: Path,
) -> pd.DataFrame:
    spatial = np.asarray(adata.obsm["spatial"])
    sns.set_theme(style="whitegrid", context="notebook")
    fig, axes = plt.subplots(2, 2, figsize=(14, 11), constrained_layout=True)

    depth = axes[0, 0].scatter(
        spatial[:, 0],
        spatial[:, 1],
        c=adata.obs["log10_total_counts_recomputed"],
        cmap="viridis",
        s=22,
    )
    axes[0, 0].invert_yaxis()
    axes[0, 0].set_aspect("equal")
    axes[0, 0].set(title="Total counts (quality context)", xlabel="x", ylabel="y")
    fig.colorbar(depth, ax=axes[0, 0], shrink=0.75, label="log10(total counts + 1)")

    fraction = axes[0, 1].scatter(
        spatial[:, 0],
        spatial[:, 1],
        c=adata.obs["rhizobia_fraction"],
        cmap="magma",
        s=22,
        vmin=0,
        vmax=np.nanquantile(adata.obs["rhizobia_fraction"], 0.99),
    )
    axes[0, 1].invert_yaxis()
    axes[0, 1].set_aspect("equal")
    axes[0, 1].set(title="Rhizobial count fraction", xlabel="x", ylabel="y")
    fig.colorbar(fraction, ax=axes[0, 1], shrink=0.75, label="Rhizobia / total")

    order = (
        adata.obs.groupby("organ", observed=True)["rhizobia_fraction"]
        .median()
        .sort_values(ascending=False)
        .index.tolist()
    )
    sns.violinplot(
        data=adata.obs,
        x="organ",
        y="rhizobia_fraction",
        order=order,
        hue="organ",
        palette="magma",
        cut=0,
        inner="quartile",
        ax=axes[1, 0],
        legend=False,
    )
    axes[1, 0].tick_params(axis="x", rotation=55)
    axes[1, 0].set(
        title="Rhizobial fraction by author annotation",
        xlabel="",
        ylabel="Rhizobia / total counts",
    )

    sns.scatterplot(
        data=adata.obs,
        x="log10_total_counts_recomputed",
        y="rhizobia_fraction",
        hue="organ",
        palette="tab10",
        s=28,
        alpha=0.75,
        ax=axes[1, 1],
        legend=False,
    )
    axes[1, 1].set(
        title="Fraction depends on bin depth and RNA composition",
        xlabel="log10(total counts + 1)",
        ylabel="Rhizobia / total counts",
    )
    fig.suptitle("Host–rhizobium co-detection: descriptive spatial analysis", fontweight="bold")
    fig.savefig(figure_path, dpi=220)
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
    if "Nitrogen fixation zone" in set(adata.obs["organ"]):
        nitrogen = adata.obs.loc[
            adata.obs["organ"] == "Nitrogen fixation zone", "rhizobia_fraction"
        ].dropna()
        other = adata.obs.loc[
            adata.obs["organ"] != "Nitrogen fixation zone", "rhizobia_fraction"
        ].dropna()
        p_value = mannwhitneyu(nitrogen, other, alternative="two-sided").pvalue
        summary["mannwhitney_p_nitrogen_vs_rest"] = p_value
    summary.to_csv(table_path, index=False, encoding="utf-8-sig")
    return summary


def write_run_metadata(
    adata: ad.AnnData,
    plant: ad.AnnData,
    metrics: dict[str, float],
    config: dict[str, Any],
    path: Path,
) -> None:
    mode_counts = {
        "X_equals_counts": bool(
            (adata.X != adata.layers["counts"]).nnz == 0
        )
    }
    metadata = {
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
        "analysis_parameters": config["analysis"],
        "input_dimensions": {
            "bins": int(adata.n_obs),
            "genes": int(adata.n_vars),
            "plant_genes": int((adata.var["species"] == "Medicago").sum()),
            "rhizobia_genes": int(
                (adata.var["species"] == "Sinorhizobium").sum()
            ),
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
        },
        "software": {
            "python": "3.12",
            "scanpy": version("scanpy"),
            "anndata": version("anndata"),
            "pandas": str(pd.__version__),
            "numpy": str(np.__version__),
        },
        "clustering_metrics": metrics,
    }
    path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")


def run_analysis(config_path: Path | None = None) -> dict[str, Any]:
    configure_logging()
    paths = Paths()
    config = load_config(config_path)
    paths.figure_dir.mkdir(parents=True, exist_ok=True)
    paths.table_dir.mkdir(parents=True, exist_ok=True)
    sample_path, annotation_path = resolve_inputs(config)
    resolve_supplementary(config)

    LOGGER.info("Loading and validating public spatial data")
    adata = load_and_validate(sample_path)
    adata = add_qc_metrics(adata)
    LOGGER.info("Writing QC figures")
    plot_qc_overview(adata, paths.figure_dir / "01_qc_overview.png")
    plot_annotation_and_depth(adata, paths.figure_dir / "02_annotation_and_depth.png")

    LOGGER.info("Preparing plant-only independent clustering")
    plant, _ = prepare_plant_analysis(adata, config)
    adata = transfer_analysis_results(adata, plant)
    LOGGER.info("Writing clustering comparison")
    metrics = plot_clustering(
        plant,
        config,
        paths.figure_dir / "03_independent_clustering.png",
        paths.table_dir / "clustering_metrics.json",
    )
    rank_region_markers(
        plant,
        load_marker_mapping(annotation_path),
        paths.table_dir / "region_marker_discovery.csv",
    )

    LOGGER.info("Resolving official Mt5.0 marker IDs")
    mapping = load_marker_mapping(annotation_path)
    marker_table = resolve_markers(
        config["marker_genes"],
        config["marker_preferred_loci"],
        mapping,
        set(plant.var_names),
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
        paths.table_dir / "run_metadata.json",
    )

    summary = {
        "clustering": metrics,
        "marker_rows": int(len(marker_summary)),
        "rhizobia_region_rows": int(len(rhizobia_summary)),
    }
    LOGGER.info("Analysis complete: %s", summary)
    return summary


if __name__ == "__main__":
    run_analysis()

