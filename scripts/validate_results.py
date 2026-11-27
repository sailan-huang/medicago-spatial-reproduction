from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"

REQUIRED_FIGURE_STEMS = [
    "01_qc_overview",
    "01b_normalization_contrast",
    "02_annotation_and_depth",
    "03_independent_clustering",
    "03b_cluster_annotation_confusion",
    "04_marker_dotplots",
    "04_marker_dotplots_by_cluster",
    "05_marker_spatial",
    "06_rhizobia_fraction",
]


def validate() -> None:
    for stem in REQUIRED_FIGURE_STEMS:
        path = FIGURES / f"{stem}.pdf"
        if not path.exists() or path.stat().st_size == 0:
            raise AssertionError(f"Missing or empty PDF figure: {path}")
        if len(PdfReader(path).pages) != 1:
            raise AssertionError(f"Figure PDF must contain one page: {path}")

    unexpected = [
        path
        for path in FIGURES.iterdir()
        if path.is_file() and path.suffix.lower() != ".pdf"
    ]
    if unexpected:
        raise AssertionError(
            "Non-PDF files remain in results/figures: "
            + ", ".join(str(path.name) for path in unexpected)
        )

    metadata = json.loads((TABLES / "run_metadata.json").read_text(encoding="utf-8"))
    assert metadata["input_dimensions"]["bins"] == 1284
    assert metadata["input_dimensions"]["genes"] == 40160
    assert metadata["integrity_checks"]["X_equals_counts"] is True
    assert metadata["integrity_checks"]["n_count_matches_reconstructed"] is True
    assert metadata["integrity_checks"]["n_feature_matches_reconstructed"] is True

    bins = pd.read_csv(TABLES / "bin_level_results.csv")
    assert len(bins) == 1284
    assert bins["rhizobia_fraction"].between(0, 1).all()
    assert (
        bins["plant_counts"] + bins["rhizobia_counts"]
        == bins["total_counts_recomputed"]
    ).all()
    sample_reference = pd.read_csv(TABLES / "sample_reference.csv")
    assert int(sample_reference.iloc[0]["Cell Number"]) == 1270
    assert metadata["public_object_bins_minus_table_s13"] == 14

    metrics = json.loads(
        (TABLES / "clustering_metrics.json").read_text(encoding="utf-8")
    )
    for key in (
        "leiden_coarse_silhouette_pca20",
        "leiden_fine_silhouette_pca20",
        "leiden_coarse_spatial_neighbor_purity",
        "leiden_fine_spatial_neighbor_purity",
    ):
        assert -1 <= metrics[key] <= 1

    markers = pd.read_csv(TABLES / "marker_specificity.csv")
    assert len(markers) == 9
    assert markers["detection_fraction_overall"].between(0, 1).all()
    assert markers["supplementary_evidence_match"].all()

    stats = pd.read_csv(TABLES / "rhizobia_statistics.csv")
    assert len(stats) == 2
    assert stats["p_value"].between(0, 1).all()

    print(
        "Validated:",
        {
            "figures": len(REQUIRED_FIGURE_STEMS),
            "figure_format": "PDF only",
            "bins": len(bins),
            "markers": len(markers),
            "statistical_tests": len(stats),
        },
    )


if __name__ == "__main__":
    try:
        validate()
    except Exception as error:
        print(f"Result validation failed: {error}", file=sys.stderr)
        raise
