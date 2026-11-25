from __future__ import annotations

import numpy as np

from analysis.reproduction import (
    Paths,
    add_qc_metrics,
    load_and_validate,
    load_config,
    load_marker_mapping,
    load_sample_reference_metadata,
    load_supplementary_marker_evidence,
    resolve_inputs,
    resolve_markers,
    resolve_supplementary,
)


def test_public_sample_contract() -> None:
    config = load_config()
    sample_path, _ = resolve_inputs(config)
    adata = load_and_validate(sample_path)
    assert adata.shape == (1284, 40160)
    assert adata.var_names.is_unique
    assert set(adata.var["species"]) == {"Medicago", "Sinorhizobium"}
    assert (adata.var["species"] == "Medicago").sum() == 34349
    assert (adata.var["species"] == "Sinorhizobium").sum() == 5811
    assert float(adata.X.sum()) == 3088044.0
    assert np.isfinite(np.asarray(adata.obsm["spatial"])).all()
    assert np.allclose(adata.X.data, np.round(adata.X.data))
    assert np.allclose(
        adata.obsm["spatial"][:, 0],
        adata.obs["array_row"].to_numpy(),
    )


def test_qc_reconstruction_adds_up() -> None:
    config = load_config()
    sample_path, _ = resolve_inputs(config)
    adata = add_qc_metrics(load_and_validate(sample_path))
    reconstructed = adata.obs["plant_counts"] + adata.obs["rhizobia_counts"]
    assert np.allclose(reconstructed, adata.obs["total_counts_recomputed"])
    assert adata.obs["rhizobia_fraction"].between(0, 1).all()


def test_marker_annotation_hashes() -> None:
    config = load_config()
    sample_path, annotation_path = resolve_inputs(config)
    supplementary_path = resolve_supplementary(config)
    assert Paths().config.exists()
    assert sample_path.stat().st_size == 11070767
    assert annotation_path.stat().st_size == 354171
    assert supplementary_path.stat().st_size == 8817432


def test_selected_markers_resolve_uniquely() -> None:
    config = load_config()
    sample_path, annotation_path = resolve_inputs(config)
    supplementary_path = resolve_supplementary(config)
    adata = load_and_validate(sample_path)
    marker_table = resolve_markers(
        config["marker_genes"],
        config["marker_preferred_loci"],
        load_marker_mapping(annotation_path),
        set(adata.var_names),
        load_supplementary_marker_evidence(supplementary_path),
    )
    available = marker_table[
        marker_table["available_in_sample"].astype(str).str.lower() == "true"
    ]
    assert len(available) == 9
    assert available["resolved_gene_id"].is_unique
    assert available["supplementary_evidence_match"].all()


def test_normalization_example_is_present() -> None:
    config = load_config()
    sample_path, _ = resolve_inputs(config)
    adata = load_and_validate(sample_path)
    gene_id = config["analysis"]["normalization_example"]["gene_id"].replace("_", "")
    assert gene_id in adata.var_names
    assert adata[:, gene_id].X.nnz > 0.8 * adata.n_obs


def test_author_sample_metadata_is_pinned() -> None:
    config = load_config()
    supplementary_path = resolve_supplementary(config)
    reference = load_sample_reference_metadata(
        supplementary_path,
        config["data"]["sample_filename"].removesuffix(".h5ad"),
    )
    assert int(reference.iloc[0]["Cell Number"]) == 1270
    assert np.isclose(float(reference.iloc[0]["Median UMI"]), 1739)
