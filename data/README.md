# Data provenance

The binary input files are intentionally excluded from Git and are downloaded by
`python scripts/download_data.py`.

| File | Source | Purpose |
|---|---|---|
| `FP200000569BR_C5_bottom4_MXP3.h5ad` | [CNGB STOmics STAMRN](https://db.cngb.org/stomics/stamrn/download/) | Public P3 spatial count matrix, coordinates and author labels |
| `Mt_gene_Acronyms_IDs_latest.xlsx` | [Medicago Mt5.0 genome portal](https://medicago.toulouse.inrae.fr/MtrunA17r5.0-ANR/) | Gene acronym to locus-ID mapping |
| `cell2026_mmc2.xlsx` | [Cell supplementary material](https://ars.els-cdn.com/content/image/1-s2.0-S0092867426008652-mmc2.xlsx) | Tables S1-1 to S1-6, including cell-type marker evidence |

All downloads are checked against the SHA-256 values in
[`config/analysis.yaml`](../config/analysis.yaml). The upstream data and
annotation are not covered by this repository's MIT license.

The Medicago portal currently serves an incomplete TLS certificate chain in some
CI environments. TLS verification is disabled only for that explicit hostname;
the downloaded file must still match the pinned SHA-256 digest.

The sample file is an author-processed object. Reproducing analyses from this
file is not equivalent to repeating FASTQ filtering, genome alignment or
count-matrix construction.

