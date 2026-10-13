# Macroscopic mechanical stress and nuclear envelope rupture in melanoma — analysis notebooks

Code for the single-cell RNA-seq and Stereo-seq spatial transcriptomics analyses in the manuscript.

## Repository layout

```
data/
  regev_lab_cell_cycle_genes.txt   cell-cycle genes used for scoring (scRNA-seq)
  genesets_used.json               mouse gene sets (MSigDB-derived) used for signature scoring
  Tiroshi_EMT.csv                  Tirosh EMT gene modules (Supplementary Fig. S7E, Tirosh EMT1)
  Peritumor_Fibroblast_DEG_*_pseudo.csv  peritumor fibroblast DEG tables (PyDESeq2) used for Fig. 3F
  Melanoma_scRNAseq_data.h5ad      (download, see below)
  stereoseq_annotated/             (download, see below)
scripts/
  make_annotated_h5ad.py           builds the FF <library>_bin30_annotated.h5ad files (provenance)
  make_annotated_h5ad_ffpe.py      builds the FFPE <library>_bin30_annotated.h5ad files (provenance)
notebooks/
  sc_utils.py                      helper functions for scRNA-seq
  spatial_utils.py                 helper functions for Stereo-seq (loading, layer plots, heatmaps)
  01_scRNAseq_celltypes.ipynb
  02_scRNAseq_fibroblasts.ipynb
  03_scRNAseq_macrophages.ipynb
  04_Stereoseq_FF_celltype_distribution.ipynb
  05_Stereoseq_FF_Lrrc15_fibroblasts.ipynb
  06_Stereoseq_FF_melanoma_signatures.ipynb
  07_Stereoseq_FF_exploratory.ipynb
  08_Stereoseq_FFPE_Lmna_KO.ipynb
```

Run the notebooks from the `notebooks/` folder. Figures are written to `notebooks/figures/`.

## Notebooks and figures

| Notebook | Data | Figure panels |
|---|---|---|
| 01_scRNAseq_celltypes | scRNA-seq | Fig. 2B, 2C, 2H; Supplementary Fig. S6A, S6B, S6C, S6F |
| 02_scRNAseq_fibroblasts | scRNA-seq | Fig. 2D, 2E, 2I; Supplementary Fig. S6D, S6G |
| 03_scRNAseq_macrophages | scRNA-seq | Fig. 2F, 2G; Supplementary Fig. S6E |
| 04_Stereoseq_FF_celltype_distribution | Stereo-seq (FF) | Fig. 3B; Supplementary Fig. S7C |
| 05_Stereoseq_FF_Lrrc15_fibroblasts | Stereo-seq (FF) | Fig. 3D, 3E, 3F; Supplementary Fig. S9 |
| 06_Stereoseq_FF_melanoma_signatures | Stereo-seq (FF) | Fig. 3C; Supplementary Fig. S7E, S8A, S10A |
| 07_Stereoseq_FF_exploratory | Stereo-seq (FF) | exploratory analyses (not shown as figure panels) |
| 08_Stereoseq_FFPE_Lmna_KO | Stereo-seq OMNI (FFPE) | Fig. 4E, 4F, 4G; Supplementary Fig. S11B |

- **Fig. 2H / 2I** (notebooks 01, 02): `sc.tl.score_genes` scores (Hallmark EMT, GO:BP cytoplasmic pattern recognition
  receptor signaling) in melanoma cells, and Acta2 / Lrrc15 log-normalized expression in fibroblasts, Mattress vs Running.
  P-values for Acta2 / Lrrc15 are computed in fibroblasts. The manuscript violins were drawn with ggplot2
  (`geom_violin(draw_quantiles = 0.5)`); `sc_utils.plot_group_violin` reproduces them, including ggplot2's median line,
  which is taken from the estimated (trimmed) density rather than from the data.
- **Supplementary Fig. S6F / S6G** (notebooks 01, 02): AUCell scores (`min_n=1`) of melanoma cells / fibroblasts,
  averaged per sample; Running (n = 5) vs Mattress (n = 4), two-sided Student's t-test of the sample means.
- **Fig. 3C / Supplementary Fig. S7E, S10A** (notebook 06): mean AUCell score of melanoma bins per 50 µm layer (-150–700 µm),
  averaged per sample; lines show the group mean and shaded areas mean ± SEM across samples (n = 3 per group).
  Four signatures that are not stored in the h5ad files (GOMF_IMMUNE_RECEPTOR_ACTIVITY, WP_CHEMOKINE_SIGNALING_PATHWAY,
  GOBP_INTERLEUKIN_1_BETA_PRODUCTION and Tirosh EMT1) are computed in the notebook with AUCell (`min_n=1`) on the same
  log-normalized data.
- **Supplementary Fig. S10A** (notebook 06): Hallmark hypoxia AUCell (default `min_n`). The per-bin scores of the
  original analysis could not be reproduced exactly (r ≈ 0.95–0.97 per bin, also when starting from the original
  objects); the layer means and the resulting plot are the same as in the manuscript.
- **Fig. 3E / Supplementary Fig. S8A** (notebooks 05, 06): mean log-normalized expression per 50 µm layer (all bins for
  Acta2 / Lrrc15; melanoma bins for Cxcl2), averaged over the layers of each region (peritumor -150–0 µm, margin
  0–150 µm, center 150–700 µm). Fig. 3E: one-way ANOVA over all region × group combinations with Tukey's post hoc test.
- **Fig. 3F** (notebook 05): pseudobulk of peritumor fibroblast bins, PyDESeq2 for each pair of groups (padj < 0.05),
  genes grouped by the group in which they are highest, and Enrichr (gseapy; needs internet access) on the genes high
  in Running / YAP5SA. The heatmap uses the DEG tables of the original analysis in `data/`; re-running PyDESeq2 with
  another version changes the genes close to padj = 0.05.
- **Fig. 4E / Supplementary Fig. S11B zoom-ins** (notebook 08): the manuscript sections were rotated and cropped by
  hand; `plot_spatial_zoom` uses a rotation angle and window that approximate those panels.
- **Fig. 4F** (notebook 08): mean AUCell score of the melanoma bins of each section, two-tailed t-test. The original
  analysis used all bins (WT 0.0790, KO 0.1005); melanoma bins give WT 0.0804, KO 0.1017.

## Data

**scRNA-seq** — `Melanoma_scRNAseq_data.h5ad` (annotated AnnData; raw data: GEO GSE289587). Place it in `data/`.

**Stereo-seq** — raw FASTQ files are deposited in GEO (accession to be added).
The notebooks use the annotated bin30 objects `<library>_bin30_annotated.h5ad` (one per slide). Place them in
`data/stereoseq_annotated/`.

- FF footpad tumors (notebooks 04–07): FF_Mattress_1–3, FF_Running_1–3, FF_YAP5SA_1–3.
  Representative slides in Fig. 3: FF_Mattress_1, FF_Running_3, FF_YAP5SA_2.
- FFPE flank skin melanomas (notebook 08): FFPE_Skin_WT1 (Lmna wt/wt), FFPE_Skin_KO1 (Lmna flox/flox).

Each file keeps the raw bin30 counts of every bin in `X` (genes in `var_names`) and the bin coordinates in
`obsm['spatial']`. Annotations are in `.obs`; colours for the categorical columns are in `.uns['<column>_colors']`,
and chip ID / bin size / resolution in `.uns['stereo_seq']`.

FF slides:

| Column | Description |
|---|---|
| `sample`, `group` | GEO library name; Mattress / Running / YAP5SA |
| `cell_type`, `cell_type_level2` | RCTD-based cell type annotation (major types / subtypes) |
| `tumor_layer` | 50 µm layer by distance from the tumor border (negative = peritumor). NaN for bins removed by QC |
| `distance_to_border_um` | signed distance from the tumor border (µm) |
| `tumor_border` | bin on the detected tumor border |
| `Lrrc15_border` | Background / Border / FB_Lrrc15+ (Fig. 3D, Supplementary Fig. S9) |
| `aucell_<GENE SET>` | AUCell score of a gene signature, for bins in the layers -150 µm to 700 µm; `aucell_YAP_3` is the YAP signature in Fig. 3C |

FFPE slides:

| Column | Description |
|---|---|
| `sample`, `group` | GEO library name; WT (Lmna wt/wt) / KO (Lmna flox/flox) |
| `cell_type` | RCTD-based cell type (pigmented melanoma merged into Melanoma). NaN for bins removed by QC |
| `cell_type_with_pigment` | same, with pigmented melanoma (`Melanoma_Pigment`) kept separate |
| `aucell_<GENE SET>` | AUCell score of a gene signature (658 signatures); `aucell_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION` is the EMT-like signature in Fig. 4F |

`var['real_gene_name']` keeps the original names of de-duplicated genes (e.g. `Pcdha11-1` → `Pcdha11`).

For the FF slides, the notebooks use QC-passed bins (`tumor_layer` not NaN) for cell type maps and proportions, and
the bins in the layers -150 µm to 700 µm (`SCORED_LAYERS` in `spatial_utils.py`) for signature scores and expression
analyses. For the FFPE slides, they use QC-passed bins (`cell_type` not NaN).

`scripts/make_annotated_h5ad.py` and `scripts/make_annotated_h5ad_ffpe.py` show how these files were assembled from
the bin30 count matrices and the processed analysis objects.

## Environment

Stereo-seq notebooks were run with stereopy 1.2.0, scanpy 1.9.6, anndata 0.9.2, numpy 1.23.5, pandas 1.5.3,
matplotlib 3.7.1, squidpy 1.2.2, decoupler 1.6.0, pydeseq2 0.4.4, gseapy 1.1.3 and statsmodels 0.14.1. The scRNA-seq notebooks were run with scanpy 1.9.6,
matplotlib 3.8.4 and seaborn 0.13.1, and additionally use harmonypy
(`sc.external.pp.harmony_integrate`).
