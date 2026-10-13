# Mechanical Stress Induces Nuclear Envelope Rupture in Melanoma Cells to Reorganize Tumor–Stromal Interactions — analysis notebooks

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
