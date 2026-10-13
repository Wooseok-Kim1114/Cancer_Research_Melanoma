"""Add the analysis annotations to the bin30 count matrices of the FFPE Stereo-seq (OMNI) slides.

Input  : <h5ad_dir>/FFPE_Skin_{WT1|KO1}_bin30_loaded.h5ad      (raw bin30 counts, all bins)
         <processed_dir>/250214_Skin_{WT2|K04}_scored.h5ad      (QC-passed bins, RCTD cell types, AUCell scores)
Output : <out_dir>/FFPE_Skin_{WT1|KO1}_bin30_annotated.h5ad
         X = raw counts for every bin (unchanged), obsm['spatial'] = bin coordinates,
         var['real_gene_name'] = gene name before de-duplication (e.g. Pcdha11-1 -> Pcdha11).
         .obs: sample, group, cell_type, cell_type_with_pigment, aucell_<gene set>
         Bins that did not pass QC have NaN annotations (cell_type is NaN).
"""
import os
import sys

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

h5ad_dir, processed_dir, out_dir = sys.argv[1:4]
os.makedirs(out_dir, exist_ok=True)

# processed sample ID -> (GEO library name, group)
SAMPLES = {'Skin_WT2': ('FFPE_Skin_WT1', 'WT'), 'Skin_K04': ('FFPE_Skin_KO1', 'KO')}
# categorical annotations: original column -> new column
CATEGORICAL = {'refine2': 'cell_type', 'refine': 'cell_type_with_pigment'}


def read_obs(path):
    """Read .obs and .uns colours only (the processed files are 16-23 GB)."""
    with h5py.File(path, 'r') as h:
        obs = ad.experimental.read_elem(h['obs']) if hasattr(ad, 'experimental') else None
        colors = {k: [c.decode() if isinstance(c, bytes) else str(c) for c in h['uns'][k][:]]
                  for k in h['uns'].keys() if k.endswith('_colors')}
        sn = h['uns']['sn']
        chip = ad.experimental.read_elem(sn)
        return obs, colors, chip


def aucell_name(col):
    """'<GENE SET>_Signature' -> 'aucell_<GENE SET>' ('YAP_Signature' -> 'aucell_YAP')."""
    return 'aucell_' + col[:-len('_Signature')]


for pid, (geo, group) in SAMPLES.items():
    raw = sc.read(f'{h5ad_dir}/{geo}_bin30_loaded.h5ad')
    pobs, colors, chip = read_obs(f'{processed_dir}/250214_{pid}_scored.h5ad')

    # consistency checks against the processed object
    assert pobs.index.isin(raw.obs_names).all()
    assert np.allclose(raw[pobs.index].obsm['spatial'], pobs[['x', 'y']].values)

    obs = pd.DataFrame(index=raw.obs_names)
    obs['sample'] = pd.Categorical([geo] * raw.n_obs)
    obs['group'] = pd.Categorical([group] * raw.n_obs)
    uns = {}
    for old, new in CATEGORICAL.items():
        cats = pobs[old].cat.categories
        obs[new] = pd.Categorical(pobs[old].reindex(raw.obs_names), categories=cats)
        uns[f'{new}_colors'] = colors[f'{old}_colors']

    # 'EMT-like Signature' (Fig. 4F) is identical to HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Signature
    assert np.array_equal(pobs['EMT-like Signature'].values,
                          pobs['HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_Signature'].values)
    sig_cols = [c for c in pobs.columns if '_Signature' in c]
    scores = pd.DataFrame({aucell_name(c): pobs[c].reindex(raw.obs_names) for c in sig_cols})
    obs = pd.concat([obs, scores], axis=1)
    assert (obs['cell_type'].notna() == raw.obs_names.isin(pobs.index)).all()

    uns['stereo_seq'] = {'chip_id': str(chip['sn'].iloc[0]) if hasattr(chip, 'iloc') else str(chip),
                         'bin_size': int(raw.uns['bin_size']), 'resolution_nm': int(raw.uns['resolution'])}
    uns['annotation_info'] = {
        'group': 'Tyr-CreERT2;BrafCA;Ptenflox/flox with Lmna wt/wt (WT) or Lmna flox/flox (KO)',
        'cell_type': 'RCTD-based cell type (Melanoma includes pigmented melanoma); NaN = bin removed by QC',
        'cell_type_with_pigment': 'same, with pigmented melanoma (Melanoma_Pigment) kept separate',
        'aucell': 'AUCell scores of gene signatures; aucell_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = EMT-like signature (Fig. 4F); '
                  'aucell_YAP = REACTOME YAP1/WWTR1(TAZ)-stimulated genes',
    }
    out_adata = ad.AnnData(X=raw.X, obs=obs, var=raw.var[['real_gene_name']],  # original names of de-duplicated genes
                           obsm={'spatial': raw.obsm['spatial']}, uns=uns)
    out = f'{out_dir}/{geo}_bin30_annotated.h5ad'
    out_adata.write(out, compression='gzip')
    print(f'{pid:9s} -> {os.path.basename(out)}  bins={raw.n_obs}  annotated={obs.cell_type.notna().sum()}  '
          f'aucell={scores.shape[1]}', flush=True)
