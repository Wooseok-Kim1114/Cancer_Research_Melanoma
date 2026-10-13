"""Add the analysis annotations to the bin30 count matrices of the FF Stereo-seq slides.

Input  : <h5ad_dir>/FF_{Sponge|Tread|YAP5SA}_{1..3}_bin30_loaded.h5ad   (raw bin30 counts, all bins)
         <processed_dir>/241103_<id>_margin_lognorm.h5ad                (cell types, tumor layers)
         <processed_dir>/241103_<id>_scored.h5ad                        (AUCell scores)
Output : <out_dir>/<GEO library name>_bin30_annotated.h5ad
         X = raw counts for every bin (unchanged), obsm['spatial'] = bin coordinates.
         .obs: sample, group, cell_type, cell_type_level2, tumor_layer, distance_to_border_um,
               tumor_border, Lrrc15_border, aucell_<gene set>
         Bins that did not pass QC have NaN annotations (tumor_layer is NaN).
         AUCell scores exist for bins in the layers -150 um to 700 um (SCORED_LAYERS).
"""
import os
import sys

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

h5ad_dir, processed_dir, out_dir = sys.argv[1:4]
os.makedirs(out_dir, exist_ok=True)

# processed sample ID -> (input file stem, GEO library name, group)
SAMPLES = {
    'Sponge_D1': ('FF_Sponge_1', 'FF_Mattress_1', 'Mattress'),
    'Sponge_E4': ('FF_Sponge_2', 'FF_Mattress_2', 'Mattress'),
    'Sponge11': ('FF_Sponge_3', 'FF_Mattress_3', 'Mattress'),
    'Tread_E2': ('FF_Tread_1', 'FF_Running_1', 'Running'),
    'Tread_G6': ('FF_Tread_2', 'FF_Running_2', 'Running'),
    'Tread8': ('FF_Tread_3', 'FF_Running_3', 'Running'),
    'YAP5SA_C2': ('FF_YAP5SA_1', 'FF_YAP5SA_1', 'YAP5SA'),
    'YAP5SA_E6': ('FF_YAP5SA_2', 'FF_YAP5SA_2', 'YAP5SA'),
    'YAP5SA12': ('FF_YAP5SA_3', 'FF_YAP5SA_3', 'YAP5SA'),
}
# categorical annotations: original column -> new column (colours copied from .uns)
CATEGORICAL = {
    'refine': 'cell_type',
    'refine_level2': 'cell_type_level2',
    'layer': 'tumor_layer',
    'Lrrc15_border': 'Lrrc15_border',
}
DNB_TO_UM = 0.5  # 'distance' is stored in DNB units (500 nm)
SCORED_LAYERS = ['-150um~-100um', '-100um~-50um', '-50um~0um'] + [f'{i}um~{i + 50}um' for i in range(0, 700, 50)]


def aucell_name(col):
    """'<GENESET>_Signature' -> 'aucell_<GENESET>', 'YAP_Signature_3' -> 'aucell_YAP_3'."""
    if col.startswith('YAP_Signature_'):
        return 'aucell_YAP_' + col.rsplit('_', 1)[1]
    return 'aucell_' + col[:-len('_Signature')]


for pid, (stem, geo, group) in SAMPLES.items():
    raw = sc.read(f'{h5ad_dir}/{stem}_bin30_loaded.h5ad')
    margin = sc.read(f'{processed_dir}/241103_{pid}_margin_lognorm.h5ad')
    scored = sc.read(f'{processed_dir}/241103_{pid}_scored.h5ad')

    # consistency checks against the processed objects
    assert margin.obs_names.isin(raw.obs_names).all() and scored.obs_names.isin(margin.obs_names).all()
    assert (raw.var_names == margin.var_names).all()
    sub = raw[margin.obs_names]
    assert abs(sub.X - margin.layers['counts']).max() == 0
    assert np.allclose(sub.obsm['spatial'], margin.obsm['spatial'])

    obs = pd.DataFrame(index=raw.obs_names)
    obs['sample'] = pd.Categorical([geo] * raw.n_obs)
    obs['group'] = pd.Categorical([group] * raw.n_obs)
    uns = {}
    for old, new in CATEGORICAL.items():
        cats = margin.obs[old].cat.categories
        obs[new] = pd.Categorical(margin.obs[old].reindex(raw.obs_names), categories=cats)
        uns[f'{new}_colors'] = list(margin.uns[f'{old}_colors'])
    obs['distance_to_border_um'] = margin.obs['distance'].reindex(raw.obs_names) * DNB_TO_UM
    obs['tumor_border'] = margin.obs['border'].reindex(raw.obs_names).astype('boolean')

    sig_cols = [c for c in scored.obs.columns if '_Signature' in c]
    scores = pd.DataFrame({aucell_name(c): scored.obs[c].reindex(raw.obs_names) for c in sig_cols})
    obs = pd.concat([obs, scores], axis=1)

    # the processed subsets are fully described by tumor_layer
    assert (obs['tumor_layer'].notna() == raw.obs_names.isin(margin.obs_names)).all()
    assert (obs['tumor_layer'].isin(SCORED_LAYERS) == raw.obs_names.isin(scored.obs_names)).all()

    sn = raw.uns['sn']
    uns['stereo_seq'] = {'chip_id': str(sn['sn'].iloc[0]) if hasattr(sn, 'iloc') else str(sn),
                         'bin_size': int(raw.uns['bin_size']), 'resolution_nm': int(raw.uns['resolution'])}
    uns['annotation_info'] = {
        'cell_type': 'RCTD-based cell type (Level 1); cell_type_level2 = subtypes',
        'tumor_layer': '50 um layers by distance from the tumor border (negative = peritumor); NaN = bin removed by QC',
        'distance_to_border_um': 'signed distance from the tumor border (um)',
        'tumor_border': 'bin on the detected tumor border',
        'Lrrc15_border': 'Background / Border / FB_Lrrc15+ (Fig. 3D, Supplementary Fig. S9)',
        'aucell': 'AUCell scores of gene signatures for bins in layers -150 um to 700 um; aucell_YAP_3 = YAP signature in Fig. 3C',
    }
    raw = ad.AnnData(X=raw.X, obs=obs, var=pd.DataFrame(index=raw.var_names), obsm={'spatial': raw.obsm['spatial']}, uns=uns)
    out = f'{out_dir}/{geo}_bin30_annotated.h5ad'
    raw.write(out, compression='gzip')
    print(f'{pid:10s} -> {os.path.basename(out)}  bins={raw.n_obs}  annotated={obs.tumor_layer.notna().sum()}  '
          f'scored={obs.tumor_layer.isin(SCORED_LAYERS).sum()}  aucell={scores.shape[1]}', flush=True)
