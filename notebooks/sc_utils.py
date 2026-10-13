"""Helper functions for the scRNA-seq notebooks (copied from the original analysis notebook)."""
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt


def toCount_HVG_Scale(adata, cell_cycle_file='../data/regev_lab_cell_cycle_genes.txt'):
    adata.X = adata.layers['counts']
    adata.layers['counts'] = adata.X.copy() 
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata.copy()

    # HVG
    sc.pp.highly_variable_genes(adata)
    print(f"detected {np.sum(adata.var['highly_variable'])} highly variable genes")
    print("="*100)
    
    # scaling
    print('scaling ...')
    sc.pp.scale(adata, max_value=10)
    print('max: ', np.max(adata.X))
    print('min: ', np.min(adata.X))
    print('mean: ', np.mean(adata.X))
    print("="*100)
    
    # cell cycle score
    cell_cycle_genes = [x.strip().lower().capitalize()
                        for x in open(cell_cycle_file)]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    print("="*100)
    
    # remove cell cycle genes in HVG
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    print('remove {} cell cycle genes in HVG ...'.format(adata.var[adata.var_names.isin(cell_cycle_genes)]['highly_variable'].value_counts()[True]))
    adata.var.loc[cell_cycle_genes, 'highly_variable'] = False
    print("="*100)

    return adata


def _bw_nrd0(x):
    # stats::bw.nrd0
    hi = np.std(x, ddof=1)
    lo = min(hi, (np.quantile(x, 0.75) - np.quantile(x, 0.25)) / 1.349)
    if not lo:
        lo = hi or abs(x[0]) or 1.0
    return 0.9 * lo * len(x) ** (-0.2)


def _r_density(x, bw, n=512):
    """stats::density(x, bw, from = min(x), to = max(x), n = 512) (gaussian kernel, linear binning + FFT)."""
    x = np.asarray(x, dtype=float)
    frm, to = x.min(), x.max()
    lo, up = frm - 4 * bw, to + 4 * bw
    # BinDist
    y = np.zeros(2 * n)
    xdelta = (up - lo) / (n - 1)
    xpos = (x - lo) / xdelta
    ix = np.floor(xpos).astype(int)
    fx = xpos - ix
    w = 1.0 / len(x)
    inside = (ix >= 0) & (ix <= n - 2)
    np.add.at(y, ix[inside], w * (1 - fx[inside]))
    np.add.at(y, ix[inside] + 1, w * fx[inside])
    np.add.at(y, np.zeros((ix == -1).sum(), dtype=int), w * fx[ix == -1])
    np.add.at(y, ix[ix == n - 1], w * (1 - fx[ix == n - 1]))
    # convolution with the gaussian kernel
    kords = np.linspace(0, 2 * (up - lo), 2 * n)
    kords[n + 1:] = -kords[n - 1:0:-1]
    kords = np.exp(-0.5 * (kords / bw) ** 2) / (bw * np.sqrt(2 * np.pi))
    conv = np.fft.ifft(np.fft.fft(y) * np.conj(np.fft.fft(kords)))
    dens = np.maximum(0, conv.real[:n])
    xords = np.linspace(lo, up, n)
    grid = np.linspace(frm, to, n)
    return grid, np.interp(grid, xords, dens)


def plot_group_violin(values, groups, title, ylabel='score', order=('Mattress', 'Running'),
                      colors=('#dcdcdc', '#3333ff'), figsize=(3, 4)):
    """Violin plot reproducing the ggplot2 panels of Fig. 2H/2I
    (geom_violin(draw_quantiles = 0.5) + theme_classic(), 3 x 4 in): density as in stat_ydensity
    (trim = TRUE, bw.nrd0, n = 512), scale = 'area', width 0.9, and the 0.5 quantile line taken from
    the cumulative estimated density, as ggplot2 does."""
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups, dtype=str)
    dens = {}
    for g in order:
        x = values[groups == g]
        dens[g] = _r_density(x, _bw_nrd0(x))
    dmax = max(d.max() for _, d in dens.values())

    lw = 0.5 * 72.27 / 25.4 * 0.75          # ggplot2 linewidth 0.5 (mm) in points
    fig, ax = plt.subplots(figsize=figsize)
    for i, (g, col) in enumerate(zip(order, colors)):
        y, d = dens[g]
        half = d / dmax * 0.45               # scale = 'area', width = 0.9
        ax.fill(np.r_[i - half, (i + half)[::-1]], np.r_[y, y[::-1]], facecolor=col, edgecolor='#333333',
                linewidth=lw, zorder=2)
        cum = np.cumsum(d) / d.sum()
        q = np.interp(0.5, cum, y)            # draw_quantiles = 0.5
        hq = np.interp(q, y, half)
        ax.plot([i - hq, i + hq], [q, q], color='#333333', linewidth=lw, zorder=3)
    ax.set_xticks(range(len(order)), order)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.grid(False)
    ax.set_facecolor('none')
    ax.spines[['top', 'right']].set_visible(False)
    for sp in ['left', 'bottom']:
        ax.spines[sp].set_color('black')
        ax.spines[sp].set_linewidth(lw)
    ax.tick_params(colors='#333333', labelcolor='#4d4d4d', labelsize=8.8, width=lw, length=2.75)
    ax.set_xlabel('')
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=14)
    return fig, ax


def plot_hbar(values, labels, color, xlabel, title='', figsize=(4, 3), top_axis=True):
    """Horizontal bar plot in the manuscript style (first value on top; Supplementary Fig. S6F/G)."""
    fig, ax = plt.subplots(figsize=figsize)
    pos = np.arange(len(values))[::-1]
    ax.barh(pos, values, color=color, edgecolor='black', linewidth=0.5, height=0.7)
    ax.set_yticks(pos, labels)
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.grid(False)
    if top_axis:
        ax.xaxis.set_label_position('top')
        ax.xaxis.tick_top()
        ax.spines[['bottom', 'right']].set_visible(False)
    else:
        ax.spines[['top', 'right']].set_visible(False)
    return fig, ax


def signature_group_test(adata, genesets, terms, case='Tread', control='Sponge'):
    """AUCell score per cell (decoupler, min_n=1; scores from adata.raw), mean per sample, and a two-sided
    Student's t-test of the sample means (case vs control), as in 241031_{Mel,Fibro}_Scoring.ipynb."""
    import pandas as pd
    import decoupler as dc
    from scipy.stats import ttest_ind
    for key in terms:
        genes = list(set(x for x in genesets[key] if x in adata.var_names))
        dc.run_aucell(mat=adata, net=pd.DataFrame({'target': genes, 'source': key}), source='source',
                      target='target', verbose=False, min_n=1)
        adata.obs[key] = adata.obsm['aucell_estimate'][key].values
    per_sample = adata.obs.groupby('sample', observed=True)[terms].mean()
    group = adata.obs.groupby('sample', observed=True)['group'].first().astype(str)
    a, b = per_sample.loc[group == case], per_sample.loc[group == control]
    return pd.DataFrame({'log2FC': np.log2(a.mean() / b.mean()),
                         '-log10P': -np.log10(ttest_ind(a, b)[1])}, index=terms)
