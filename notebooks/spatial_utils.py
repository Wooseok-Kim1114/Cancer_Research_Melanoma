"""Helper functions for the Stereo-seq notebooks.

Function bodies are copied from the original analysis notebook. The two versions of
plot_layer_distribution in the original are kept as plot_layer_distribution (Fig. 3E style,
y-label + Mean_expression_*.pdf) and plot_layer_distribution_trend (title + Decreasing_tendency_*.pdf).
"""
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

# GEO library names (groups: Mattress, Running, YAP5SA).
SAMPLE_LIST = ['FF_Mattress_1', 'FF_Mattress_2', 'FF_Mattress_3',
               'FF_Running_1', 'FF_Running_2', 'FF_Running_3',
               'FF_YAP5SA_1', 'FF_YAP5SA_2', 'FF_YAP5SA_3']


# Layers with AUCell scores (from 150 um outside to 700 um inside the tumor border)
SCORED_LAYERS = ['-150um~-100um', '-100um~-50um', '-50um~0um'] + [f'{i}um~{i + 50}um' for i in range(0, 700, 50)]


def read_annotated(annot_dir, sample):
    """bin30 counts (X) with analysis annotations in .obs (<sample>_bin30_annotated.h5ad)."""
    return sc.read(f"{annot_dir}/{sample}_bin30_annotated.h5ad")


def load_margin_lognorm(annot_dir, sample_list=SAMPLE_LIST):
    """QC-passed bins (tumor_layer not NaN) with cell_type / tumor_layer annotations,
    log-normalised on the gene set shared by all samples (as in the original analysis)."""
    adataL = {}
    for sample in sample_list:
        adata = read_annotated(annot_dir, sample)
        adata = adata[adata.obs['tumor_layer'].notna()].copy()
        adata.layers['counts'] = adata.X.copy()
        adataL[sample] = adata
    combined_adata = ad.concat(adataL)
    for sample in sample_list:
        adataL[sample] = adataL[sample][:, combined_adata.var_names]
        del adataL[sample].X
        adataL[sample].X = adataL[sample].layers['counts']
        sc.pp.normalize_total(adataL[sample], target_sum=1e4)
        sc.pp.log1p(adataL[sample])
        adataL[sample].raw = adataL[sample].copy()
    return adataL


def load_scored(annot_dir, sample_list=SAMPLE_LIST):
    """Bins in SCORED_LAYERS (AUCell scores in obs.aucell_*),
    log-normalised on the genes shared by all samples (as in the original analysis)."""
    lognorm4 = {}
    for sample in sample_list:
        adata = read_annotated(annot_dir, sample)
        lognorm4[sample] = adata[adata.obs['tumor_layer'].isin(SCORED_LAYERS)].copy()
    shared_genes = ad.concat(lognorm4).var_names
    for sample in sample_list:
        adata = lognorm4[sample][:, shared_genes].copy()
        adata.layers['counts'] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata.copy()
        lognorm4[sample] = adata
    return lognorm4


FFPE_SAMPLES = ['FFPE_Skin_WT1', 'FFPE_Skin_KO1']


def load_ffpe(annot_dir, sample_list=FFPE_SAMPLES):
    """FFPE slides: QC-passed bins (cell_type not NaN), log-normalised per sample (as in the original analysis)."""
    adataL = {}
    for sample in sample_list:
        adata = read_annotated(annot_dir, sample)
        adata = adata[adata.obs['cell_type'].notna()].copy()
        adata.layers['counts'] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata.copy()
        adataL[sample] = adata
    return adataL


def combine_adatas(adataL, sample_list):
    combined = adataL[sample_list[0]].copy()
    combined.obs['sample'] = sample_list[0]
    
    for sample in sample_list[1:]:
        temp = adataL[sample].copy()
        temp.obs['sample'] = sample
        combined = combined.concatenate(temp, join='outer')
    
    return combined

def plot_combined_group_cell_type_proportions(combined_adata, anno_col='cell_type', group_col='group'):
    # 색상 딕셔너리 정의
    Level1_colors_dict = {
        'B': '#95D5B2',
        'NK': '#2D6A4F',
        'Melanoma': '#DCD694',
        'Endothelial': '#E09399',
        'Fibroblast': '#590D22',
        'Pericyte': '#C26146',
        'DC': '#03045E',
        'Macrophage': '#593670',
        'Monocyte': '#8020FF',
        'Neutrophil': '#5F8E99',
        'T': '#4ED96B'
    }
    
    cell_types = combined_adata.obs[anno_col].unique()
    layers = ['-150um~-100um', '-100um~-50um', '-50um~0um', '0um~50um', 
              '50um~100um', '100um~150um', '150um~200um', '200um~250um', '250um~300um', 
              '300um~350um', '350um~400um', '400um~450um', '450um~500um', '500um~550um', 
              '550um~600um', '600um~650um', '650um~700um']
    groups = combined_adata.obs[group_col].unique()
    
    for group in groups:
        # 각 그룹별로 새로운 figure 생성
        fig, ax = plt.subplots(figsize=(6, 6))
        
        group_adata = combined_adata[combined_adata.obs[group_col] == group]
        samples = group_adata.obs['sample'].unique()
        
        all_proportions = {cell_type: np.zeros((len(samples), len(layers))) for cell_type in cell_types}
        for i, sample in enumerate(samples):
            sample_adata = group_adata[group_adata.obs['sample'] == sample]
            
            for j, layer in enumerate(layers):
                layer_cells = sample_adata.obs[sample_adata.obs['tumor_layer'] == layer]
                total_cells = len(layer_cells)
                for cell_type in cell_types:
                    cell_count = sum(layer_cells[anno_col] == cell_type)
                    all_proportions[cell_type][i, j] = cell_count / total_cells if total_cells > 0 else 0
        
        x_values = range(len(layers))
        for cell_type in cell_types:
            y_values = np.mean(all_proportions[cell_type], axis=0)
            
            # 95% 신뢰 구간 계산
            ci_low = np.percentile(all_proportions[cell_type], 2.5, axis=0)
            ci_high = np.percentile(all_proportions[cell_type], 97.5, axis=0)
            
            # 색상 딕셔너리에서 색상 가져오기
            color = Level1_colors_dict.get(cell_type, '#000000')
            
            # 라인 플롯과 신뢰구간에 같은 색상 적용
            ax.plot(x_values, y_values, label=cell_type, marker='o', color=color)
            ax.fill_between(x_values, ci_low, ci_high, alpha=0.2, color=color)
        
        ax.set_xticks(x_values[::1])
        ax.set_xticklabels([layers[i] for i in range(0, len(layers), 1)], rotation=60, ha='right')
        ax.set_xlabel('')
        ax.set_ylabel('Cell Type Proportion')
        ax.set_title(f'{group}')
        ax.grid(True, linestyle='--', alpha=0.7)
        
        # 범례 추가
        ax.legend(title='Cell Type', bbox_to_anchor=(1.05, 0.5), loc='center left')
        
        # 각 그룹별로 개별 저장
        plt.tight_layout()
        plt.savefig(f"./figures/{group}_celltype_proportions_across_layers.pdf", 
                   bbox_inches='tight')  # bbox_inches='tight'를 추가하여 범례가 잘리지 않도록 함
        plt.show()
        plt.close()  # 메모리 관리를 위해 figure 닫기


def extract_group(sample_name):
    match = re.search(r'(Mattress|Running|YAP5SA)', sample_name)
    return match.group(1) if match else 'Unknown'

def plot_cell_type_proportions_line(proportion_dict, cell_types=None):

    # Define custom colors for groups
    group_dict = {
        'Mattress': '#000000',
        'Running': '#3333FF',
        'YAP5SA': '#77B3FC'
    }

    # Define the correct order of layers
    layers = ['-150um~-100um',
              '-100um~-50um', '-50um~0um', '0um~50um', '50um~100um', '100um~150um',
              '150um~200um', '200um~250um', '250um~300um', '300um~350um',
              '350um~400um', '400um~450um', '450um~500um', '500um~550um',
              '550um~600um', '600um~650um', '650um~700um']

    # Prepare data for plotting
    plot_data = []
    for sample, df in proportion_dict.items():
        group = extract_group(sample)
        # Filter cell types if specified
        if cell_types is not None:
            if isinstance(cell_types, str):
                cell_types = [cell_types]
            df = df.loc[df.index.intersection(cell_types)]
            
        for cell_type in df.index:
            for layer in layers:
                proportion = df.loc[cell_type, layer] if layer in df.columns else np.nan
                plot_data.append({
                    'Cell Type': cell_type,
                    'Layer': layer,
                    'Group': group,
                    'Proportion': proportion,
                    'Sample': sample
                })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Get unique cell types and groups
    plot_cell_types = plot_df['Cell Type'].unique()
    groups = sorted(plot_df['Group'].unique())  # Sort groups to ensure consistent order

    # Plotting
    for cell_type in plot_cell_types:
        fig, ax = plt.subplots(figsize=(6, 6))
        
        for group in groups:
            group_data = plot_df[(plot_df['Cell Type'] == cell_type) & (plot_df['Group'] == group)]
            
            # Calculate mean and standard error
            mean_data = group_data.groupby('Layer')['Proportion'].mean().reset_index()
            sem_data = group_data.groupby('Layer')['Proportion'].sem().reset_index()
            
            # Ensure data is in the correct order
            mean_data['Layer_num'] = mean_data['Layer'].apply(lambda x: layers.index(x))
            mean_data = mean_data.sort_values('Layer_num')
            sem_data['Layer_num'] = sem_data['Layer'].apply(lambda x: layers.index(x))
            sem_data = sem_data.sort_values('Layer_num')
            
            # Line plot with custom colors
            ax.plot(mean_data['Layer_num'], mean_data['Proportion'], 
                   color=group_dict[group], label=group, marker='o')
            ax.fill_between(mean_data['Layer_num'], 
                          mean_data['Proportion'] - sem_data['Proportion'],
                          mean_data['Proportion'] + sem_data['Proportion'],
                          color=group_dict[group], alpha=0.2)

        ax.set_title('')
        ax.set_xlabel('')
        ax.set_ylabel('Cell Type Proportion')
        ax.set_xticks(range(len(layers)))
        ax.set_xticklabels(layers, rotation=45, ha='right')
        
        # Add legend
        ax.legend(title='Group', loc='best')

        # Perform Mann-Whitney U tests and add significance bars
        for layer_idx, layer in enumerate(layers):
            layer_data = plot_df[(plot_df['Cell Type'] == cell_type) & 
                               (plot_df['Layer'] == layer)]
            
            comparisons = [('Mattress', 'Running'), ('Running', 'YAP5SA'), 
                         ('Mattress', 'YAP5SA')]
            y_max = layer_data['Proportion'].max()
            
            for i, (group1, group2) in enumerate(comparisons):
                data1 = layer_data[layer_data['Group'] == group1]['Proportion'].dropna().values
                data2 = layer_data[layer_data['Group'] == group2]['Proportion'].dropna().values
                
                if len(data1) > 0 and len(data2) > 0:
                    stat, p = stats.mannwhitneyu(data1, data2)
                    if p < 0.05:
                        y = y_max + 0.02 + i * 0.02
                        x1 = layer_idx
                        x2 = layer_idx
                        ax.plot([x1, x1, x2, x2], [y, y + 0.01, y + 0.01, y], 
                               lw=1.5, c='k')
                        ax.text((x1 + x2) * 0.5, y + 0.015, f'p = {p:.2e}', 
                               ha='center', va='bottom', fontsize=8)

        plt.tight_layout()
        plt.savefig(f"./figures/{cell_types[0]}_proportion.pdf")
        plt.show()


def plot_layer_distribution(exp_df, gene = None):
    # Define custom colors for groups
    group_dict = {
        'Mattress': '#000000',
        'Running': '#3333FF',
        'YAP5SA': '#77B3FC'
    }
    
    # Get all unique layers (columns)
    layers = exp_df.columns.tolist()
    
    # Prepare data for plotting
    plot_data = []
    for sample in exp_df.index:
        group = extract_group(sample)
        for layer in layers:
            value = exp_df.loc[sample, layer]
            plot_data.append({
                'Layer': layer,
                'Group': group,
                'Value': value,
                'Sample': sample
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Get unique groups
    groups = sorted(plot_df['Group'].unique())
    
    # Create plot
    fig, ax = plt.subplots(figsize=(6, 6))
    
    for group in groups:
        group_data = plot_df[plot_df['Group'] == group]
        
        # Calculate mean and standard error for each layer
        mean_data = group_data.groupby('Layer')['Value'].mean().reset_index()
        sem_data = group_data.groupby('Layer')['Value'].sem().reset_index()
        
        # Ensure data is in the correct order
        mean_data['Layer_num'] = range(len(layers))
        sem_data['Layer_num'] = range(len(layers))
        
        # Line plot with custom colors
        ax.plot(mean_data['Layer_num'], mean_data['Value'], 
               color=group_dict[group], label=group, marker='o')
        ax.fill_between(mean_data['Layer_num'], 
                      mean_data['Value'] - sem_data['Value'],
                      mean_data['Value'] + sem_data['Value'],
                      color=group_dict[group], alpha=0.2)
    
    # Customize plot
    ax.set_xlabel('')
    ax.set_ylabel(f'{gene} mean expression')
    ax.set_xticks(range(len(layers)))
    ax.set_xticklabels([''] * len(layers))  # Empty strings for x-tick labels
    
    # Remove legend
    ax.legend().remove()  # This removes the legend
    
    # Perform Mann-Whitney U tests and add significance bars
    for layer_idx, layer in enumerate(layers):
        layer_data = plot_df[plot_df['Layer'] == layer]
        
        comparisons = [('Mattress', 'Running'), ('Running', 'YAP5SA'), 
                     ('Mattress', 'YAP5SA')]
        y_max = layer_data['Value'].max()
        
        for i, (group1, group2) in enumerate(comparisons):
            data1 = layer_data[layer_data['Group'] == group1]['Value'].dropna().values
            data2 = layer_data[layer_data['Group'] == group2]['Value'].dropna().values
            
            if len(data1) > 0 and len(data2) > 0:
                stat, p = stats.mannwhitneyu(data1, data2)
                if p < 0.05:
                    y = y_max + 0.02 + i * 0.02
                    x1 = layer_idx
                    x2 = layer_idx
                    ax.plot([x1, x1, x2, x2], [y, y + 0.01, y + 0.01, y], 
                           lw=1.5, c='k')
                    ax.text((x1 + x2) * 0.5, y + 0.015, f'p = {p:.2e}', 
                           ha='center', va='bottom', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(f"./figures/Mean_expression_{gene}.pdf")
    plt.show()


def calculate_group_signature_scores(adata_dict, layers, signature, celltype='Melanoma'):
    # Define the correct order of layers
    layer_order = layers
    
    # Validate that all input layers are in the defined order
    if not all(layer in layer_order for layer in layers):
        raise ValueError("Some input layers are not in the predefined layer order")
    
    results = []
    for sample, adata in adata_dict.items():
        # Filter for melanocytes and specified layers
        mask = (adata.obs['cell_type'] == celltype) & (adata.obs['tumor_layer'].isin(layers))
        
        if mask.sum() > 0:
            group = adata.obs.loc[mask, 'group'].iloc[0]  # Assuming group is consistent for a sample
            for layer in layers:
                layer_mask = mask & (adata.obs['tumor_layer'] == layer)
                if layer_mask.sum() > 0:
                    mean_score = adata.obs.loc[layer_mask, signature].mean()
                    results.append({
                        'Group': group,
                        'Layer': layer,
                        f'Mean {signature} Score': mean_score
                    })
    
    df = pd.DataFrame(results)
    
    # Convert Layer column to categorical with specific order
    df['Layer'] = pd.Categorical(df['Layer'], categories=layer_order, ordered=True)
    
    # Calculate the average score for each group and layer
    group_avg = df.groupby(['Group', 'Layer'])[f'Mean {signature} Score'].mean().reset_index()
    
    return group_avg

def plot_group_yap_scores(df, signature, layer_order, path = None):
    # Define custom colors for groups
    group_dict = {
        'Mattress': '#000000',
        'Running': '#3333FF',
        'YAP5SA': '#77B3FC'
    }
    
    plt.figure(figsize=(6, 6))
    
    # Create line plot with custom colors
    ax = sns.lineplot(
        data=df, 
        x='Layer', 
        y=f'Mean {signature} Score', 
        hue='Group',
        marker='o',
        palette=group_dict  # Use custom color palette
    )
    
    # Set x-axis range (from 0 to number of layers - 1)
    ax.set_xlim(-0.2, len(layer_order) - 0.8)
    
    # Set titles and labels
    plt.title(f'{signature}')
    plt.xlabel('')
    plt.ylabel('Mean AUCell scores')
    
    # Set x-axis tick labels
    ax.set_xticks(range(len(layer_order)))
    ax.set_xticklabels(layer_order, rotation=45, ha='right')
    
    # Add grid and legend
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(title='', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Save the figure if path is provided
    if path is not None:
        plt.savefig(path, bbox_inches='tight', dpi=300)
    
    plt.show()


def plot_layer_distribution_trend(exp_df, gene = None):
    # Define custom colors for groups
    group_dict = {
        'Mattress': '#000000',
        'Running': '#3333FF',
        'YAP5SA': '#77B3FC'
    }
    
    # Get all unique layers (columns)
    layers = exp_df.columns.tolist()
    
    # Prepare data for plotting
    plot_data = []
    for sample in exp_df.index:
        group = extract_group(sample)
        for layer in layers:
            value = exp_df.loc[sample, layer]
            plot_data.append({
                'Layer': layer,
                'Group': group,
                'Value': value,
                'Sample': sample
            })
    
    plot_df = pd.DataFrame(plot_data)
    
    # Get unique groups
    groups = sorted(plot_df['Group'].unique())
    
    # Create plot
    fig, ax = plt.subplots(figsize=(6, 6))
    
    for group in groups:
        group_data = plot_df[plot_df['Group'] == group]
        
        # Calculate mean and standard error for each layer
        mean_data = group_data.groupby('Layer')['Value'].mean().reset_index()
        sem_data = group_data.groupby('Layer')['Value'].sem().reset_index()
        
        # Ensure data is in the correct order
        mean_data['Layer_num'] = range(len(layers))
        sem_data['Layer_num'] = range(len(layers))
        
        # Line plot with custom colors
        ax.plot(mean_data['Layer_num'], mean_data['Value'], 
               color=group_dict[group], label=group, marker='o')
        ax.fill_between(mean_data['Layer_num'], 
                      mean_data['Value'] - sem_data['Value'],
                      mean_data['Value'] + sem_data['Value'],
                      color=group_dict[group], alpha=0.2)
    
    # Customize plot
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(range(len(layers)))
    ax.set_xticklabels([''] * len(layers))  # Empty strings for x-tick labels
    
    # Remove legend
    ax.legend().remove()  # This removes the legend
    
    # Perform Mann-Whitney U tests and add significance bars
    for layer_idx, layer in enumerate(layers):
        layer_data = plot_df[plot_df['Layer'] == layer]
        
        comparisons = [('Mattress', 'Running'), ('Running', 'YAP5SA'), 
                     ('Mattress', 'YAP5SA')]
        y_max = layer_data['Value'].max()
        
        for i, (group1, group2) in enumerate(comparisons):
            data1 = layer_data[layer_data['Group'] == group1]['Value'].dropna().values
            data2 = layer_data[layer_data['Group'] == group2]['Value'].dropna().values
            
            if len(data1) > 0 and len(data2) > 0:
                stat, p = stats.mannwhitneyu(data1, data2)
                if p < 0.05:
                    y = y_max + 0.02 + i * 0.02
                    x1 = layer_idx
                    x2 = layer_idx
                    ax.plot([x1, x1, x2, x2], [y, y + 0.01, y + 0.01, y], 
                           lw=1.5, c='k')
                    ax.text((x1 + x2) * 0.5, y + 0.015, f'p = {p:.2e}', 
                           ha='center', va='bottom', fontsize=8)
    
    plt.title(f'{gene}', fontsize = 12)
    plt.tight_layout()
    plt.savefig(f"./figures/Decreasing_tendency_{gene}.pdf")
    plt.show()


def plot_gene_group_heatmaps(adata, gene_list, group_col='group', sample_col='sample', 
                           term=None, genes_per_plot=50):
    """
    여러 히트맵을 하나의 페이지에 가로로 배치하여 PDF로 저장하는 함수
    
    Parameters
    -----------
    adata : AnnData object
        분석할 데이터
    gene_list : list
        표시할 유전자 목록
    group_col : str
        그룹 정보가 있는 컬럼명
    sample_col : str
        샘플 정보가 있는 컬럼명
    term : str
        히트맵 제목 및 파일명에 사용될 문자열
    genes_per_plot : int
        각 부분 히트맵당 표시할 유전자 수
    """
    
    # 기본 term 설정
    if term is None:
        term = 'gene_heatmaps'
    
    # 파일명 설정
    output_pdf = f'Signature_gene_expression_{term}.pdf'
    
    # 존재하는 유전자만 필터링
    if not all(gene in adata.var_names for gene in gene_list):
        missing_genes = [gene for gene in gene_list if gene not in adata.var_names]
        print(f"Warning: Following genes not found: {missing_genes}")
        gene_list = [gene for gene in gene_list if gene in adata.var_names]
    
    # 전체 유전자 발현량 계산
    gene_exp = pd.DataFrame(adata.raw[:, gene_list].X.toarray(), 
                          columns=gene_list, 
                          index=adata.obs.index)
    
    # 그룹별 평균 발현량 계산
    group_means = gene_exp.join(adata.obs[group_col]).groupby(group_col).mean()
    
    # Z-score 변환
    z_scored = pd.DataFrame(data=stats.zscore(group_means, axis=0),
                          columns=gene_list,
                          index=group_means.index)
    
    # Tread와 YAP5SA 그룹의 발현량을 기준으로 유전자 정렬
    if 'Running' in z_scored.index and 'YAP5SA' in z_scored.index:
        gene_order = z_scored.T.sort_values(by=['Running', 'YAP5SA'], ascending=False).index.tolist()
        z_scored = z_scored[gene_order]
    
    # 정렬된 유전자 리스트를 지정된 크기로 분할
    gene_chunks = [gene_order[i:i + genes_per_plot] 
                  for i in range(0, len(gene_order), genes_per_plot)]
    
    # 전체 figure 크기 계산
    n_chunks = len(gene_chunks)
    width_per_chunk = 3  # 각 히트맵의 기본 너비
    total_width = n_chunks * width_per_chunk
    height = len(gene_chunks[0]) / 2  # 유전자 수에 비례한 높이
    
    # PDF 생성
    with PdfPages(output_pdf) as pdf:
        # 하나의 큰 figure 생성
        fig = plt.figure(figsize=(total_width + 1, height))  # 여백을 위해 1 추가
        
        # GridSpec을 사용하여 subplot 배치
        gs = gridspec.GridSpec(1, n_chunks, width_ratios=[1] * n_chunks)
        gs.update(wspace=0.3)  # subplot 간 간격 조정
        
        # 각 청크에 대한 히트맵 그리기
        for i, chunk in enumerate(gene_chunks):
            # subplot 생성
            ax = plt.subplot(gs[i])
            
            # 데이터 준비
            chunk_data = z_scored[chunk].T
            
            # 히트맵 그리기
            sns.heatmap(chunk_data,
                       cmap='RdBu_r',
                       center=0,
                       annot=False,
                       fmt='.2f',
                       cbar=True if i == n_chunks-1 else False,  # 마지막 히트맵에만 컬러바 표시
                       cbar_kws={'label': 'Z-score'} if i == n_chunks-1 else None,
                       linewidths=0,
                       square=True,
                       ax=ax)
            
            plt.grid(False)
        
        # 전체 제목 추가
        plt.suptitle(term, y=1.02)
        plt.tight_layout()
        
        # PDF로 저장
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    print(f"Heatmap saved to {output_pdf}")
    return z_scored


def plot_gene_group_heatmaps_raw(adata, gene_list, group_col='group', sample_col='sample', 
                           term=None, genes_per_plot=50):
    """
    여러 히트맵을 하나의 페이지에 가로로 배치하여 PDF로 저장하는 함수 (Raw expression values)
    
    Parameters
    -----------
    adata : AnnData object
        분석할 데이터
    gene_list : list
        표시할 유전자 목록
    group_col : str
        그룹 정보가 있는 컬럼명
    sample_col : str
        샘플 정보가 있는 컬럼명
    term : str
        히트맵 제목 및 파일명에 사용될 문자열
    genes_per_plot : int
        각 부분 히트맵당 표시할 유전자 수
    """
    
    # 기본 term 설정
    if term is None:
        term = 'gene_heatmaps'
    
    # 파일명 설정
    output_pdf = f'Signature_raw_gene_expression_{term}.pdf'
    
    # 존재하는 유전자만 필터링
    if not all(gene in adata.var_names for gene in gene_list):
        missing_genes = [gene for gene in gene_list if gene not in adata.var_names]
        print(f"Warning: Following genes not found: {missing_genes}")
        gene_list = [gene for gene in gene_list if gene in adata.var_names]
    
    # 전체 유전자 발현량 계산
    gene_exp = pd.DataFrame(adata.raw[:, gene_list].X.toarray(), 
                          columns=gene_list, 
                          index=adata.obs.index)
    
    # 그룹별 평균 발현량 계산
    group_means = gene_exp.join(adata.obs[group_col]).groupby(group_col).mean()
    
    # Tread와 YAP5SA 그룹의 발현량을 기준으로 유전자 정렬
    if 'Running' in group_means.index and 'YAP5SA' in group_means.index:
        gene_order = group_means.T.sort_values(by=['Running', 'YAP5SA'], ascending=False).index.tolist()
        group_means = group_means[gene_order]
    
    # 정렬된 유전자 리스트를 지정된 크기로 분할
    gene_chunks = [gene_order[i:i + genes_per_plot] 
                  for i in range(0, len(gene_order), genes_per_plot)]
    
    # 전체 figure 크기 계산
    n_chunks = len(gene_chunks)
    width_per_chunk = 3  # 각 히트맵의 기본 너비
    total_width = n_chunks * width_per_chunk
    height = len(gene_chunks[0]) / 2  # 유전자 수에 비례한 높이
    
    # PDF 생성
    with PdfPages(output_pdf) as pdf:
        # 하나의 큰 figure 생성
        fig = plt.figure(figsize=(total_width + 1, height))  # 여백을 위해 1 추가
        
        # GridSpec을 사용하여 subplot 배치
        gs = gridspec.GridSpec(1, n_chunks, width_ratios=[1] * n_chunks)
        gs.update(wspace=0.3)  # subplot 간 간격 조정
        
        # 각 청크에 대한 히트맵 그리기
        for i, chunk in enumerate(gene_chunks):
            # subplot 생성
            ax = plt.subplot(gs[i])
            
            # 데이터 준비
            chunk_data = group_means[chunk].T
            
            # 히트맵 그리기
            sns.heatmap(chunk_data,
                       cmap='RdBu_r',
                       center=None,  # z-score가 아니므로 center 제거
                       annot=False,
                       fmt='.2f',
                       cbar=True if i == n_chunks-1 else False,  # 마지막 히트맵에만 컬러바 표시
                       cbar_kws={'label': 'Expression'} if i == n_chunks-1 else None,
                       linewidths=0,
                       square=True,
                       ax=ax)
            
            plt.grid(False)
        
        # 전체 제목 추가
        plt.suptitle(term, y=1.02)
        plt.tight_layout()
        
        # PDF로 저장
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()
    
    print(f"Heatmap saved to {output_pdf}")
    return group_means


# Line plot used for Fig. 3C / Supplementary Fig. S7E (241223_Tiroshi_Module_scoring):
# group mean of the per-sample layer means, shaded = mean +/- SEM across samples
def plot_layer_profile(exp_df, gene=None, figpath=None):
    # Define custom colors for groups
    group_dict = {
        'Mattress': '#000000',
        'Running': '#3333FF',
        'YAP5SA': '#77B3FC'
    }
    
    # Get layers and sort them properly
    def extract_first_number(x):
        return int(x.split('~')[0].replace('um',''))
    
    # Sort layers based on their starting number
    layers = sorted(exp_df.columns.tolist(), key=extract_first_number)
    
    # Rest of the code remains the same
    plot_data = []
    for sample in exp_df.index:
        group = extract_group(sample)
        for layer in layers:  # Now using sorted layers
            value = exp_df.loc[sample, layer]
            plot_data.append({
                'Layer': layer,
                'Group': group,
                'Value': value,
                'Sample': sample
            })
    
    plot_df = pd.DataFrame(plot_data)
    groups = sorted(plot_df['Group'].unique())
    
    fig, ax = plt.subplots(figsize=(6, 6))
    
    for group in groups:
        group_data = plot_df[plot_df['Group'] == group]
        
        # Calculate mean and standard error for each layer
        mean_data = group_data.groupby('Layer')['Value'].mean().reset_index()
        sem_data = group_data.groupby('Layer')['Value'].sem().reset_index()
        
        # Ensure data is in the correct order by mapping to sorted layers
        mean_data['Layer_num'] = mean_data['Layer'].map({layer: i for i, layer in enumerate(layers)})
        sem_data['Layer_num'] = sem_data['Layer'].map({layer: i for i, layer in enumerate(layers)})
        
        # Sort by Layer_num to ensure correct plotting order
        mean_data = mean_data.sort_values('Layer_num')
        sem_data = sem_data.sort_values('Layer_num')
        
        ax.plot(mean_data['Layer_num'], mean_data['Value'], 
               color=group_dict[group], label=group, marker='o')
        ax.fill_between(mean_data['Layer_num'], 
                       mean_data['Value'] - sem_data['Value'],
                       mean_data['Value'] + sem_data['Value'],
                       color=group_dict[group], alpha=0.2)
    
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_xticks(range(len(layers)))
    ax.set_xticklabels(layers, rotation=45, ha='right')
    
    ax.legend().remove()
    
    plt.title(f'{gene}', fontsize=12)
    plt.tight_layout()
    if figpath is not None:
        plt.savefig(figpath)
    plt.show()


import matplotlib

# ---------------------------------------------------------------- manuscript-style panels added for the revision
GROUP_COLORS = {'Mattress': '#000000', 'Running': '#3333FF', 'YAP5SA': '#77B3FC'}
REGION_LAYERS = {
    'Peritumor': ['-150um~-100um', '-100um~-50um', '-50um~0um'],
    'Margin': ['0um~50um', '50um~100um', '100um~150um'],
    'Center': ['150um~200um', '200um~250um', '250um~300um', '300um~350um', '350um~400um', '400um~450um',
               '450um~500um', '500um~550um', '550um~600um', '600um~650um', '650um~700um'],
}


def region_mean_expression(lognorm4, gene, cell_type=None, sample_list=SAMPLE_LIST):
    """Mean log-normalised expression per sample and 50 um layer (all bins, or bins of one cell type),
    then the unweighted mean of the layer means per region (Spatial_Revision2.ipynb / CXCL.ipynb).
    Returns a long data frame with columns Sample, Group, Region, Value."""
    layers = sum(REGION_LAYERS.values(), [])
    exp_df = pd.DataFrame(index=sample_list, columns=layers, dtype=float)
    for sample in sample_list:
        adata = lognorm4[sample]
        for layer in layers:
            mask = adata.obs['tumor_layer'] == layer
            if cell_type is not None:
                mask &= adata.obs['cell_type'] == cell_type
            exp_df.loc[sample, layer] = np.mean(adata[mask].raw[:, gene].X)
    rows = []
    for sample in sample_list:
        for region, region_layers in REGION_LAYERS.items():
            rows.append({'Sample': sample, 'Group': extract_group(sample), 'Region': region,
                         'Value': exp_df.loc[sample, region_layers].mean()})
    return pd.DataFrame(rows)


def tukey_region_group(box_df):
    """One-way ANOVA over all Region x Group combinations followed by Tukey's HSD (as in the original analysis)."""
    from scipy.stats import f_oneway
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    combined = box_df['Region'] + '_' + box_df['Group']
    f_stat, p_val = f_oneway(*[box_df.loc[combined == g, 'Value'] for g in combined.unique()])
    print(f'One-way ANOVA: F = {f_stat:.3f}, p = {p_val:.3g}')
    tukey = pairwise_tukeyhsd(endog=box_df['Value'], groups=combined, alpha=0.05)
    res = pd.DataFrame(tukey._results_table.data[1:], columns=tukey._results_table.data[0])
    same_region = res['group1'].str.split('_').str[0] == res['group2'].str.split('_').str[0]
    return res[same_region].reset_index(drop=True)


def plot_region_box(box_df, title, ylabel='Mean expression', figsize=(4, 4)):
    """Box plots per region and group (Fig. 3E, Supplementary Fig. S8A): median, IQR, whiskers 1.5 x IQR."""
    order = list(REGION_LAYERS)
    hue_order = ['Mattress', 'Running', 'YAP5SA']
    palette = {'Mattress': '#FFFFFF', 'Running': '#3333FF', 'YAP5SA': '#77B3FC'}
    fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(data=box_df, x='Region', y='Value', hue='Group', order=order, hue_order=hue_order,
                palette=palette, saturation=1, fliersize=2, linewidth=1, ax=ax,
                boxprops={'edgecolor': 'black'}, medianprops={'color': 'black'},
                whiskerprops={'color': 'black'}, capprops={'color': 'black'})
    ax.legend(title=None, frameon=False, bbox_to_anchor=(1, 1), loc='upper left')
    ax.set_title(title, style='italic')
    ax.set_ylabel(ylabel)
    ax.set_xlabel('')
    ax.grid(False)
    return fig, ax


def rotate_section(adata, angle):
    """Coordinates of the section as displayed by sc.pl.spatial (y axis pointing down), rotated counter-clockwise
    by `angle` degrees; returns (x, y) with y pointing up."""
    xy = adata.obsm['spatial'].astype(float)
    x, y = xy[:, 0], -xy[:, 1]
    t = np.deg2rad(angle)
    return x * np.cos(t) - y * np.sin(t), x * np.sin(t) + y * np.cos(t)


def plot_spatial_zoom(adata, color, angle, xlim, ylim, spot_size=50, title='', **kwargs):
    """Zoom-in of a rotated section (Fig. 4E, Supplementary Fig. S11B). The sections in the manuscript were rotated
    and cropped by hand; `angle` / `xlim` / `ylim` (rotated coordinates) approximate those panels."""
    x, y = rotate_section(adata, angle)
    keep = (x >= xlim[0]) & (x <= xlim[1]) & (y >= ylim[0]) & (y <= ylim[1])
    bdata = adata[keep].copy()
    bdata.obsm['spatial'] = np.column_stack([x[keep], -y[keep]])
    fig, ax = plt.subplots(figsize=(4 * (xlim[1] - xlim[0]) / (ylim[1] - ylim[0]), 4))
    sc.pl.spatial(bdata, color=color, spot_size=spot_size, frameon=False, title=title, ax=ax, show=False, **kwargs)
    ax.set_xlim(*xlim)
    ax.set_ylim(-ylim[0], -ylim[1])
    ax.set_aspect('equal')
    ax.axis('on')                      # frame around the zoom-in, as in the manuscript
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('black')
    return fig, ax


def plot_hbar(values, labels, color, xlabel, title='', xerr=None, figsize=(4, 3), top_axis=True):
    """Horizontal bar plot in the manuscript style (first value on top)."""
    fig, ax = plt.subplots(figsize=figsize)
    pos = np.arange(len(values))[::-1]
    ax.barh(pos, values, xerr=xerr, color=color, edgecolor='black', linewidth=0.5, height=0.7,
            error_kw={'capsize': 3, 'elinewidth': 0.8})
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


def plot_group_zscore_heatmap(expr, sample_groups, figsize=(8, 2.2), vmax=2):
    """Genes (columns) x samples (rows) heatmap of per-gene z-scores (Fig. 3F), with a group colour bar."""
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1, ddof=0), axis=0)
    fig = plt.figure(figsize=figsize)
    ax_bar = fig.add_axes([0.0, 0.1, 0.012, 0.8])
    ax = fig.add_axes([0.015, 0.1, 0.85, 0.8])
    cax = fig.add_axes([0.88, 0.3, 0.012, 0.4])
    sns.heatmap(z.T, cmap='RdBu_r', vmin=-vmax, vmax=vmax, xticklabels=False, yticklabels=False, ax=ax,
                cbar_ax=cax, cbar_kws={'label': 'z-score'})
    colors = [GROUP_COLORS[g] for g in sample_groups]
    ax_bar.imshow(np.array([matplotlib.colors.to_rgb(c) for c in colors])[:, None, :], aspect='auto')
    ax_bar.set_xticks([])
    ax_bar.set_yticks([])
    ax_bar.set_ylabel('Peritumor')
    return fig, ax


def significance_stars(p):
    return '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
