import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def spatial(adata):
    img = plt.imread(adata.uns['lung_image'])
    fig, ax = plt.subplots()
    ax.imshow(img, extent=adata.uns['extent'])

    plt.scatter(adata.obsm['spatial'][:,0],
                adata.obsm['spatial'][:,1],
                s=10,
                c='cornflowerblue')
    ax.invert_yaxis()
    ax.axis('off')


def spatial_macrophage (adata):
    macrophage_markers = ['Cd14', 'Cd68', 'Adgre1']
    mask1 = np.isin(adata.var_names, macrophage_markers)
    total_counts = np.array(adata.X[:, mask1].sum(axis=1))
    mask = (total_counts > 0)
    custom_palette = []
    for c in mask:
        if c:
            custom_palette.append('green')
        else:
            custom_palette.append('cornflowerblue')
    img = plt.imread(adata.uns['lung_image'])
    fig, ax = plt.subplots()
    ax.imshow(img, extent=adata.uns['extent'])
    plt.scatter(adata.obsm['spatial'][:,0],
                adata.obsm['spatial'][:,1],
                s=10,
                c=custom_palette)
    ax.invert_yaxis()
    ax.axis('off')


def total_plot(adata, lstGenes=[]):
    if len(lstGenes) > 0:
        mask = np.isin(adata.var_names, lstGenes)
        total_counts = np.array(adata.X[:, mask].sum(axis=1))
    else:
        total_counts = adata.obs['total_counts']
        
    img = plt.imread(adata.uns['lung_image'])
    fig, ax = plt.subplots()
    ax.imshow(img, extent=adata.uns['extent'])
    plt.scatter(adata[total_counts > 0].obsm['spatial'][:,0],
                adata[total_counts > 0].obsm['spatial'][:,1],
                s=10,
                c=total_counts[total_counts > 0],
                cmap='Reds')
    ax.invert_yaxis()

    ax.axis('off')
    ax.set_xlabel('Spatial_1')
    ax.set_ylabel('Spatial_2')
    plt.colorbar(label='Expression')


def gene_count_plot(adata):
    img = plt.imread(adata.uns['lung_image'])
    fig, ax = plt.subplots()
    ax.imshow(img, extent=adata.uns['extent'])
    plt.scatter(adata.obsm['spatial'][:,0],
                adata.obsm['spatial'][:,1],
                s=10,
                c=adata.obs['n_genes_by_counts'],
                cmap='Reds')
    ax.invert_yaxis()

    ax.axis('off')
    ax.set_xlabel('Spatial_1')
    ax.set_ylabel('Spatial_2')
    plt.colorbar(label='Gene Counts')


def QC_plot(adata):
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))

    p1 = sns.histplot(adata.obs['total_counts'], kde=True, ax=axes[0, 0])
    p1.set_xlabel('Total Counts')
    p1.set_ylabel('Frequency')

    p2 = sns.histplot(adata.obs['n_genes_by_counts'], kde=True, ax=axes[0, 1])
    p2.set_xlabel('Number of Genes')
    p2.set_ylabel('Frequency')

    count_data = adata.obs['total_counts'].copy()
    count_data.sort_values(inplace=True, ascending=False)
    order = range(1, len(count_data)+1)
    axes[1, 0].semilogy(order, count_data, 'b-')
    axes[1, 0].set_xlabel('Spot Rank')
    axes[1, 0].set_ylabel('Total Counts')

    p4 = sns.scatterplot(x='total_counts', y ='n_genes_by_counts', hue='pct_counts_mt', data=adata.obs, ax=axes[1, 1])
    p4.set_xlabel('Total Counts')
    p4.set_ylabel('Number of Genes')

    plt.show()


def plt_umap_cluster (adata):
    color_code = dict(zip(range(0, 10), plt.cm.tab10(range(0, 10))))
    plt.figure(figsize=(10,8))
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_ticks([])
    frame1.axes.get_yaxis().set_ticks([])
    for g in np.unique(adata.obs['clusters']):
        plt.scatter(adata[adata.obs['clusters'] == g].obsm['X_umap'][:, 0],
                    adata[adata.obs['clusters'] == g].obsm['X_umap'][:, 1],
                    color=adata[adata.obs['clusters'] == g].obs['clusters'].astype(int).map(color_code),
                    label=g)
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.legend()
    plt.show()


def plt_spatial_cluster (adata, sample_code : str, sample_name : str):
    '''
    This function takes 3 arguments:
    adata : the annotated data container
    sample_code : one of these sample codes: M24, F31, M63, F62
    sample_name : the name of the sample: 0_Per_1_M24, 0_Per_2_F31, CTL_1_M63, CTL_2_F62
    '''
    color_code = dict(zip(range(0, 10), plt.cm.tab10(range(0, 10))))
    img = plt.imread(adata.uns['lung_image_' + sample_code])
    fig, ax = plt.subplots()
    ax.imshow(img, extent=adata.uns['extent_' + sample_code])
    plt.scatter(adata[adata.obs['sample'] == sample_name].obsm['spatial'][:,0],
                adata[adata.obs['sample'] == sample_name].obsm['spatial'][:,1],
                s=10,
                c=adata[adata.obs['sample'] == sample_name].obs['clusters'].astype(int).map(color_code))
    ax.invert_yaxis()
    for axis in ['bottom', 'left', 'top', 'right']:
            ax.spines[axis].set_visible(False)

    plt.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=False,
        labelbottom=False)

    plt.tick_params(
        axis='y',
        which='both',
        left=False,
        right=False,
        labelleft=False)

    ax.set_xlabel('Spatial_1')
    ax.set_ylabel('Spatial_2')


def plt_umap_sample (adata):
    plt.figure(figsize=(10,8))
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_ticks([])
    frame1.axes.get_yaxis().set_ticks([])
    sample_color_code = {'0_Per_1_M24' : plt.cm.tab20b(range(0,20)[12]),
                         '0_Per_2_F31' : plt.cm.tab20b(range(0,20)[15]),
                         'CTL_1_M63' : plt.cm.tab20b(range(0,20)[8]),
                         'CTL_2_F62' : plt.cm.tab20b(range(0,20)[11])}
    for s in np.unique(adata.obs['sample']):
        plt.scatter(adata[adata.obs['sample'] == s].obsm['X_umap'][:, 0],
                    adata[adata.obs['sample'] == s].obsm['X_umap'][:, 1],
                    color=sample_color_code[s], label=s)
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.legend()
    plt.show()


def plt_umap_sex (adata):
    plt.figure(figsize=(10,8))
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_ticks([])
    frame1.axes.get_yaxis().set_ticks([])
    sex_color_code = {'female' : plt.cm.Accent(range(0,20)[5]),
                      'male' : plt.cm.tab20b(range(0,20)[6])}
    for s in np.unique(adata.obs['sex']):
        plt.scatter(adata[adata.obs['sex'] == s].obsm['X_umap'][:, 0],
                    adata[adata.obs['sex'] == s].obsm['X_umap'][:, 1],
                    color=sex_color_code[s], label=s)
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.legend()
    plt.show()


def plt_umap_condition (adata):
    plt.figure(figsize=(10,8))
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_ticks([])
    frame1.axes.get_yaxis().set_ticks([])
    lung_color_code = {'Z_Per' : plt.cm.tab20b(range(0,20)[17]),
                      'CTL' : plt.cm.tab20c(range(0,20)[9])}
    for l in np.unique(adata.obs['lung']):
        plt.scatter(adata[adata.obs['lung'] == l].obsm['X_umap'][:, 0],
                    adata[adata.obs['lung'] == l].obsm['X_umap'][:, 1],
                    color=lung_color_code[l], label=l)
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.legend()
    plt.show()


def plt_proportion (adata):
    dict1 = dict(adata[adata.obs['sample'] == '0_Per_1_M24'].obs['clusters'].value_counts(normalize=True) * 100)
    dict2 = dict(adata[adata.obs['sample'] == '0_Per_2_F31'].obs['clusters'].value_counts(normalize=True) * 100)
    dict3 = dict(adata[adata.obs['sample'] == 'CTL_1_M63'].obs['clusters'].value_counts(normalize=True) * 100)
    dict4 = dict(adata[adata.obs['sample'] == 'CTL_2_F62'].obs['clusters'].value_counts(normalize=True) * 100)

    df = pd.DataFrame({'Vaped_Male' : dict1,
                   'Vaped_Female' : dict2,
                   'Non-Vaped_Male' : dict3,
                   'Non-Vaped_Female' : dict4}).T
    df = df.fillna(0)
    df = df[['0', '1', '2', '3',
            '4', '5', '6', '7']]
    ax = df.plot.barh(
        stacked = True,
        title = 'Proportion of Clusters')
    ax.legend(bbox_to_anchor=(1,1), loc="upper left")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_xaxis().set_visible(False)


def plt_volcano_clt (adata, p_value=0.05):
  cats = adata.obs['clusters'].cat.categories
  full_de_res = adata.uns['cluster_markers'].copy()
  full_de_res['log10_pscore'] = -np.log10(full_de_res['proba_not_de'])
  full_de_res = full_de_res.join(adata.var, how='inner')

  for c in cats:
    de_per_cluster = full_de_res.loc[full_de_res['group1'] == c]
    plt.figure(figsize = (6,6))

    ax = sns.scatterplot(data=de_per_cluster, x='lfc_mean', y='log10_pscore',
                         color=plt.cm.tab10(int(c)),
                         edgecolor=plt.cm.tab10(int(c)))
    FDR_cutoff = -np.log10(p_value)
    ax.axhline(FDR_cutoff, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(1, zorder = 0, c = 'k', lw = 2, ls = '--')
    ax.axvline(-1, zorder = 0, c = 'k', lw = 2, ls = '--')

    for axis in ['bottom', 'left']:
        ax.spines[axis].set_linewidth(2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.tick_params(width = 0)

    plt.xticks(size = 10, weight = 'bold')
    plt.yticks(size = 10, weight = 'bold')

    plt.xlabel('$log_{2}$ fold change', size = 30)
    plt.ylabel('$-log_{10}$ FDR', size = 30)
    plt.title('cluster_' + c, fontsize=25)

    ax.xaxis.set_tick_params(labelsize=30, which='major')
    ax.yaxis.set_tick_params(labelsize=30, which='major')

    plt.show()



