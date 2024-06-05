import pandas as ps
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


def plt_umap (adata, obs):
    color_code = dict(zip(range(0, 10), plt.cm.tab10(range(0, 10))))
    plt.figure(figsize=(10,8))
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_ticks([])
    frame1.axes.get_yaxis().set_ticks([])
    plt.scatter(adata.obsm['X_umap'][:, 0],
                adata.obsm['X_umap'][:, 1],
                c=adata.obs[obs].astype(int).map(color_code))
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.legend()
    plt.show()


def plt_umap_sample (adata):
    plt.figure(figsize=(10,8))
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_ticks([])
    frame1.axes.get_yaxis().set_ticks([])
    sample1 = plt.scatter(adata[adata.obs['sample'] == '0_Per_1_M24'].obsm['X_umap'][:, 0],
                          adata[adata.obs['sample'] == '0_Per_1_M24'].obsm['X_umap'][:, 1],
                          color=plt.cm.tab20b(range(0,20)[12]), label='0_Per_1_M24')
    sample2 = plt.scatter(adata[adata.obs['sample'] == '0_Per_2_F31'].obsm['X_umap'][:, 0],
                          adata[adata.obs['sample'] == '0_Per_2_F31'].obsm['X_umap'][:, 1],
                          color=plt.cm.tab20b(range(0,20)[15]), label='0_Per_2_F31')
    sample3 = plt.scatter(adata[adata.obs['sample'] == 'CTL_1_M63'].obsm['X_umap'][:, 0],
                          adata[adata.obs['sample'] == 'CTL_1_M63'].obsm['X_umap'][:, 1],
                          color=plt.cm.tab20b(range(0,20)[8]), label='CTL_1_M63')
    sample4 = plt.scatter(adata[adata.obs['sample'] == 'CTL_2_F62'].obsm['X_umap'][:, 0],
                          adata[adata.obs['sample'] == 'CTL_2_F62'].obsm['X_umap'][:, 1],
                          color=plt.cm.tab20b(range(0,20)[11]), label='CTL_2_F62')
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.legend()
    plt.show()


def plt_umap_sex (adata):
    plt.figure(figsize=(10,8))
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_ticks([])
    frame1.axes.get_yaxis().set_ticks([])
    sample1 = plt.scatter(adata[adata.obs['sex'] == 'female'].obsm['X_umap'][:, 0],
                          adata[adata.obs['sex'] == 'female'].obsm['X_umap'][:, 1],
                          color=plt.cm.Accent(range(0,20)[5]), label='Female')
    sample2 = plt.scatter(adata[adata.obs['sex'] == 'male'].obsm['X_umap'][:, 0],
                          adata[adata.obs['sex'] == 'male'].obsm['X_umap'][:, 1],
                          color=plt.cm.tab20b(range(0,20)[6]), label='Male')
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.legend()
    plt.show()


def plt_umap_condition (adata):
    plt.figure(figsize=(10,8))
    frame1 = plt.gca()
    frame1.axes.get_xaxis().set_ticks([])
    frame1.axes.get_yaxis().set_ticks([])
    sample1 = plt.scatter(adata[adata.obs['lung'] == 'Z_Per'].obsm['X_umap'][:, 0],
                          adata[adata.obs['lung'] == 'Z_Per'].obsm['X_umap'][:, 1],
                          color=plt.cm.tab20b(range(0,20)[17]), label='Z_Per')
    sample2 = plt.scatter(adata[adata.obs['lung'] == 'CTL'].obsm['X_umap'][:, 0],
                          adata[adata.obs['lung'] == 'CTL'].obsm['X_umap'][:, 1],
                          color=plt.cm.tab20c(range(0,20)[9]), label='CTL')
    plt.xlabel('UMAP_1')
    plt.ylabel('UMAP_2')
    plt.legend()
    plt.show()
