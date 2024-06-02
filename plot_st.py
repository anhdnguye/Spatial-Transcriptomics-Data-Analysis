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


def total_plot(adata, lstGenes):
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
    ax.colorbar(label='Expression')
