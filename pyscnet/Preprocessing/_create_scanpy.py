import scanpy.api as sc
import anndata
import pandas as pd


def scanpy_pipeline(expr):
    """ scanpy pipeline for QC, clustering, marker gene identification
        """
    adata = anndata.AnnData(expr.T)
    adata.var = pd.DataFrame(index=expr.index)

    # filtering
    sc.pp.filter_cells(adata, min_genes=500)
    sc.pp.filter_genes(adata, min_counts=50)

    #        mito_genes = adata.var_names.str.startswith('MT-')
    #        adata.obs['percent_mito'] = np.sum(
    #                        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    #
    #        adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    #        adata = adata[adata.obs['n_genes'] < 3000, :]
    #        adata = adata[adata.obs['percent_mito'] < 0.0ls 5, :]

    # normalize
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)

    # gene selection
    sc.pp.highly_variable_genes(adata)

    adata_filter = adata[:, adata.var['highly_variable']]

    sc.pp.scale(adata_filter, max_value=10, copy=True)

    # pca and clustering
    sc.tl.pca(adata_filter)
    sc.pp.neighbors(adata_filter, n_neighbors=200, n_pcs=20)
    sc.tl.louvain(adata_filter)
    sc.tl.tsne(adata_filter)
    sc.tl.umap(adata_filter)

    # find marker gene
    sc.tl.rank_genes_groups(adata_filter, 'louvain', method='wilcoxon')
    #        adata.write(result_file)

    return adata_filter

