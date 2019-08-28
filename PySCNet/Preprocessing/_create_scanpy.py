import scanpy.api as sc
import numpy as np
import anndata






def _create_scanpy(expr, feature, **kwargs):
        """ scanpy pipeline for QC, clustering, marker gene identification
        """
        adata = anndata.AnnData(expr.T)
        adata.var = feature

        #filtering
        sc.pp.filter_cells(adata, **kwargs)
        sc.pp.filter_genes(adata, **kwargs)

        mito_genes = adata.var_names.str.startswith('MT-')
        adata.obs['percent_mito'] = np.sum(
                        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1

        adata.obs['n_counts'] = adata.X.sum(axis=1).A1
        adata = adata[adata.obs['n_genes'] < 3000, :]
        adata = adata[adata.obs['percent_mito'] < 0.05, :]

        #normalize
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)

        #gene selection
        sc.pp.highly_variable_genes(adata,
                              min_mean = 0.025, max_mean = 2, min_disp = 0.4)

        adata_filter = adata [:,adata.var['highly_variable']]
        sc.pp.scale(adata_filter, max_value=10)

        #pca and clustering
        sc.tl.pca(adata)
        sc.pp.neighbors(adata_filter, n_neighbors=50, n_pcs=20)
        sc.tl.louvain(adata_filter)
        sc.tl.tsne(adata_filter)
        sc.tl.umap(adata_filter)

        #find marker gene
        sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
#        adata.write(result_file)

        return(adata)
