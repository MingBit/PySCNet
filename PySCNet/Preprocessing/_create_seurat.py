import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter


def _create_seurat(expr, design, result_file, **kwargs):

        """ seurat pipeline for QC, clustering and marker gene detection
        """
        utils = rpackages.importr('utils')
        utils.chooseCRANmirror(ind = 16)
        package_names = ('M3Drop', 'scater', 'tibble', 'ggpubr', 'dplyr', 'data.table',
                        'SingleCellExperiment', 'openxlsx', 'Seurat', 'base', 'Matrix')

        names_to_install = [x for x in package_names if not rpackages.isinstalled(x)]

        if len(names_to_install) > 0:
                utils.install_packages(StrVector(names_to_install))

        for package in package_names:
                vars()[package] = rpackages.importr(package)


        #Seurat pipeline
        r('''
          Seurat_pipe <- function(sc_df, Design){
                          Seurat_obj = CreateSeuratObject(raw.data = sc_df, normalization.method = 'LogNormalize', meta.data =
                                                          data.frame(row.names = colnames(sc_df),
                                                                     cell_type = Design[colnames(sc_df),]))


                         SC3_Sureat@scale.data <- t(scale(t(SC3_Sureat@data)))
                         SC3_Sureat <- RunPCA(object = SC3_Sureat, pc.genes = rownames(sc_df), pcs.compute = 50,
                                              do.print = TRUE, pcs.print = 1:5, genes.print = 5)

                         SC3_Sureat <- FindClusters(SC3_Sureat, reduction.type = 'pca', dims.use = 1:8, k.param =150,
                                                    n.iter =1000)
                         SC3_Sureat <- RunTSNE(object = SC3_Sureat, reduction.use = "pca", dims.use = 1:8, nthreads = 4,
                                               reduction.name = "FItSNE", reduction.key = "FItSNE_", max_iter = 2000)

                         tmp <- as.factor(as.numeric(SC3_Sureat@ident))
                         names(tmp) <- SC3_Sureat@cell.names
                         SC3_Sureat@ident <- tmp


                        }
          ''')
        Seurat_pipe = ro.globalenv['Seurat_pipe']
        with localconverter(ro.default_converter + pandas2ri.converter):
                  Seurat_Obj = Seurat_pipe(expr, design)

        return(Seurat_Obj)