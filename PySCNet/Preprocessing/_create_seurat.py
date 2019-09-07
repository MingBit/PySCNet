import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter


def seurat_pipeline(expr, design, **kwargs):

        """ seurat pipeline for QC, clustering and marker gene detection
        """
        utils = rpackages.importr('utils')
        utils.chooseCRANmirror(ind = 16)
        package_names = ('M3Drop', 'tibble', 'ggpubr', 'dplyr', 'data.table',
                        'SingleCellExperiment', 'openxlsx', 'Seurat', 'Rtsne', 'scater')

        names_to_install = [x for x in package_names if not rpackages.isinstalled(x)]

        if len(names_to_install) > 0:
                utils.install_packages(StrVector(names_to_install))

        for package in package_names:
                print(package)
                vars()[package] = rpackages.importr(package)


        #Seurat pipeline
        r('''
          set.seed(1)
          Seurat_pipe <- function(sc_df, Design){
                         Seurat_Obj = CreateSeuratObject(raw.data = sc_df, normalization.method = 'LogNormalize', meta.data =
                                                          data.frame(row.names = colnames(sc_df),
                                                                     cell_type = Design[colnames(sc_df),]))


                         Seurat_Obj <- FindVariableGenes(object = Seurat_Obj, mean.function = ExpMean, dispersion.function = LogVMR,
                                                         x.low.cutoff = 0.0125, x.high.cutoff = 3,
                                                         y.cutoff = 0.5, do.plot = F)
                         Seurat_Obj <- ScaleData(Seurat_Obj)
                         Seurat_Obj <- RunPCA(object = Seurat_Obj, pc.genes = rownames(sc_df), pcs.compute = 20,
                                              do.print = TRUE, pcs.print = 1:5, genes.print = 5)

                         Seurat_Obj <- FindClusters(Seurat_Obj, reduction.type = 'pca', dims.use = 1:8, k.param =150,
                                                    n.iter =1000)
                         Seurat_Obj <- RunTSNE(object = Seurat_Obj, reduction.use = "pca", dims.use = 1:8, nthreads = 4,
                                               reduction.name = "FItSNE", reduction.key = "FItSNE_", max_iter = 2000, check_duplicates = FALSE)

                         tmp <- as.factor(as.numeric(Seurat_Obj@ident))
                         names(tmp) <- Seurat_Obj@cell.names
                         Seurat_Obj@ident <- tmp

                         annotation = data.frame(row.names = Seurat_Obj@cell.names,
                                                 cell_type = Seurat_Obj@meta.data$cell_type,
                                                 seurat_cluster = Seurat_Obj@ident)
                         return(list(cluster = annotation, obj = Seurat_Obj))

                        }
          ''')

        Seurat_pipe = ro.globalenv['Seurat_pipe']
        print('start seurat pipeline')
        with localconverter(ro.default_converter + pandas2ri.converter):
                Seurat_res = Seurat_pipe(expr, design)

        return(Seurat_res)