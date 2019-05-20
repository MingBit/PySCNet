from __future__ import absolute_import
import sys
#sys.path.append('/home/mwu/MING_V9T/PhD_Pro/SCNetEnrich/BuildNet')


from collections import OrderedDict
import rpy2.robjects as ro
from BuildNet.rlib.imports import base, wgcna, rsnippets, stats, dynamicTreeCut, grdevices

#from r.manager import RManager


class PyWGCNA():
    '''
    WGCNA python object to run the general pipleline for WGCNA. Main functions are
    adjacency(), TOM_Dist(), TOM_Simlarity, picksoftthershold, Generate_GeneTree,
    cutreeDynamic(), module_EigenGene() and plots fucntions.
    '''
    def __init__(self, data, params = None, debug=False):
        self.Expr = data.T
        self.debug = debug
        self.param = params
        self.adjacencyMatrix = None
        self.TOM = None
        self.dissimilarityMatrix = None
        self.geneTree = None
        self.moduleColors = None
        self.softthreshold = None
        self.dynamicModule = None
        self.moduleEigenGene = None
        return None

    def collect_garbage(self):
        wgcna().collectGarbage()

    def pickSoftThreshold(self, method = 'default', Nettype = 'unsigned',
                      powerVector = list(range(1,10,2))):
        """scale free topology for softthreshold
         data,
          dataIsExpr = TRUE,
          weights = NULL,
          RsquaredCut = 0.85,
          powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
          removeFirst = FALSE, nBreaks = 10, blockSize = NULL,
          corFnc = cor, corOptions = list(use = 'p'),
          networkType = "unsigned",
          moreNetworkConcepts = FALSE,
          gcInterval = NULL,
          verbose = 0, indent = 0)
        """
        subPara = {}
        subPara['powerVector'] = powerVector
        subPara['networkType'] = Nettype

        if(method == 'default'):
            sft = wgcna().pickSoftThreshold(self.Expr, **subPara)
        elif(method == 'similarity'):
            self.TOM = wgcna().TOMsimilarityFromExpr(datExpr = self.Expr)
            sft = wgcna().pickSoftThreshold.fromSimilarity(self.TOM, **subPara)
        else: print('only default and similarity can be used for method parameter')

    #        index <- which(diff(sign(diff(-sign(sft$fitIndices[,3])*sft$fitIndices[,2])))==-2)+1
    #        estimateBeta <- sft$fitIndices[index,1][1]

        self.softthreshold = sft
        self.collect_garbage()



    def adjacency(self, Nettype = 'unsigned', power = 6):
        """build adjacency matrix
            datExpr,
          selectCols = NULL,
          type = "unsigned",
          power = if (type=="distance") 1 else 6,
          corFnc = "cor", corOptions = list(use = "p"),
          weights = NULL,
          distFnc = "dist", distOptions = "method = 'euclidean'",
          weightArgNames = c("weights.x", "weights.y"))
        """
        adjParams = {}
        adjParams['power'] = power
        adjParams['corFnc'] = 'cor'
        adjParams['corOptions'] = "use='p'"
        adjParams['type'] = Nettype
        adjParams['datExpr'] = self.Expr

        self.adjacencyMatrix = wgcna().adjacency(**adjParams)
        self.collect_garbage()



    def Tom_dist_similarity(self, method = 'dist'):
        """calculate topological overlap matrix distance/similarity
              datExpr,
              weights = NULL,
              corType = "pearson",
              networkType = "unsigned",
              power = 6,
              TOMType = "signed",
              TOMDenom = "min",
              maxPOutliers = 1,
              quickCor = 0,
              pearsonFallback = "individual",
              cosineCorrelation = FALSE,
              replaceMissingAdjacencies = FALSE,
              suppressTOMForZeroAdjacencies = FALSE,
              suppressNegativeTOM = FALSE,
              useInternalMatrixAlgebra = FALSE,
              nThreads = 0,
              verbose = 1, indent = 0)
        """
        if(method == 'dist'):
            self.TOM = wgcna().TOMdist(self.adjacencyMatrix)
            self.dissimilarityMatrix = rsnippets.dissMatrix(self.adjacencyMatrix)
            self.collect_garbage()
        elif(method == 'similarity'):
            TomParam = {}
            self.TOM = wgcna().TOMsimilarityFromExpr(datExpr = self.Expr, **TomParam)
            self.dissimilarityMatrix = rsnippets.dissMatrix(self.TOM)
            self.dissimilarityMatrix = rsnippets.powerWeightMatrix(self.TOM, 7)
            self.collect_garbage()
        else: print('only dist and similarity can be used for method')


    def generate_geneTree(self):
        """
        """
#        import pdb
#        pdb.set_trace()
#        distTom = 1 - np.matrix(self.TOM)
        distMatrix = stats().as_dist(self.dissimilarityMatrix)

        self.geneTree = stats().hclust(distMatrix, method="average")
        self.collect_garbage()


    def cutreeDynamic(self, minModuleSize = 5, maxTreeHeight=5):

        """Detect clusters in a hierarchical dendrogram using a variable cut height approach.

        """
        subpara = {}
        subpara['maxTreeHeight'] = maxTreeHeight
        subpara['minModuleSize'] = minModuleSize
        self.dynamicModule = dynamicTreeCut().cutreeDynamicTree(self.geneTree, **subpara)
        self.collect_garbage()



    def set_module_colors(self, modules):
        '''
        set module colors from dict of module properties
        expects dict to have a 'color' sub-dict
        '''
        self.moduleColors = OrderedDict((module, values['color']) \
                                        for module, values in modules.items())


    def module_eigen_gene(self):
        """detect eigen gene for each module. need to debug
        expr,
         colors,
         impute = TRUE,
         nPC = 1,
         align = "along average",
         excludeGrey = FALSE,
         grey = if (is.numeric(colors)) 0 else "grey",
         subHubs = TRUE,
         trapErrors = FALSE,
         returnValidOnly = trapErrors,
         softPower = 6,
         scale = TRUE,
         verbose = 0, indent = 0
        """
        subPara = {}
        subPara['colors'] = wgcna().labels2colors(self.dynamicModule)
        self.module_eigen_gene = wgcna().moduleEigengenes(self.Expr, **subPara)


    def plot_eigengene_network(self, filepath = None, width = 512, height = 512):
        '''
        wrapper for plotEigengeneNetworks function
        plots an eigengene network
        '''
#        import pdb
#        pdb.set_trace()
        params = {}
#        params['multiME'] = base().as_data_frame(np.array(self.module_eigen_gene)[2])
        params['multiME'] = base().as_data_frame(self.Expr)
        params['setLabels'] = ''
        params['marDendro'] = ro.IntVector([0, 4, 1, 2])
        params['marHeatmap'] = ro.IntVector([3, 4, 1, 2])
        params['cex.lab'] = 0.8
        params['xLabelsAngle'] = 90
        params['colorLabels'] = False
        params['signed'] = True

        grdevices().png(filepath, width=width, height=height)
        wgcna().plotEigengeneNetworks(**params)
        grdevices().dev_off()


    def plotDendroHeatmap(self, filepath = None, width = 512, height = 512):

        """plot Dendrogram heatmap with color annotations of objects
          dendro,
          colors,
          groupLabels = NULL,
          rowText = NULL,
          rowTextAlignment = c("left", "center", "right"),
          rowTextIgnore = NULL,
          textPositions = NULL,
          setLayout = TRUE,
          autoColorHeight = TRUE,
          colorHeight = 0.2,
          colorHeightBase = 0.2,
          colorHeightMax = 0.6,
          rowWidths = NULL,
          dendroLabels = NULL,
          addGuide = FALSE, guideAll = FALSE,
          guideCount = 50, guideHang = 0.2,
          addTextGuide = FALSE,
          cex.colorLabels = 0.8, cex.dendroLabels = 0.9,
          cex.rowText = 0.8,
          marAll = c(1, 5, 3, 1), saveMar = TRUE,
          abHeight = NULL, abCol = "red"
        """

        subpara = {}
        subpara['dendro'] = self.geneTree
        subpara['colors'] = wgcna().labels2colors(self.dynamicModule)
        grdevices().png(filepath, width=width, height=height)
        wgcna().plotDendroAndColors(**subpara)
        grdevices().dev_off()











