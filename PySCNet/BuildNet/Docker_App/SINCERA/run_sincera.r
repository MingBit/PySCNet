
library(SINCERA)

args <- commandArgs(TRUE)

create_sincera <- function(Expr, clusterInfo, Mms_TF){
  
  #build R S4 object
  # Expr = read.csv(Expr, sep = '\t', row.names = 1)
  # clusterInfo = (read.csv(clusterInfo, sep = '\t'))$x
  # Mms_TF = as.character(read.csv(Mms_TF, sep = '\t')$x)
  
  sc <- construct(exprmatrix = Expr, samplevector = c(rep('sample_1', ncol(Expr))))
  
  # Add cluster info
  sc@data$CLUSTER = as.factor(clusterInfo)
  sc@data$GROUP = as.factor(clusterInfo)
  sc <- normalization.zscore(sc, pergroup=FALSE)
  #build genes info dataframe
  genes = data.frame(row.names = rownames(Expr))
  genes$TF = ifelse(toupper(rownames(Expr)) %in% toupper(Mms_TF), 1, 0)
  
  #add TF knowledge
  sc <- setTFs(sc, value = rownames(genes)[which(genes$TF == 1)])
  
  return (sc)
  
}


run_sincera <- function(sc, clu, output){
  
  #Add TF and selected genes
  select_gene = sc@genes.forclustering
  sc@df[[clu]]$tgs = select_gene
  sc@df[[clu]]$tfs = select_gene[select_gene %in% sc@tfs]
  
  #infer a TRN 
  sc <- drivingforce.inferTRN(sc, groups = clu)
  
  write.csv(sc@df[[clu]]$edges , output, quote = FALSE, row.names = FALSE)
  
}


sc <- create_sincera(args[1], args[2], args[3])

Expr = read.csv(file = '~/MING_V9T/PhD_Pro/Test/Expr.txt', sep = '\t', row.names = 1, header = T)
ClusterInfo = data.frame(row.names = colnames(Expr))
ClusterInfo['clusterid'] = c(rep('1', 241))
Mms_TF = read.csv('MING_V9T/PhD_Pro/PySCNet/Mus_TF_and_TFcofactor/Mus_musculus_TF.txt', sep = '\t')

sc <- create_sincera(Expr, ClusterInfo$clusterid, Mms_TF$Symbol)

for(i in levels(sc@data$GROUP)){
  
  run_sincera(sc, clu = i, select_gene, output = paste0("links_", i, "txt"))
}
