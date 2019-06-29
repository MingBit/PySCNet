
library(SINCERA)

args <- commandArgs(TRUE)

create_sincera <- function(Expr, clusterInfo, Mms_TF){
  
  #build R S4 object
  Expr = read.csv(Expr, sep = ',', row.names = 1)
  clusterInfo = (read.csv(clusterInfo, sep = ','))$x
  Mms_TF = as.character(read.csv(Mms_TF, sep = ',')$x)
  
  sc <- construct(exprmatrix = Expr, samplevector = colnames(Expr))
  
  # Add cluster info
  sc@data$CLUSTER = as.factor(clusterInfo)
  sc@data$GROUP = as.factor(clusterInfo)
  sc <- normalization.zscore(sc, pergroup=FALSE)
  #build genes info dataframe
  genes = data.frame(row.names = rownames(Expr))
  genes$TF = ifelse(rownames(Expr) %in% Mms_TF, 1, 0)
  
  #add TF knowledge
  sc <- setTFs(sc, value = rownames(genes)[which(genes$TF == 1)])
  
  return (sc)
  
}


run_sincera <- function(sc, clu, select_gene, output){
  
  #Add TF and selected genes
  sc@df[[clu]]$tgs = select_gene
  sc@df[[clu]]$tfs = select_gene[select_gene %in% sc@tfs]
  
  #infer a TRN 
  sc <- drivingforce.inferTRN(sc, groups = clu)
  
  write.csv(sc@df[[clu]]$edges , output, quote = FALSE, row.names = FALSE)
  
}

sc <- create_sincera(args[1], args[2], args[3])


for(i in levels(sc@data$GROUP)){
  select_gene = as.character(read.csv(args[4], sep = ',')$x)
  run_sincera(sc, clu = i, select_gene, output = paste0("linkage_", i, "txt"))
}
