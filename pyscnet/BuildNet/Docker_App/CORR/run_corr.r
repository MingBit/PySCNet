library(ppcor)
library(lsa)

args = commandArgs(TRUE)
run_corr <- function(Expr, pcor_type = c('pcor', 'spcor', 'cosine'), method, output){
  
  if(pcor_type == 'pcor'){
    method = c("pearson", "kendall", "spearman")
    tmp = pcor(t(Expr), method)
    cc = tmp$estimate
    cc[tmp$p.value > 0.05] = 0
  }
  else if(pcor_type == 'spcor'){
    method = c("pearson", "kendall", "spearman")
    tmp = spcor(t(Expr), method)
    cc = tmp$estimate
    cc[tmp$p.value > 0.05] = 0
  }
  
  else{
    cc = cosine(as.matrix(t(Expr)))
  }
  colnames(cc) = rownames(Expr)
  rownames(cc) = rownames(Expr)
  cc[lower.tri(cc, diag = TRUE)] = NA
  cc = as.data.frame(as.table(cc))
  cc = na.omit(cc)
  colnames(cc) = c('source', 'target', 'weight')
  write.table(cc, output, row.names = F, quote = F, sep = '\t', col.name = F)
} 

Expr = read.csv(args[1], sep ='\t', header = T, row.names = 1)
run_corr(Expr, pcor_type = 'pcor', method = 'pearson', output = 'links.txt')
