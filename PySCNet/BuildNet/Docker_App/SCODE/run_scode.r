install.packages('igraph')
library(igraph)

matrix = as.matrix(read.csv('MING_V9T/PhD_Pro/PySCNet/BuildNet/Docker_App/SCODE/out/A.txt', sep = '\t', header = F))
g <- graph.adjacency(matrix)

### what is the output of W???