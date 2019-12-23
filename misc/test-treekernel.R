library(ape)
library(igraph)
source('treekernel.R')

set.seed(123)
t1 <- rtree(8)
t2 <- rtree(8)



x <- tree.kernel(t1,t2)
y <- tree.kernel(t1,t1)
