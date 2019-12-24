library(ape)
library(igraph)
source('treekernel.R')

set.seed(123)
t1 <- rtree(8)
t2 <- rcoal(13)

unroot(t1)
unroot(t2)

x  <- tree.kernel(t1,t2, normalize = TRUE)

#y1 <- tree.kernel(t1,t1, normalize = TRUE)
