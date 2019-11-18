library(ggplot2)
library(ape)
tc <- require(treeCentrality)
if(!tc){
  devtools::install_github('Leonardini/treeCentrality')  
}

message("\n\nReconstruction analysis started...")

# ---- Data ----

# Retrieve the "true" benchmark tree:
tstar <- ape::read.tree('trees/sim.nwk')


# Retrieve all trees calculated upstream:

software <- 'raxml'   # raxml  fasttree

if(software=='fasttree') fname <- 'tree-mc-'
if(software=='raxml')    fname <- 'RAxML_bestTree.tree-raxml-mc-'

treename <- system(paste0('ls trees/',
                          fname,
                          '*.out'), 
                   intern = TRUE)
n.mc = length(treename)
x <- list()
for(i in 1:n.mc){
  x[[i]] <- ape::read.tree(file = treename[i])
  if(i%%10==0) message(paste('Read tree',i,'/',n.mc))
}


ape::dist.topo(x[[1]], x[[2]] , method='PH85')
ape::dist.topo(tstar, x[[2]] , method='PH85')


plot.phylo(tstar, use.edge.length = T, main = 'True')
edgelabels(round(tstar$edge.length,3), 
           cex=0.6)
is.rooted(x[[1]])
plot.phylo(x[[1]], use.edge.length = T)
edgelabels(round(x[[1]]$edge.length,3), 
           cex=0.6)



for(i in seq_along(x)){
  print(i)
  print(ape::dist.topo(tstar, x[[i]], method='PH85'))
}

# ---- Rerooting ----

# Reroot the trees 
# (a rooted tree is needed for treeCentrality summary stats)
# TODO: 
# Check if that makes sense. I'm not sure at all this is correct!
# https://phylobotanist.blogspot.com/2015/01/how-to-root-phylogenetic-tree-outgroup.html
# https://cabbagesofdoom.blogspot.com/2012/06/how-to-root-phylogenetic-tree.html

# Generate random dummy dates for tips:
n.tips <- length(x[[1]]$tip.label)
dummy.tip.dates <- runif(n = n.tips, 
                         min = 0, 
                         max = 1)

xrr <- list()
for(i in 1:n.mc){
  if(i%%10==0 || i==1) print(paste('Re-rooting tree #',i,'/',n.mc))
  xrr[[i]] <- ape::rtt(t = x[[i]], 
                       tip.dates = dummy.tip.dates,
                       objective = 'rms') #correlation rms rsquared
}


# ---- Summary stats ----

# Calculate statistics for each MC tree:
ss.list <- list()
for(i in 1:n.mc){
  ss.list[[i]]  <- c(
    treeCentrality::computeBasicStats(xrr[[i]]),
    treeCentrality::computeSpectralStats(xrr[[i]]),
    'diameter' = treeCentrality::computeDiameter(xrr[[i]])
  )
  
}


#' Return the coefficient of variation (cv)
#' of all summary statistics calculated, 
#' across all MC iterations.
cv_ss <- function(ss.list, varname) {
  y <- sapply(ss.list, '[[', varname)
  m <- mean(y)
  stdv <- sd(y)
  cv <- stdv / m
  if(m==0) cv <- NA
  return(cv)
}

# ---- Plots ----

# plot some trees:
pdf('plot-trees.pdf', width = 12, height = 12)
n.plot <- min(n.mc, 8)
par(mfrow=c(3,3))
plot(tstar, main='True Phylogeny')
for(i in 1:n.plot){
  plot(x[[i]], main=paste('Raw FastTree',i))
  plot.phylo(xrr[[i]], type = 'phyl', 
             main=paste('Rerooted FastTree',i))
}
dev.off()


# Plot CV and histogram of summary statistics:
nm  <- names(ss.list[[1]])
nm2 <- ceiling(sqrt(length(nm)))

cv <- numeric(length(nm))
for(i in 1:length(nm)){
  cv[i] <- cv_ss(ss.list, nm[i])
}
df <- data.frame(nm = nm, cv=cv)

pdf('plot-ss.pdf', width = 12, height = 10)
g <- ggplot(df, aes(x=nm,y=cv))+
  geom_bar(stat='identity')+
  coord_flip()+
  ggtitle('coeff. of variation')
plot(g)

par(mfrow=c(nm2,nm2))
for(i in seq_along(nm)){
  y <- sapply(ss.list, '[[', nm[i])
  hist(y, breaks = 20,
       col='grey',
       main =nm[i],
       yaxt='n', xlab='', ylab='')  
}
dev.off()


# ----- OLD STUFF
# treeCentrality::computeLMStats(x[[1]])
# treeCentrality::computeNetworkStats(x[[1]])
# 
# treeCentrality::computeBetweenness(x[[1]])
# treeCentrality::computeCloseness(x[[1]])
# treeCentrality::computeDiameter(x[[1]])
# 
# treeCentrality::computeWienerIndex(x[[1]])
# treeCentrality::computeEigenvector(x[[1]])
# x[[1]]$node.label
# 
#   treeCentrality::rankDiscriminatoryStats(tList1 = x[[1:2]],
#                                           tList2 = x[[3:4]],
#                                           basic = TRUE, 
#                                           network = F, 
#                                           spectral = F)  