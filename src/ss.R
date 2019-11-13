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
treename <- system('ls trees/tree-mc-*.nwk', intern = TRUE)
n.mc = length(treename)
x <- list()
for(i in 1:n.mc){
    x[[i]] <- ape::read.tree(file = treename[i])
    if(i%%10==0) message(paste('Read tree',i,'/',n.mc))
}

par(mfrow=c(2,2))
plot(tstar, main='True Phylogeny')
for(j in 1:3) plot(x[[j]], main=paste('FastTree',j))


# ---- Rerooting ----

# Reroot the trees 
# (a rooted tree is needed for summary stats)

# TODO: 
# Check if that makes sense. I'm not sure at all this is correct!
# https://phylobotanist.blogspot.com/2015/01/how-to-root-phylogenetic-tree-outgroup.html
# https://cabbagesofdoom.blogspot.com/2012/06/how-to-root-phylogenetic-tree.html

# Generate random dummy dates for tips:
n.tips <- length(x[[1]]$tip.label)
dummy.tip.dates <- runif(n = n.tips, 
                         min = 0, 
                         max = 1)

for(i in 1:n.mc){
  if(i%%10==0 || i==1) print(paste('Re-rooting tree #',i,'/',n.mc))
  x[[i]] <- ape::rtt(t = x[[i]], 
                     tip.dates = dummy.tip.dates,
                     objective = 'rms') #correlation rms rsquared
}


# ---- Summary stats ----

# Calculate statistics for each MC tree:
ss.list <- list()
for(i in 1:n.mc){
  ss.list[[i]]  <- c(
    treeCentrality::computeBasicStats(x[[i]]),
    treeCentrality::computeSpectralStats(x[[i]]),
    'diameter' = treeCentrality::computeDiameter(x[[i]])
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


plot(tstar, type = 'phylo')
ape::cherry(tstar)
ape::cherry(x[[1]])

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
    
# ---- Plots ----

# plot a random selection of 9 trees:
pdf('plot-trees.pdf', width = 8, height = 8)
n.plot <- min(n.mc, 9)
par(mfrow=c(3,3))
for(i in seq_along(x))
  plot.phylo(x[[i]], type = 'phyl', 
       main=paste('Rerooted FastTree',i))
dev.off()


# Plot CV of summary statistics
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
