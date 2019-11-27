library(ggplot2)
library(ape)


message('Inference assessment...')

args <- commandArgs(trailingOnly = TRUE)

# args <- 1

software <- 'raxml'   # 'raxml' or 'fasttree'

# ---- Data ----

# Retrieve the "true" benchmark tree:
tstar  <- ape::read.tree('trees/sim.nwk')
tstaru <- ape::unroot(tstar)

# Retrieve tree infered from "certain" sequence:
if(software=='raxml')    
  t.certn <- ape::read.tree('trees/RAxML_bestTree.tree-raxml-certain.out')


# Retrieve all trees calculated upstream:
if(software=='fasttree') fname <- 'tree-mc-'
if(software=='raxml')    fname <- paste0('RAxML_bestTree.tree-raxml-prm-',
                                         args[1],'-mc-')

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


# ---- Distances ----

# Robinson-Foulds ("edit") distance:

d.star <- numeric(length = n.mc)
d.star <- numeric(length = n.mc)
d      <- numeric(n.mc*(n.mc-1)/2) ; k=1

rf.normalize <- TRUE # if the RF distance is normalized.

for(i in seq_along(x)){
  
  if(i%%10==0) message(paste('Distances',i,'/',n.mc))
  
  # Distance normalization:
  n.tips <- length(x[[i]]$tip.label)
  tmp <- 2 * (n.tips - 3)
  rf.norm <- ifelse(rf.normalize, tmp, 1.0)
  
  # Distance b/w `x` and `tstar`:
  d.star[i] <- ape::dist.topo(tstaru, x[[i]], 
                              method='PH85') / rf.norm
  
  # Distance b/w the `x`s:
  for(j in 1:i){
    if(j != i) {
      d[k] <- ape::dist.topo(x[[i]], x[[j]], 
                             method='PH85') / rf.norm
      k <- k+1
    }
  }
}

# ---- Save ----

dist.list <- list(d.star = d.star, d=d, prmset = args[1])
save(list = 'dist.list', 
     file = paste0('treedist-',args[1],'.RData'))

# ---- Plots ----

message('Plotting assessment...')

plot_tree <- function(p, title='') {
  plot(p, main=title)
  edgelabels(round(p$edge.length,3), cex=0.65)
}

# plot some trees:

pdf(paste0('plot-trees-prm-',args[1],'.pdf'), 
    width = 12, height = 10)
n.plot <- min(n.mc, 7)
par(mfrow=c(3,3))

plot_tree(tstaru, 'True Tree')
plot_tree(t.certn, 'Tree from certain seq.')
for(i in 1:n.plot){
  plot_tree(x[[i]], paste(software,i))
}
dev.off()

pdf(paste0('plot-ss-prm-',args[1],'.pdf'))
par(mfrow=c(1,2))
hist(d.star, col = 'lightgrey')
hist(d, col = 'lightgrey')
dev.off()

message(paste('Assessment done for parameter set', args[1]))

# ----- OLD STUFF -----

# tc <- require(treeCentrality)
# if(!tc){
#   devtools::install_github('Leonardini/treeCentrality')  
# }
# 
# message("\n\nReconstruction analysis started...")

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

if(0){
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
  
}
