library(ggplot2)
library(ape)


source('dist-fcts.R')

args <- commandArgs(trailingOnly = TRUE)

# args <- 3

software <- 'raxml'   # 'raxml' or 'fasttree'

message(paste('Inference assessment for prmset',args[1],'...'))

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

# ---- _benchmark ----

d.rf.star <- sapply(x, dist.RF, tstaru)
d.sh.star <- sapply(x, dist.shared, tree.ref = tstar)

# ---- _b/w inferred ----

d.rf <- numeric(n.mc*(n.mc-1)/2) 
d.sh <- numeric(n.mc*(n.mc-1)/2) 
k    <- 1

for(i in seq_along(x)){
  if(i%%10==0) message(paste('Calculating distances',i,
                             '/',n.mc,'...'))
  
  # Distance b/w the `x`s:
  for(j in 1:i){
    if(j != i) {
      d.rf[k] <- dist.RF(x[[i]], x[[j]])
      d.sh[k] <- dist.shared(x[[i]], x[[j]])
      k <- k+1
    }
  }
}


# ---- Save ----

dist.list <- list(d.rf.star = d.rf.star, 
                  d.sh.star = d.sh.star, 
                  d.rf      = d.rf,
                  d.sh      = d.sh,
                  prmset    = args[1])
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
hist(d.rf.star, col = 'lightgrey')
hist(d.rf, col = 'lightgrey')
dev.off()

message(paste('Assessment done for parameter set', args[1]))


