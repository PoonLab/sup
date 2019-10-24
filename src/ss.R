library(ggplot2)
library(ape)
tc <- require(treeCentrality)
if(!tc){
devtools::install_github('Leonardini/treeCentrality')  
}

message("\n\nReconstruction analysis started...")

# ---- Data ----

# Retrieve all trees calculated upstream:
treename <- system('ls trees/tree-mc-*.out', intern = TRUE)
n = length(treename)
x <- list()
for(i in 1:n){
    x[[i]] <- ape::read.tree(file = treename[i])
    if(i%%10==0) message(paste('Read tree',i,'/',n))
}

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
# dummy.tip.dates <- rnorm(n.tips, mean = 100, sd = 11)

for(i in 1:n){
  if(i%%10==0) print(paste('Re-rooting tree #',i,'/',n))
  x[[i]] <- ape::rtt(t = x[[i]], 
                     tip.dates = dummy.tip.dates,
                     objective = 'rms') #correlation rms rsquared
}

plot(z, type = 'phyl')
plot(x[[1]], type = 'phyl')

# ---- Summary stats ----

# Calculate statistics for each MC tree:
ss.list <- list()
for(i in 1:n){
  ss.list[[i]]  <- c(
    treeCentrality::computeBasicStats(x[[i]]),
    treeCentrality::computeSpectralStats(x[[i]]),
    'diameter' = treeCentrality::computeDiameter(x[[i]])
  )
  
}
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
par(mfrow=c(3,3))
for(i in sort(sample(1:n, 9)))
    plot.phylo(x[[i]], 
               main = paste0('tree-',i), 
               type = 'phylo')
dev.off()


# plot summary stats:

extract_ss <- function(ss.list, varname) {
  y <- sapply(ss.list, '[[', varname)
  m <- mean(y)
  stdv <- sd(y)
  cv <- stdv / m
  if(m==0) cv <- NA
  return(cv)
}

nm  <- names(ss.list[[1]])
nm2 <- ceiling(sqrt(length(nm)))

par(mfrow=c(nm2,nm2))
cv <- numeric(length(nm))
for(i in 1:length(nm)){
  cv[i] <- extract_ss(ss.list, nm[i])
}
df <- data.frame(nm = nm, cv=cv)



y <- sapply(ss.list, '[[', 'cherries')
hist(y, col='grey')
plot(x=1:length(y), y)

pdf('plot-ss.pdf', width = 12, height = 10)
g <- ggplot(df, aes(x=nm,y=cv))+
  geom_bar(stat='identity')+
  coord_flip()+
  ggtitle('coeff. of variation')
plot(g)
dev.off()
