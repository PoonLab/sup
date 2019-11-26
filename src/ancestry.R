library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())

library(ape)
library(phytools)
library(phylobase)

nwks <- system('ls ./trees/RAxML_bestTree*.out', intern = TRUE)
n <- length(nwks)




descendants_list <- function(t, do.root=F, do.plot = F) {
    
    if(do.plot) plot(t)
    
    phy <- phylo4(t)
    
    if(do.root){
        tr <- phytools::midpoint.root(t) %>%
            phylo4()
        if(do.plot) plot(tr)
        phy <- tr
    }
    
    #print(phy)
    internal.nodes <- nodeId(phy, 'internal')
    
    a <- list() ; k=1
    for (i in internal.nodes) {  # i=12
        tmp <- phylobase::descendants(phy, i)
        a[[k]] <- sort(names(tmp))
        k <- k+1
    }
    return(a)
}


# Proportion shared against the first sim (mc=1, prmset=1)
i = 1

mb <- numeric(length = n)
for(j in 1:n){
    
    ti <- ape::read.tree(nwks[i])
    tj <- ape::read.tree(nwks[j])
    
    xi <- descendants_list(ti)
    xj <- descendants_list(tj)
    
    b <- vector(length = length(xi))
    for(k in seq_along(xi)){
        b[k] <- xi[k] %in% xj
    }
    mb[j] <- mean(b)
}

# Need to merge them by prmset
# quick and dirty:
n.mc = 30 # TODO: read from file
n.prmset <- n / n.mc
prmset <- rep(1:n.prmset, each=n.mc)

df <- data.frame(prmset, mb)

g <- ggplot(df) +
    geom_violin(aes(x=factor(prmset), y=mb))
plot(g)

dfs <- df %>%
    group_by(prmset) %>%
    summarize(m = mean(mb), 
              md = median(mb),
              s = sd(mb),
              qlo = quantile(mb, probs = 0.05),
              qhi = quantile(mb, probs = 0.95))
dfs

gs <- dfs %>%
    ggplot(aes(x=prmset, y=m))+
    geom_pointrange(aes(ymin=qlo, ymax=qhi))+
    geom_line()+
    geom_point(aes(y=md), shape=4, size=3)+
    coord_cartesian(ylim=c(0,1))+
    ylab('mean ancestry shared')


gsv <- dfs %>%
    ggplot(aes(x=prmset, y=s))+
    geom_point()+
    geom_line()

plot(gsv)
plot(gs)


# Against each other:
M <- matrix(nrow = n, ncol=n)
for(i in 1:n){
    print(paste(i,'/',n))
    for(j in 1:i){
        ti <- ape::read.tree(nwks[i])
        tj <- ape::read.tree(nwks[j])
        
        xi <- descendants_list(ti)
        xj <- descendants_list(tj)
        
        b <- vector(length = length(xi))
        for(k in seq_along(xi)){
            b[k] <- xi[k] %in% xj
        }
        M[i,j] <- mean(b)
    }
}
image(M, axes=T, col = terrain.colors(12))





