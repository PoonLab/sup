library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)

library(ape)
library(phytools)
library(phylobase)

# ---- FUNCTIONS ----

#' Get all internal nodes and their descendants labels.
#' @param t phylo tree.
#' @param do.root Boolean. Root the tree if TRUE. Default = TRUE.
#' @param do.plot Boolean. Plot trees. Default = TRUE.
descendants_list <- function(t, do.root=F, do.plot = F) {
    
    if(do.plot) plot(t)
    phy <- phylobase::phylo4(t)
    
    if(do.root){
        tr <- phytools::midpoint.root(t) %>%
            phylobase::phylo4()
        if(do.plot) plot(tr)
        phy <- tr
    }
    
    #print(phy)
    internal.nodes <- nodeId(phy, 'internal')
    
    res <- list() ; k=1
    for (i in internal.nodes) {  # i=12
        tmp <- phylobase::descendants(phy, i)
        res[[k]] <- sort(names(tmp))
        k <- k+1
    }
    return(res)
}


# tree     <- ape::read.tree(nwks[149])

#' Calculate the proportion of nodes in a tree that have
#' the same descendance as a reference tree. Descendance
#' is calculated from ordered tips labels.
shared_descendance <- function(tree, tree.ref ) {
    d.ref <- descendants_list(tree.ref)
    d     <- descendants_list(tree)
    b     <- d %in% d.ref
    return(mean(b))
}


# ---- RUN ----

nwk.ref <- system('ls ./trees/RAxML_bestTree*certain*.out', intern = TRUE)
nwks    <- system('ls ./trees/RAxML_bestTree*prm*mc*.out', intern = TRUE)
n       <- length(nwks)


trees    <- lapply(nwks, ape::read.tree)
tree.ref <- ape::read.tree(nwk.ref)
sh.dsc   <- sapply(trees, shared_descendance, tree.ref = tree.ref)


# Need to merge them by prmset
# quick and dirty:
n.prmset <- as.numeric(system('wc -l < prm-btshp.csv', intern=TRUE))
n.mc     <- n / n.prmset
prmset   <- rep(1:n.prmset, each=n.mc)
df       <- data.frame(prmset, sh.dsc)



dfs <- df %>%
    group_by(prmset) %>%
    summarize(m = mean(sh.dsc), 
              md = median(sh.dsc),
              s = sd(sh.dsc),
              qlo = quantile(sh.dsc, probs = 0.05),
              qhi = quantile(sh.dsc, probs = 0.95))
dfs

g <- ggplot(df) +
    geom_violin(aes(x=factor(prmset), y=sh.dsc))


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

pdf('plot-analysis-ancestry.pdf')
grid.arrange(gs, gsv, ncol=1)
dev.off()


