###
###  DEFINE THE DISTANCE FUNCTIONS
### 

suppressPackageStartupMessages({library(ape);
library(phytools);
library(phylobase)})
source('treekernel.R')

#' Robinson-Foulds distance. 
#' @param tree1 `ape::phylo` object. First tree. 
#' @param tree2 `ape::phylo` object. Second tree. 
#' @param normalize Boolean. Normalize distance by 2(n-3) where n is the maximum number of tips between tree1 and tree2. Default=TRUE.
dist.RF <- function(tree1, tree2, normalize=TRUE) {
    
    # Distance normalization:
    rf.norm <- 1.0
    if(normalize){
        n.tips1 <- length(tree1$tip.label)
        n.tips2 <- length(tree2$tip.label)
        n.tips  <- max(n.tips1, n.tips2)
        rf.norm <- 2 * (n.tips - 3)
    }
    
    # Distance b/w tree1 and tree2:
    res <- ape::dist.topo(tree1, tree2, method='PH85') / rf.norm
    
    return(res)
}


#' Get all internal nodes and their descendants labels.
#' @param t phylo tree.
#' @param do.root Boolean. Root the tree if TRUE. Default = TRUE.
#' @param do.plot Boolean. Plot trees. Default = TRUE.
.descendants_list <- function(t, do.root=F, do.plot = F) {
    
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


#' Calculate the proportion of nodes in a tree that 
#' DO NOT have the same descendance as a reference tree. 
#' Descendance is calculated from ordered tips labels.
#' @param tree `phylo` tree.
#' @param tree.ref `phylo` reference tree.
dist.shared <- function(tree, tree.ref ) {
    d.ref <- .descendants_list(tree.ref)
    d     <- .descendants_list(tree)
    b     <- d %in% d.ref
    return(1.0 - mean(b))
}

#' @param tree1 `ape::phylo` object. First tree. 
#' @param tree2 `ape::phylo` object. Second tree. 
#' @param normalize Boolean. Normalize kernel distance.
dist.kernel <- function(tree1, tree2, 
                        lambda = 0.5, 
                        sigma = 1.0, 
                        rho = TRUE,
                        normalize = TRUE) {
    chk <- try(d  <- tree.kernel(tree1, 
                                 tree2, 
                                 lambda = lambda, 
                                 sigma = sigma, 
                                 rho = rho,
                                 normalize = normalize), 
               silent = TRUE)
    if(class(chk)=='try-error')
        res <- NA
    else{
        # "kernel" is a similarity measure,
        # so calculate inverse for distance.
        res <- 1.0/d
        }
    
    return(res)
}

