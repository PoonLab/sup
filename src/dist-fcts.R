###
###  DEFINE THE DISTANCE FUNCTIONS
### 

suppressPackageStartupMessages({
    library(dplyr);
    library(ggplot2);
    library(ape);
    library(phytools);
    library(phylobase)})
source('treekernel.R')

# ---- Tree distances ----

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

#' Kuhner-Felsenstein distance. 
#' @param tree1 `ape::phylo` object. First tree. 
#' @param tree2 `ape::phylo` object. Second tree. 
dist.KF <- function(tree1, tree2) {
    # Distance b/w tree1 and tree2:
    res <- ape::dist.topo(tree1, tree2, method='score')
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


# ---- TN93 distance ----


.extract_mc <- function(s) {
    return(stringr::str_extract(s, "mc-\\d+") %>%
        stringr::str_extract("\\d+") %>%
        as.numeric())
}

#' Retrieve the TN93 distances (b/w sequences of one single tree)
#'
#' @param prmset Integer. Parameter set used.
.dist.tn93.prmset <- function(prmset) {
    # prmset = 1 ; mc = 7
    
    # WARNING: file name template hard-coded!
    fastafiles <- system(paste0('ls seqs/seqs-prm-', 
                                prmset,
                                '*.fasta.out'),
                         intern = TRUE)    
    n <- length(fastafiles)
    if(n==0){
        stop(paste("CANNOT FIND seqs/*.fasta.out FOR PRM SET #",
                   prmset))
    }
    tmp <- list()
    for(i in 1:n){  # i=2
        mc <- .extract_mc(fastafiles[i])
        tmp[[i]] <- read.csv(fastafiles[i])
        tmp[[i]]$mc <- mc
        tmp[[i]]$prmset <- prmset
    }
    df <- do.call('rbind', tmp)
    return(df)
}

#' Retrieve all TN93 distances calculated
#' by the script `dist-tn93.sh` and store
#' them in a dataframe
dist.tn93 <- function() {
    n <- system('wc -l < prm-btshp.csv', intern = TRUE) %>% 
        as.numeric()
    tmp <- list()
    for(i in 1:n){
        tmp[[i]] <- .dist.tn93.prmset(i)
    }
    df <- do.call('rbind', tmp)
    return(df)
}



.num_of_clusters <- function(df, p, m, dist.thresh.mean) {
    # m=2; p=3
    dfi <- df %>%
        filter(mc==m, prmset==p) %>%
        filter(Distance < dist.thresh.mean * mean(Distance)) %>%
        select(ID1,ID2)
    g <- igraph::graph.data.frame(dfi, directed = FALSE)
    n <- igraph::no.clusters(g)
    return(n)
}

#' Calculate the number of clusters for each 
#' reconstructed tree across all parameter sets 
#' and MC iterations.
#' @param df Data frame as output by the function `dist.tn93()`
#' @param dist.thresh.mean Numeric. Threshold TN93 distance for defining a cluster, expressed as a proportion of the mean TN93 distances of one reconstructed tree (i.e., 0.1 will set the threshold at 0.1*mean).
clstr_num <- function(df, 
                      dist.thresh.mean) {
    
    mc         <- numeric()
    prmset     <- numeric()
    n.clusters <- numeric()
    k = 1
    for(p in unique(df$prmset)){
        print(paste('Calculating number of clusters', 
                    'for prmset #',
                    p))
        for(m in unique(df$mc)){
            mc[k] <- m
            prmset[k] <- p
            n.clusters[k] <- .num_of_clusters(df, 
                                              p, m, 
                                              dist.thresh.mean)
            k <- k+1   
        }
    }
    dfclst <- data.frame(prmset, mc, n.clusters)
    return(dfclst)
}


