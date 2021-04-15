library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)
library(ape)

set.seed(123456)

# ---- Functions ----

#' Swap 2 tips of a given phylogeny.
#' @param phy Phylo object
#' @param i Integer. Tip number.
#' @param j Integer. Tip number.
swap_tips <- function(phy, i, j) {
    old.lab <- phy$tip.label
    new.lab <- replace(old.lab, c(i,j), old.lab[c(j,i)])
    phy.new <- phy
    phy.new$tip.label <- new.lab
    return(phy.new)
}

#' Swap all possible tip pairs.
swap_all_tip_pairs <- function(phy) {
    n <- length(phy$tip.label)
    ivec <- numeric(n*(n+1)/2)
    jvec <- numeric(n*(n+1)/2)
    d.rf <- numeric(n*(n+1)/2)
    k <- 1
    for(i in 1:n){
        for(j in 1:i){
            tmp <- swap_tips(phy, i, j)
            ivec[k] <- i
            jvec[k] <- j
            d.rf[k] <- dist.topo(phy, tmp)
            k <- k+1
        }
    }
    df <- data.frame(ivec, jvec, d.rf) %>%
        mutate(dij = (ivec-jvec)) 
    
    g1 <- df %>%
        ggplot(aes(x=ivec, y=jvec))+
        geom_tile(aes(fill=d.rf)) + 
        scale_fill_gradient(low='lightgreen', high='red')
    
    g2 <- df %>%
        ggplot(aes(x=dij, y=d.rf))+
        geom_point(alpha=0.3, size=3) 
    
    grid.arrange(g1, g2)
}

#' Swap p tip pairs randomly chosen.
#' @param phy Phylo object
#' @param p Integer. Number of tip pairs to swap.
random_p_swaps <- function(phy, p) { # p = 3
    n <- length(phy$tip.label)
    if(2*p > n) stop("ERROR: too many swaps asked!")
    idx <- sample(1:n, 2*p)
    ivec <- idx[1:p]
    jvec <- idx[(p+1):(2*p)]
    
    tmp <- phy
    for(k in 1:p){
        tmp <- swap_tips(tmp, ivec[k], jvec[k])
    }
    return(tmp)
}

#' Calculate the RF distances b/w a benchmark phylo
#' and phylogenies that had some tip pairs swapped.
#' 
#' @param phy Phylo object
#' @param iter.per.swap Integer. Number of times the tip pairs are swapped
#' @param swapmax Integer. The number of swaps will range from 1 to `swapmax`
multi_swap <- function(phy, 
                       iter.per.swap, 
                       swapmax ) {
    m    <- iter.per.swap * swapmax
    d.rf <- numeric(m)
    sw   <- numeric(m)
    iter <- numeric(m)
    k = 1
    for(p in 1:swapmax){
        print(paste(p,'/',swapmax))
        for(i in 1:iter.per.swap){
            tmp     <- random_p_swaps(phy, p)
            d.rf[k] <- dist.topo(phy, tmp)
            sw[k]   <- p
            iter[k] <- i
            k = k+1  
        }
    }
    df <- data.frame(sw, iter, d.rf)
    return(df)
}


# ---- RUN ----

n   <- 75   # number of tips
phy <- unroot(rtree(n))

df <- multi_swap(phy, 
                 iter.per.swap = 30, 
                 swapmax = round(n/2)-1) 

g <- df %>%
    ggplot()+
    geom_boxplot(aes(x=factor(sw), y=d.rf),fill='lightgrey')+
    geom_hline(yintercept=2*(n-1)-3, linetype='dashed')+
    xlab('number of pair swaps')+
    ylab('RF distance')+
    ggtitle(paste('Tree with',n,'tips'))
plot(g)


