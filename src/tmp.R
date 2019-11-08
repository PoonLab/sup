library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)

source('utils.R')



simul_entropy <- function(beta_shapes, sq) {
    n <- length(sq)
    p <- rbeta(n = n, 
               shape1 = beta_shapes[1], 
               shape2=beta_shapes[2])
    
    df <- data.frame(base = sq,
                     p     = p, 
                     stringsAsFactors = FALSE)
    ps <- pmaxseq_to_probseq(df)
    
    sim.e <- apply(ps, 1, entropy)
    return(sim.e)
}


    
error_fct <- function(x, zan.entropy, sq, doplot = FALSE) {
    
    brks  <- seq(0,2,by = 0.1)
    
    sim.e <- simul_entropy(beta_shapes = abs(x),  # <-- FIXME abs()
                           sq = sq)
    h     <- hist(sim.e, plot = FALSE, breaks = brks)
    h.zan <- hist(zan.entropy, plot = FALSE, breaks = brks)
    
    ds <- h$density
    dz <- h.zan$density
    
    if(doplot){
        plot(x=h.zan$mids, y=dz, 
             typ='b', las=1, log='y',
             xlab = 'entropy',ylab='density',
             ylim=c(1e-5, max(ds,dz)))
        lines(x=h$mids, y=ds)
    }
    
    ds <- ds[ds>0]
    dz <- dz[dz>0]
    ns <- length(ds)    
    nz <- length(dz)    
    if(ns<nz) dz <- dz[1:ns]
    if(nz<ns) ds <- ds[1:nz]
    
    return( sqrt(sum((log(ds/dz))^2)) / min(nz,ns) )
}


# ---- RUN ----

load('../data/zanini/zanini-entropy.RData')

pat = 9
tpt = 2
zan <- zanini.entropy %>%
    filter(patient==pat & tp == tpt)


nucleo <- c('A','C','G','T')
sq <- nucleo[sample(1:4, size = nrow(zan), replace = TRUE)]


error_fct(x = c(20, .05), 
          zan.entropy = zan$entropy, sq, doplot = T)

thefit <- optim(par = c(20,0.05), 
                fn = error_fct, 
                zan.entropy = zan$entropy, 
                sq = sq, 
                control = list(trace  = 2, 
                               maxit  = 20,
                               reltol = 1e-4))

pfit <- thefit$par

error_fct(x=pfit, zan.entropy = zan$entropy, 
          sq, doplot=TRUE)

