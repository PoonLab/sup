###
###  Fit shape parameters of Beta distribution 
###  for $p$, which is the probability of 
###  observing a given nucleotide. 
###  The Fit is made on the entropy of Zanini's data. 
###


library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)

source('utils.R')

set.seed(1234)

# ---- Function ----

#' simulate entropy for all position of a sequence
#' assuming the base call probability has a 
#' Beta distribution.
simul_entropy <- function(beta_shapes, sq) {
    n <- length(sq)
    p <- rbeta(n = n, 
               shape1 = beta_shapes[1], 
               shape2=beta_shapes[2])
    
    df <- data.frame(base = sq,
                     p     = p, 
                     stringsAsFactors = FALSE)
    ps <- pmaxseq_to_probseq(df, 
                             test.mode = TRUE)
    
    sim.e <- apply(ps, 1, entropy)
    return(sim.e)
}

error_fct <- function(x, zan.entropy, sq, 
                      doplot = FALSE, 
                      plot.title='') {
    
    brks  <- seq(0,2,by = 0.1)
    
    sim.e <- simul_entropy(beta_shapes = exp(x),  # <-- For numerical stability
                           sq = sq)
    h     <- hist(sim.e, plot = FALSE, breaks = brks)
    h.zan <- hist(zan.entropy, plot = FALSE, breaks = brks)
    
    ds <- h$density
    dz <- h.zan$density
    
    if(doplot){
        plot(x=h.zan$mids, y=dz, 
             typ='b', las=1, log='y',
             xlab = 'entropy',ylab='density',
             ylim=c(1e-5, max(ds,dz)),
             main = plot.title)
        lines(x=h$mids, y=ds, 
              col='blue', typ='o', pch=16)
    }
    
    ds <- ds[ds>0]
    dz <- dz[dz>0]
    ns <- length(ds)    
    nz <- length(dz)    
    if(ns<nz) dz <- dz[1:ns]
    if(nz<ns) ds <- ds[1:nz]
    
    return( sqrt(sum((log(ds/dz))^2)) / min(nz,ns) )
}

#' Extract the entropy from Zanini's data, 
#' for given a patient and time point.
extract_pat_tp_e <- function(pat, tpt, zanini.entropy) {
    tmp <- zanini.entropy %>%
        filter(patient==pat & tp == tpt)
    return(tmp$entropy)
}

fit_beta_entropy <- function(pat, tp, sq, zan.entropy) {
    
    thefit <- optim(par = c(3,-3), # <--- Exp() 
                    fn = error_fct, 
                    zan.entropy = zan.entropy, 
                    sq = sq, 
                    method = 'Nelder-Mead',
                    control = list(trace  = 0, 
                                   maxit  = 100,
                                   reltol = 1e-4))
    return(exp(thefit$par))
}

fit_one <- function(pat, tpt, zanini.entropy) {
    # Extract entropy associated 
    # to patient & timepoint: 
    zan.entropy <- extract_pat_tp_e(pat, tpt, zanini.entropy)
    
    # Randomly generate a sequence 
    # of the same length:
    nucleo <- c('A','C','G','T')
    sq <- nucleo[sample(1:4, 
                        size = length(zan.entropy), 
                        replace = TRUE)]
    
    # Fit shape parameters of Beta distribution
    # for the observation probability of nucleotides:
    pf <- fit_beta_entropy(pat, tp, sq, zan.entropy)
    
    # Return result:
    return(list(pat = pat,
                tpt = tpt, 
                betashape = pf))
}

#' Fit shape parameters of the base call probability 
#' which is Beta distributed for several patients and timepoints.
fit_all <- function(pat.vec, tpt.vec, zanini.entropy) {
    res <- list() ; k=1
    for(pat in pat.vec){
        for(tpt in tpt.vec){
            print(paste('Beta entropy fit: pat =',pat,'; t =',tpt,
                        ';',k,'/',length(pat.vec)*length(tpt.vec)))
            res[[k]] <- fit_one(pat, tpt, zanini.entropy)
            k=k+1
        }
    }
    return(res)
}

# ---- RUN ----

# Load `zanini.entropy`:
load('../data/zanini/zanini-entropy.RData')

pat.vec <- c(1,5, 9)  # c(1:11)
tpt.vec <- c(1, 0)    # c(1:5, 0)

betafit <- fit_all(pat.vec, tpt.vec, zanini.entropy)

x <- lapply(betafit, '[[', 'betashape')
beta_1 <- sapply(x,'[[',1)
beta_2 <- sapply(x,'[[',2)
print('\nFitted Beta_1:')
print(summary(beta_1))
print('\nFitted Beta_2:')
print(summary(beta_2))

save.image(file = 'fit-entropy-beta.RData')

# - - - Plots 

pdf('plot-fit-entropy-beta-zanini.pdf', 
    width=20, height = 20)

par(mfrow=c(1,2))
hist(beta_1, col='grey')
abline(v=mean(beta_1), lty=2)
hist(beta_2, col='grey')
abline(v=mean(beta_2), lty=2)

par(mfrow=c(length(pat.vec),length(tpt.vec)))

for(i in 1:length(betafit)){ #i=1
    pat <- betafit[[i]]$pat
    tpt <- betafit[[i]]$tpt
    zan.entropy <- extract_pat_tp_e(pat = pat, 
                                    tpt = tpt, 
                                    zanini.entropy)
    nucleo <- c('A','C','G','T')
    sq <- nucleo[sample(1:4, 
                        size = length(zan.entropy), 
                        replace = TRUE)]
    
    error_fct(x = log(betafit[[i]]$betashape), 
              zan.entropy,
              sq, 
              doplot = TRUE ,
              plot.title = paste('patient:',pat, 'Timepoint:',tpt))
}
dev.off()
