### 
###  Generate Probabilistic Sequences from existing sequence
###

suppressPackageStartupMessages({
library(seqinr);
library(tidyr);
library(dplyr);
library(ggplot2) ; theme_set(theme_bw())
library(stringr);
})


source('utils.R')

args <- commandArgs(trailingOnly = TRUE)  # args=4

# ---- Functions ----

add_uncertainty <- function(prm, fasta.file,
                            beta.shape.p,
                            prmset,
                            do.plot = FALSE,
                            fname.out = 'entropy-prmset.out') {
    
    # Read the simulated phylogeny:
    seqs.sim <- read.fasta(file = fasta.file, #'seqs/sim.fasta', 
                           forceDNAtolower = FALSE)
    
    # Retrieve the sequences names:
    seqs.names <- sapply(seqs.sim, attr, which='Annot') %>%
        sapply(stringr::str_extract, pattern='\\w+_\\w+')
    
    n <- get_prm(prm, 'phylosim.root.seq.length')
    
    # Error probability of base call.
    s1 <- beta.shape.p[1]
    s2 <- beta.shape.p[2] 
    p  <- rbeta(n=n, shape1 = s1, shape2 = s2)
    
    # Add uncertainty by transforming the 
    # (certain) sequence string into a 
    # probabilistic sequence. 
    # Also calculate entropy.
    PS.list     <- list()
    seq.entropy <- list()
    for(i in 1:length(seqs.sim)){  #i=1
        seq <- seqs.sim[[i]]
        seq <- seq[seq!=' '] 
        PS.list[[i]]     <- build_from_seq(seq, p)
        seq.entropy[[i]] <- apply(PS.list[[i]] , 
                                  MARGIN=2, 
                                  FUN=entropy)
    }
    
    # Save the mean entropy 
    # (summed over the sequence length):
    sum.entropy <- sapply(seq.entropy, sum, na.rm=T)
    write.table(x = t(c(prmset,s1,s2,mean(sum.entropy))), 
                file = fname.out,
                append = TRUE, sep = ',', 
                row.names = FALSE, 
                col.names = FALSE, quote = FALSE)
    
    if(do.plot){
        plotname <- paste0('plot-add-uncertainty-',args[1],'.pdf')
        pdf(plotname, 
            width = 10)
        par(mfrow = c(2,2))
        hist(p, breaks = 30, col='grey',
             main = paste('Base call probability\nbeta shape =',
                          paste(beta.shape.p,collapse = ' ; ')))
        boxplot(sum.entropy, 
                las=1,
                main = paste('Total entropy.\nSequence length =',n))
        hist(unlist(seq.entropy), 
             col = 'lightgray',
             breaks = 40, 
             main = 'Distribution of position entropy',
             xlab='entropy')
        jj = sample(1:length(seq.entropy), size = 1)
        plot(seq.entropy[[jj]], 
             las=1,
             main = paste('Entropy by position\nfor seq #',jj),
             xlab = 'position', ylab = 'entropy')
        dev.off()
        
    }
    return(list(PS = PS.list, 
                seqs.names = seqs.names))
}

plot_prmset_distrib <- function(fname) { #fname='prm-btshp.csv'
    
    b <- read.csv(fname, header = F)
    
    dfl <- list()
    for(i in 1:nrow(b)){
        x <- seq(0,1,length=1e3)
        y <- dbeta(x, shape1 = b[i,1], shape2 = b[i,2])
        dfl[[i]] <- data.frame(x=x, y=y, 
                               prmset=paste('prmset',i,
                                            ': s1 =',b[i,1],
                                            '; s2 =',b[i,2]))
    }
    df <- do.call('rbind',dfl) %>%
        mutate(prmset = factor(prmset))
    
    g <- df %>%
        ggplot(aes(x=x, y=y))+
        geom_line(aes(colour=prmset), size=2, alpha=0.8)+
        scale_y_log10(limits = c(1e-12, 1e3))+
        coord_cartesian(xlim = c(0.5,1))+
        ggtitle('Beta Density of Base Call Probability')+
        xlab('Base call probability') + ylab('density')
    plot(g)
 }


# ---- RUN ----

prm <- read.csv('prm.csv')

if(length(args) == 0){
    btshp <- c(get_prm(prm, 'beta.shape.p1'), 
               get_prm(prm, 'beta.shape.p2'))
}

if(length(args) == 1){
    fname = 'prm-btshp.csv'
    idx <- as.numeric(args[1])
    btshp.file <- read.csv(fname, header = F)
    btshp <- as.numeric(btshp.file[idx,])
    pdf('plot-proba-basecall-beta.pdf')
    plot_prmset_distrib(fname)
    dev.off()
}

print(paste('prm set =', args[1],'-->',
            paste('beta shape =',btshp, 
            collapse = ' ; ')))

prob_seqs <- add_uncertainty(prm = prm, 
                             fasta.file = 'seqs/sim.fasta', 
                             beta.shape.p = btshp,
                             prmset = args[1],
                             do.plot = TRUE)

fname <- paste0('prob_seqs_',args[1],'.RData')
save(list = 'prob_seqs', file = fname)
