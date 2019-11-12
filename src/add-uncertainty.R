### 
###  Generate Probabilistic Sequences from existing sequence
###

library(seqinr)
source('utils.R')
set.seed(1234)

add_uncertainty <- function(prm, fasta.file,
                            beta.shape.p,
                            do.plot = FALSE) {
    
    # Read the simulated phylogeny:
    seqs.sim <- read.fasta(file = fasta.file, #'seqs/sim.fasta', 
                           forceDNAtolower = FALSE)
    
    n <- get_prm(prm, 'phylosim.root.seq.length')
    
    # Error probability of base call.
    s1 = beta.shape.p[1]
    s2 = beta.shape.p[2] 
    p <- rbeta(n=n, shape1 = s1, shape2 = s2)
    
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
        seq.entropy[[i]] <- apply(PS.list[[i]] , MARGIN=2, FUN=entropy)
    }
    
    if(do.plot){
        pdf('plot-add-uncertainty.pdf', width = 10)
        par(mfrow = c(2,2))
        hist(p, breaks = 30, col='grey')
        sum.entropy <- sapply(seq.entropy, sum, na.rm=T)
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
    return(PS.list)
}

prm <- read.csv('prm.csv')
prob_seqs <- add_uncertainty(prm = prm, 
                             fasta.file = 'seqs/sim.fasta', 
                             beta.shape.p = c(29, 0.1), 
                             do.plot = TRUE)

save(list = 'prob_seqs', file = 'prob_seqs.RData')
