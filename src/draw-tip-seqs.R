### 
###  Sample sequences at the tips
###  of a simulated phylogeny.
###  Sampling is made according to
###  probabilistic sequences generated
###  here from an input error probability. 
###

library(seqinr)

source('utils.R')

set.seed(1234)
prm <- read.csv('prm.csv')

# Error probability of base call.
# Fixed and the same across 
# all sequences and positions
err.prob <- get_prm(prm, 'err.proba')
p <- 1 - err.prob

# Read the simulated phylogeny:
M.list <- list()
a <- read.fasta(file = 'seqs/sim.fasta', 
                forceDNAtolower = FALSE)

# Add uncertainty by transforming the 
# (certain) sequence string into a 
# probabilistic sequence:
for(i in 1:length(a)){
    seq <- a[[i]]
    seq <- seq[seq!=' '] 
    M.list[[i]] <- build_from_seq(seq, p)
}


# Number of times we sample
# from the probabilistic sequences
# at the tips (Monte Carlo):
n.mc <- get_prm(prm, 'sample.tips.n.mc')
path.seqs <- 'seqs/seqs-mc-'
for(i in 1:n.mc){
    draw_multiple_seq(M.list,
                      filename = paste0(path.seqs,i,'.fasta')) 
    message(paste('Tip seqs sampling:',i,'/',n.mc))
}

message("Sampled tips sequences completed.")
