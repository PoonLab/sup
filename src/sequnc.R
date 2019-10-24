
library(seqinr)

source('utils.R')


# ---- RUN ----

set.seed(1234)

# Number of times we are going to 
# reconstruct the phylogenic tree (Monte Carlo):
n.trees <- 25

# Fixed probability 
p <- 0.9


# From a simulated phylogeny:
M.list.sim <- list()
a <- read.fasta(file = 'seqs/sim.fasta', 
                forceDNAtolower = FALSE)
attr(a[[1]],'Annot')
for(i in 1:length(a)){
    seq <- a[[i]]
    seq <- seq[seq!=' '] 
    M.list.sim[[i]] <- build_from_seq(seq = seq, p)
}

# Choose here if you want randomly built
# or from simulations:
M.list = M.list.sim    # M.list.sim  M.list.rand
for(i in 1:n.trees){
    draw_multiple_seq(M.list,
                      filename = paste0('seqs/seqs-',i,'.fasta')) 
}


