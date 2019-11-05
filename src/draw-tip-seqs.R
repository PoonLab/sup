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

# Read the simulated phylogeny:
M.list <- list()
a <- read.fasta(file = 'seqs/sim.fasta', 
                forceDNAtolower = FALSE)

# Error probability of base call.
# Fixed and the same across 
# all sequences and positions
err.prob <- get_prm(prm, 'err.proba')
s1 = 0.1
s2 = s1/err.prob-s1
# s1 = s2 / (1-err.prob)
n = get_prm(prm, 'phylosim.root.seq.length')
err.prob <- rbeta(n=n, shape1 = s1, shape2 = s2)
hist(err.prob)
p <- 1 - err.prob


# Add uncertainty by transforming the 
# (certain) sequence string into a 
# probabilistic sequence:
seq.entropy <- list()
for(i in 1:length(a)){  #i=1
    seq <- a[[i]]
    seq <- seq[seq!=' '] 
    M.list[[i]]<- build_from_seq(seq, p)
    seq.entropy[[i]] <- apply(M.list[[i]] , MARGIN=2, FUN=entropy)
}
sum.entropy <- sapply(seq.entropy, sum, na.rm=T)
sum.entropy
plot(seq.entropy[[1]])

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
