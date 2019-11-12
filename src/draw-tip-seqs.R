### 
###  Sample sequences at the tips
###  of a simulated phylogeny.
###  Sampling is made according to
###  probabilistic sequences generated
###  here from an input error probability. 
###

library(seqinr)

source('utils.R')
source('add-uncertainty.R')
set.seed(1234)

# Load the probabilistic sequences object: `prob_seqs`
load('prob_seqs.RData')


# Number of times we sample
# from the probabilistic sequences
# at the tips (Monte Carlo).
# Each MC draw saves a FASTA file 
# for the associated sequence.
n.mc      <- get_prm(prm, 'sample.tips.n.mc')
path.seqs <- 'seqs/seqs-mc-'
for(i in 1:n.mc){
    draw_multiple_seq(prob_seqs,
                      filename = paste0(path.seqs,i,'.fasta')) 
    message(paste('Tip seqs sampling:',i,'/',n.mc))
}

message("Sampled tips sequences completed.")
