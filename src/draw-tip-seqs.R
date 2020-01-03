### 
###  Sample sequences at the tips
###  of a simulated phylogeny.
###  Sampling is made according to
###  probabilistic sequences generated
###  here from an input error probability. 
###

suppressPackageStartupMessages(library(seqinr))
source('utils.R')

set.seed(1234)
args <- commandArgs(trailingOnly = TRUE)   # args=4

# Load the probabilistic sequences object: `prob_seqs`
load(paste0('prob_seqs_',args[1],'.RData'))


# Number of times we sample
# from the probabilistic sequences
# at the tips (Monte Carlo).
# Each MC draw saves a FASTA file 
# for the associated sequence.
prm       <- read.csv('prm.csv')
n.mc      <- get_prm(prm, 'sample.tips.n.mc')
path.seqs <- paste0('seqs/seqs-prm-',args[1],'-mc-')

for(i in 1:n.mc){
    draw_multiple_seq(prob_seqs,
                      filename = paste0(path.seqs,i,'.fasta')) 
    if(i%%10==0) message(paste('prm set #',args[1], 
                               '--> Tip seqs sampling:',i,'/',n.mc))
}

message(paste("Sampled tips sequences completed. (prm set =",args[1],
              ')'))
