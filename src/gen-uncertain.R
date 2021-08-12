suppressPackageStartupMessages({
    library(ggplot2); theme_set(theme_bw())
    library(sung)
})

message("gen-uncertain.R")

system("figgy=`whereis figlet`; figgylen=${#figgy}; if [ $figgylen > 10 ]; then figlet gen-uncertain.R; fi")

source("utils.R")

# Read which parameter set to use:
args <- commandArgs(trailingOnly = TRUE) # args = 4 # DEBUG

# Observed sequence:
fastafile <- "seqs/sim.fasta"

# Parameter shape for the Beta-Uniform uncertainty model:
btshp <- load_beta_shapes(fname.prm = "prm-btshp.csv", args = args)


# Load other parameters:
prm <- read.csv("prm.csv")
n_repl <- get_prm(prm, "sample.tips.n.mc")


# Deletions specifications:
# NEW - 2021-04-27, copied from sung/test-for-sung.R
prm_del <- list(alpha = 1,
                beta  = 6,
                pos   = list(c(5:15, 200:250)))

# Draw replicates from the probabilistic sequences
# defined by the Beta-Uniform uncertainty model:
seqs <- sung::draw_fasta_beta_unif(fastafile = fastafile,
                                   prm.beta  = btshp,
                                   n.repl    = n_repl,
                                   prm.del   = prm_del,
                                   alphabet.type = "nucleotide")

# Save the replicates to FASTA files:
path_seqs <- paste0("seqs/seqs-prm-", args[1], "-mc")
sung::export_fasta(seqs.replicates = seqs,
                   fasta.name = path_seqs)


# Save the total entropy
# (summed over the sequence length):
entr <- lapply(seqs, "[[", 3)
sum.entropy <- sapply(entr, sum, na.rm = T)
fname_out <- "entropy-prmset.out"
write.table(x = t(c(prmset = args[1],
                    s1 = btshp[1],
                    s2 = btshp[2],
                    mean(sum.entropy))),
            file = fname_out,
            append = TRUE, sep = ",",
            row.names = FALSE,
            col.names = FALSE, quote = FALSE)
