# open S matrix and sample from it
# TODO: read available lineages from a directory.
    # Will need to account for the reference sequence (possibly manually)
    # Need to access the actual sequence as well
library(ape)

asc <- "ERR4363387"
S <- readRDS(paste0("data/unc_covid/open_seSAMe-S-", asc, ".RDS"))
alph <- names(S)

hist(apply(S, 1, sum))

# Make columns add to 1
S2sum <- as.vector(apply(S, 1, sum))
S2mat <- matrix(rep(S2sum, ncol(S)), ncol = ncol(S), byrow = FALSE)
S2 <- S / S2mat

# Sample the sequences
n <- 5000
set.seed(1000)
sampleseq_mat <- apply(S2, 1, function(x) {
    if(any(is.na(x))){
        return(rep("N", n))
    } else {
        return(sample(alph, size = n, prob = x, replace = TRUE))
    }
})

if(FALSE){
    seqbin <- as.DNAbin(sampleseq_mat)
    seqdist <- dist.dna(seqbin, model = "TN93")
    
    uppers <- as.numeric(seqdist)[which(upper.tri(seqdist))]
    hist(uppers, freq = FALSE)
    mean(uppers, na.rm = TRUE)
    sd(uppers, na.rm = TRUE)
    curve(dnorm(x, mean(uppers, na.rm = TRUE), sd(uppers, na.rm = TRUE)), add = TRUE, col = 2)
    # Okay that is super weird how normal that is.
}

# Convert sample letters to single string
sampleseq <- apply(sampleseq_mat, 1, paste, collapse = "")

# Look up the reference sequence
refseq <- readLines("data/MN908947-3.fasta")
# Add reference sequence as the first row (easy to compare to others)
refseq <- paste(c(refseq[1], paste0(refseq[-1], collapse = "")), collapse = "\n")

# Create well-formatted fasta file
con <- "data/SRR13020989_small.sam"
name <- paste0("> ", strsplit(readLines(con, n=4)[2], split = "\t")[[1]][2], 1:length(sampleseq))
fasta <- paste(refseq, name, sampleseq, sep = "\n", collapse = "\n")
writeLines(fasta, con = "data/SRR13020989_sampled.fasta")

# saveRDS(sampleseq, file = "data/SRR13020989_sampled.rds")

# # Commands for shell (can't be run from RStudio because RStudio uses the wrong path)
# conda activate pangolin
# pangolin data/SRR13020989_sampled.fasta --outfile data/SRR13020989_pangolineages.csv
# # for 1000 sequences, took 30 seconds






