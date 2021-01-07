# open S matrix and sample from it
library(ape)

asc <- "ERR4085809"

S <- readRDS(paste0("data/S/open_seSAMe-S-", asc, ".RDS"))

alph <- names(S)

hist(apply(S, 1, sum))

S2sum <- as.vector(apply(S, 1, sum))
S2mat <- matrix(rep(S2sum, ncol(S)), ncol = ncol(S), byrow = FALSE)

S2 <- S / S2mat

n <- 5000
set.seed(1000)
sampleseq_mat <- apply(S2, 1, function(x) {
    if(any(is.na(x))){
        return(rep("N", n))
    } else {
        return(sample(alph, size = n, prob = x, replace = TRUE))
    }
})

'
seqbin <- as.DNAbin(sampleseq_mat)
seqdist <- dist.dna(seqbin, model = "TN93")

uppers <- as.numeric(seqdist)[which(upper.tri(seqdist))]
hist(uppers, freq = FALSE)
mean(uppers, na.rm = TRUE)
sd(uppers, na.rm = TRUE)
curve(dnorm(x, mean(uppers, na.rm = TRUE), sd(uppers, na.rm = TRUE)), add = TRUE, col = 2)
# Okay that is super weird how normal that is.
'

sampleseq <- apply(sampleseq_mat, 1, paste, collapse = "")

# TODO: Automatic lookup

refseq <- readLines("data/MN908947-3.fasta")
refseq <- paste(c(refseq[1], paste0(refseq[-1], collapse = "")), collapse = "\n")

con <- "data/SRR13020989_small.sam"
name <- paste0("> ", strsplit(readLines(con, n=4)[2], split = "\t")[[1]][2], 1:length(sampleseq))
fasta <- paste(refseq, name, sampleseq, sep = "\n", collapse = "\n")
writeLines(fasta, con = "data/SRR13020989_sampled.fasta")

# saveRDS(sampleseq, file = "data/SRR13020989_sampled.rds")

# conda activate pangolin
# pangolin data/SRR13020989_sampled.fasta --outfile data/SRR13020989_pangolineages.csv

# for 1000, took 30 seconds






