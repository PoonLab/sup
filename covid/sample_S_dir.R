# open S matrix and sample from it
# TODO: read available lineages from a directory.
# Will need to account for the reference sequence (possibly manually)
# Need to access the actual sequence as well
library(ape)
set.seed(2112) # \m/

rds_names <- list.files("data/unc_covid/", pattern = "*.RDS")
rds_names <- rds_names[!grepl("-1", rds_names)]
S_list <- vector(mode = "list", length = length(rds_names))
asc_names <- c()
header_list <- S_list
for(i in 1:length(rds_names)){
    asc <- substr(rds_names[i], 15, 100)
    asc_names <- c(asc_names, sub(".RDS", "", asc))
    
    S_list[[i]] <- readRDS(paste0("data/unc_covid/", rds_names[i]))
    names(S_list)[i] <- asc_names[i]
}
asc_names

hist(apply(S_list[[6]], 1, sum), breaks = 30)

sampled_files <- list()

t0 <- Sys.time()# Outer timer
nloops <- length(S_list)
collapsed <- estlapseds <- double(nloops)
for(i in 1:nloops){
    if(is.null(dim(S_list[[i]]))) next
    
    t1 <- Sys.time() # Inner timer
    
    S <- S_list[[i]]
    alph <- colnames(S)
    # Make columns add to 1
    S2sum <- as.vector(apply(S, 1, sum))
    S2mat <- matrix(rep(S2sum, ncol(S)), ncol = ncol(S), byrow = FALSE)
    S2 <- S / S2mat
    
    # Sample the sequences
    n <- 1000 # number of sampled genomes
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
    #refseq <- readLines("data/MN908947-3.fasta")
    # Add reference sequence as the first row (easy to compare to others)
    #refseq <- paste(c(refseq[1], paste0(refseq[-1], collapse = "")), collapse = "\n")
    
    # Create well-formatted fasta file
    con <- "data/SRR13020989_small.sam"
    name <- paste0("> ", asc_names[i], ".", 1:length(sampleseq))
    fasta <- paste(name, sampleseq, sep = "\n", collapse = "\n")
    sampled_files[[i]] <- fasta
    #writeLines(fasta, con = paste0("data/sampled_covid/", asc_names[i], "_sampled.fasta"))
    
    # Loop Timing
    elapsed <- difftime(Sys.time(), t1, units = "mins")
    collapsed[i] <- elapsed
    estlapsed <- round((nloops - i)*mean(collapsed), 3)
    estlapseds[i] <- estlapsed
    if(FALSE){ # Optional
        par(mfrow = c(1,2))
        plot(collapsed[!is.na(collapsed)])
        plot(estlapseds[!is.na(estlapsed)]/60)
    }
    print(paste0("loop ", i, " of ", nloops, ", ", round(elapsed, 3),
        " mins | Approx. ", estlapsed, " mins (", round(estlapsed/60, 2),
        " hours) remaining"))
}

# saveRDS(sampleseq, file = "data/SRR13020989_sampled.rds")

# # Commands for shell (can't be run from RStudio because RStudio uses the wrong path)
# conda activate pangolin
# pangolin data/SRR13020989_sampled.fasta --outfile data/SRR13020989_pangolineages.csv
# # for 1000 sequences, took 30 seconds

for(i in 1:length(asc_names)){
    cat(paste0(
        "pangolin data/sampled_covid/", asc_names[i], "_sampled.fasta",
        " --outfile data/pangolineages/", asc_names[i], "_pangolineages.csv"
    ))
    cat("\n")
}
cat('spd-say "Hey devan your code is done"')




