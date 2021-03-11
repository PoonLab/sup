# open S matrix and sample from it
# TODO: 
#   Need to access the actual sequence 
#       Search for sam->fasta conversion
#   Switch to sampling from beta posterior
#       Assumes Dirichlet prior for base probability
#       Should sometimes include Ns (current method does not)
library(ape)
library(gtools) # rdirichlet
set.seed(2112) # \m/

# Read list of files ending with .RDS
rds_names <- list.files("data/unc_covid/", pattern = "*.RDS")
# The "-1" indicates that it's a copy. Remove it.
rds_names <- rds_names[!grepl("-1", rds_names)]

# Avoid re-tracing my steps
rds_done <- sapply(list.files("data/pangolineages/"), 
    function(x) strsplit(x, "\\_")[[1]][1])
rds_todo <- sapply(strsplit(rds_names, "-"),
    function(x){
        strsplit(rev(x)[1], "\\.")[[1]][1]
    })
rds_names <- rds_names[which(!rds_todo %in% rds_done)]

# Prepare empty lists
S_list <- vector(mode = "list", length = length(rds_names))
asc_names <- c()
header_list <- S_list

# Read in uncertainty matrices and record accession names
for(i in 1:length(rds_names)){
    asc <- strsplit(rev(strsplit(rds_names[i], "-")[[1]])[1], "\\.")[[1]][1]
    asc_names <- c(asc_names, sub(".RDS", "", asc))
    
    S_list[[i]] <- readRDS(paste0("data/unc_covid/", rds_names[i]))
    names(S_list)[i] <- asc_names[i]
}
asc_names

# Check coverage of each read
#hist(apply(S_list[[1]], 1, sum), breaks = 30)

# Prepare empty list
sampled_files <- list()

# Set up timing system
t0 <- Sys.time() # Outer timer
nloops <- length(S_list)
collapsed <- estlapseds <- double(nloops)
for(i in 1:nloops){
    # Error checking
    if(is.null(dim(S_list[[i]]))) {
        print("Empty file")
        next
        }
    
    t1 <- Sys.time() # Inner timer
    
    # Normalize matrix
    S <- S_list[[i]]
    if(ncol(S) == 6) {
        S[,5] <- S[,5] + S[,6]
        S <- S[,1:5]
    }
    alph <- toupper(colnames(S))
    #print(alph)
    # Make columns add to 1
    S2sum <- as.vector(apply(S, 1, sum))
    S2mat <- matrix(rep(S2sum, ncol(S)), ncol = ncol(S), byrow = FALSE)
    S2 <- S / S2mat
    
    # Sample the sequences
    n <- 1000 # number of sampled genomes
    sampleseq_mat <- apply(S, 1, function(x) {
        if(any(is.na(x))){
            return(rep("N", n))
        } else {
            newx <- rdirichlet(1, x + rep(1/4, length(x)))
            return(sample(alph, size = n, prob = newx, replace = TRUE))
        }
    })
    
    conseq <- apply(S, 1, function(x) {
        if(any(is.na(x)) | sum(x) < 10) {
            return("N")
        } else {
            return(alph[which.max(x)])
        }
    })
    
    # Optional code for personal interest
    if(FALSE){
        seqbin <- as.DNAbin(sampleseq_mat)
        seqdist <- dist.dna(seqbin, model = "TN93")
        
        uppers <- as.numeric(seqdist)[which(upper.tri(seqdist))]
        hist(uppers, freq = FALSE)
        mean(uppers, na.rm = TRUE)
        sd(uppers, na.rm = TRUE)
        curve(dnorm(x, 
                mean(uppers, na.rm = TRUE), 
                sd(uppers, na.rm = TRUE)), 
            add = TRUE, col = 2)
        # Okay that is super weird how normal that is.
    }
    
    # Convert sample letters to single string
    conseq <- paste(conseq, collapse = "", sep = "")
    sampleseq <- c(conseq, apply(sampleseq_mat, 1, paste, collapse = ""))
    
    # Create well-formatted fasta file
    name <- paste0("> ", asc_names[i], ".", 0:length(sampleseq))
    fasta <- paste(name, sampleseq, sep = "\n", collapse = "\n")
    sampled_files[[i]] <- fasta
    writeLines(fasta, con = paste0("data/sampled_covid/", asc_names[i], "_sampled.fasta"))
    
    # Loop Timing
    elapsed <- difftime(Sys.time(), t1, units = "mins")
    collapsed[i] <- elapsed
    estlapsed <- round((nloops - i)*mean(collapsed[collapsed > 0]), 3)
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






