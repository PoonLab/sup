library(here)
library(gtools) # rdirichlet
source("aux_funk.R")

sams <- list.files("samDONE")
sams <- unique(sapply(strsplit(sams, split = "-|_|\\."), function(x) x[2]))

args <- commandArgs(TRUE)
N <- 10
if ("-N" %in% args) {
    N <- args[which(args == "-N") + 1]
    N <- as.numeric(N)
}
dirich <- "-d" %in% args

metadl <- read.csv("sequences_descr_mt_downloaded.csv")

rds_names <- list.files("samDONE")
rds_names <- rds_names[!grepl("insertions", rds_names)]
asc_names <- parse_accession(rds_names) # aux_funk.R

if (!dir.exists("sampled_trees")) dir.create("sampled_trees")
for (i in seq_len(N)) {
    sample_file <- paste0("sampled_trees/sampled_tree_", i, ".fasta")
    if (file.exists(sample_file)) {
        file.remove(sample_file)
    }
}




# Set up timing system --------------------------
t0 <- Sys.time() # Outer timer
nloops <- length(rds_names)
collapsed <- estlapseds <- double(nloops)

for (i in seq_len(nloops)) { 
    t1 <- Sys.time()

	S <- readRDS(here("RTT", "samDONE", rds_names[i]))
	S <- fix_unc(S) # aux_funk.R
    if ("character" %in% class(S)) {
        print(paste("S", asc_names[i], sep = " - "))
        next
    }

    alph <- toupper(colnames(S))

    # Sample the sequences
    sampleseq_mat <- apply(S, 1, function(x) {
        if (any(is.na(x)) | (sum(x) < 10)) {
            return(rep("N", N))
        } else {
            if (dirich) {
                newx <- rdirichlet(1,
                    as.numeric(x + c(rep(1 / 4, 4),
                    rep(0, length(x) - 4))))
                return(sample(alph, size = N, prob = newx, replace = TRUE))
            } else {
                return(sample(alph, size = N, prob = x, replace = TRUE))
            }

        }
    })

    sampleseq <- apply(sampleseq_mat, 1, paste, collapse = "")
    
    asc_descr <- metadl$def[which(metadl$SRAAccession == asc_names[i])]
    first_colon <- gregexpr(text = asc_descr, pattern = ":")[[1]][1]
    asc_descr <- substr(asc_descr, first_colon + 1, nchar(asc_descr))
    asc_descr <- gsub(x = asc_descr, pattern = ":", replacement = "_")
    asc_descr <- gsub(x = asc_descr, pattern = ",", replacement = "_")
    asc_descr <- gsub(x = asc_descr, pattern = " ", replacement = "_")
    

    # Append to each file
    for (ii in seq_len(N)) {
        sample_file <- paste0("sampled_trees/sampled_tree_", ii, ".fasta")
        if (!file.exists(sample_file)) file.create(sample_file)
        fasta_string <- paste(asc_descr, sampleseq[ii], "\n", sep = "\n", collapse = "\n")
        cat(fasta_string, append = TRUE, file = sample_file)
    }

    # Loop Timing
    elapsed <- difftime(Sys.time(), t1, units = "mins")
    collapsed[i] <- elapsed
    estlapsed <- round((nloops - i) * mean(collapsed[collapsed > 0]), 3)
    estlapseds[i] <- estlapsed
    print(paste0("loop ", i, " of ", nloops, ", ", round(elapsed, 3),
        " mins | Approx. ", estlapsed, " mins (", round(estlapsed / 60, 2),
        " hours) remaining"))
}


