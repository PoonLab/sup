# open S matrix and sample from it
library(ape)
library(gtools) # rdirichlet
library(here)
#set.seed(2112) # \m/
source(here("src", "aux_funk.R"))
library(dplyr)

args <- commandArgs(TRUE)
# --overwrite: Resample ALL
# -N number of resamples
# -d dirichlet prior
N <- 10
if ("-N" %in% args) {
    N <- args[which(args == "-N") + 1]
    N <- as.numeric(N)
}

append <- FALSE
if (!"--overwrite" %in% args) {
    append <- TRUE
}

dirich <- "-d" %in% args

print(args)

# Read list of files ending with .RDS
rds_names <- list.files(here("data", "output"), pattern = "*csv")
rds_names <- rds_names[!grepl(rds_names, pattern = "_insertions")]

# Record accession names
asc_names <- unique(parse_accession(rds_names)) # aux_funk.R
print(asc_names)

# Prepare empty lists
S_list <- vector(mode = "list", length = length(rds_names))
header_list <- S_list

# Read in uncertainty matrices
for (i in seq_along(rds_names)) {
    S_list[[i]] <- read.csv(here("data", "output", rds_names[i]))
    names(S_list)[i] <- asc_names[i]
}


# Set up timing system --------------------------
t0 <- Sys.time() # Outer timer
nloops <- length(S_list)
collapsed <- estlapseds <- double(nloops)
for (i in 1:nloops) {
    t1 <- Sys.time() # Inner timer
    print(asc_names[i])

    # Normalize matrix
    S <- S_list[[i]]
    S <- fix_unc(S) # aux_funk.R
    if ("character" %in% class(S)) {
        print(paste("S-", asc_names[i], " had characters in the matrix", sep = ""))
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
                    as.numeric(x))
                return(sample(alph, size = N, prob = newx, replace = TRUE))
            } else {
                return(sample(alph, size = N, prob = x, replace = TRUE))
            }

        }
    })

    sampleseq <- lapply(seq_len(nrow(sampleseq_mat)), function(x) sampleseq_mat[x, ])

    conseq <- apply(S, 1, function(x) {
        if (any(is.na(x)) | sum(x) < 10) {
            return("N")
        } else {
            return(alph[which.max(x)])
        }
    })

    filename <- paste0("data/", "sampled_covid/", asc_names[i],
            "_sampled", ifelse(dirich, "_d", ""), ".fasta")

    # Convert sample letters to single string
    if (!append) {
        conseq <- paste(conseq, collapse = "", sep = "")
        sampleseq <- sapply(sampleseq, paste, collapse = "")
        sampleseq <- c(conseq, sampleseq)

        # Create well-formatted fasta file
        name <- paste0("> ", asc_names[i], ".", 0:(length(sampleseq) - 1))
        fasta <- paste(name, sampleseq, sep = "\n", collapse = "\n")
        file_line_count <- length(sampleseq)
    } else {
        # When appending, conseq is not needed.
        file_line_count <- system(paste0("wc -l ", filename), intern = TRUE)
        file_line_count <- as.numeric(strsplit(file_line_count, " ")[[1]][1])
        sampleseq <- sapply(sampleseq, paste, collapse = "")
        name <- paste0("> ", asc_names[i], ".", 1:length(sampleseq) + file_line_count)
        fasta <- paste(name, sampleseq, sep = "\n", collapse = "\n")
    }
    cat(fasta, append = append,
        file = filename)

    # Loop Timing
    elapsed <- difftime(Sys.time(), t1, units = "mins")
    collapsed[i] <- elapsed
    estlapsed <- round((nloops - i) * mean(collapsed[collapsed > 0]), 3)
    estlapseds[i] <- estlapsed
    print(paste0("loop ", i, " of ", nloops, ", ",
        file_line_count, " lines per file, ", 
        round(elapsed, 3), " mins | Approx. ", estlapsed, 
        " mins (", round(estlapsed / 60, 2), " hours) remaining"))
}

#print(asc_names)
writeLines(asc_names, con = here("data", "sampled_covid", "ascnames.txt"))
