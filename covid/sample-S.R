# open S matrix and sample from it
# TODO:
#   Need to access the actual sequence
#       Search for sam->fasta conversion
#   Switch to sampling from beta posterior
#       Assumes Dirichlet prior for base probability
#       Should sometimes include Ns (current method does not)
library(ape)
library(gtools) # rdirichlet
library(here)
#set.seed(2112) # \m/
source(here("covid", "aux_funk.R"))

args <- commandArgs(TRUE)
# --overwrite: Resample ALL
# -N number of resamples
# -d dirichlet prior
N <- 10
if ("-N" %in% args) {
    N <- args[which(args == "-N") + 1]
    N <- as.numeric(N)
}
print(args)
dirich <- "-d" %in% args

# Read list of files ending with .RDS
rds_names <- list.files("data/unc_covid/", pattern = "*RDS")
# The "-1" indicates that it's a copy. Remove it.
rds_names <- rds_names[!grepl("-1", rds_names)]

# Avoid re-tracing my steps
rds_done <- sapply(list.files("data/pangolineages/",
    pattern = ifelse(dirich, "*_d.csv", "*.csv")),
    function(x) strsplit(x, "\\_")[[1]][1])
rds_todo <- sapply(strsplit(rds_names, "-"),
    function(x) {
        strsplit(rev(x)[1], "\\.")[[1]][1]
    })

# Record accession names
asc_names <- parse_accession(rds_names)
if (FALSE) print(asc_names)

append <- FALSE
if (!"--overwrite" %in% args) {
    #rds_names <- rds_names[which(!rds_todo %in% rds_done)]
    append <- TRUE
}

# Prepare empty lists
S_list <- vector(mode = "list", length = length(rds_names))
header_list <- S_list

# Read in uncertainty matrices
for (i in seq_along(rds_names)) {
    S_list[[i]] <- readRDS(paste0("data/unc_covid/", rds_names[i]))
    names(S_list)[i] <- asc_names[i]
}

# Set up timing system
t0 <- Sys.time() # Outer timer
nloops <- length(S_list)
collapsed <- estlapseds <- double(nloops)
for (i in 1:nloops) {
    # Error checking

    t1 <- Sys.time() # Inner timer

    # Normalize matrix
    S <- S_list[[i]]
    S <- fix_unc(S)
    if (class(S) == "character") {
        print(paste(S, acc, sep = " - "))
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

    conseq <- apply(S, 1, function(x) {
        if (any(is.na(x)) | sum(x) < 10) {
            return("N")
        } else {
            return(alph[which.max(x)])
        }
    })

    # Convert sample letters to single string
    if (!append) {
        conseq <- paste(conseq, collapse = "", sep = "")
        sampleseq <- c(conseq, apply(sampleseq_mat, 1, paste, collapse = ""))

        # Create well-formatted fasta file
        name <- paste0("> ", asc_names[i], ".", 0:length(sampleseq))
        fasta <- paste(name, sampleseq, sep = "\n", collapse = "\n")
    } else {
        sampleseq <- apply(sampleseq_mat, 1, paste, collapse = "")
        name <- paste0("> ", asc_names[i], ".", 1:length(sampleseq))
        fasta <- paste(name, sampleseq, sep = "\n", collapse = "\n")
    }
    cat(fasta, append = append,
        file = paste0("data/sampled_covid/", asc_names[i],
            "_sampled", ifelse(dirich, "_d", ""), ".fasta"))

    # Loop Timing
    elapsed <- difftime(Sys.time(), t1, units = "mins")
    collapsed[i] <- elapsed
    estlapsed <- round((nloops - i) * mean(collapsed[collapsed > 0]), 3)
    estlapseds[i] <- estlapsed
    print(paste0("loop ", i, " of ", nloops, ", ", round(elapsed, 3),
        " mins | Approx. ", estlapsed, " mins (", round(estlapsed / 60, 2),
        " hours) remaining"))
}

#print(asc_names)
writeLines(asc_names, con = "data/sampled_covid/ascnames.txt")
