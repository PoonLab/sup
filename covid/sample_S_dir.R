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

args <- commandArgs(TRUE)
# --overwrite: Resample ALL
# -N number of resamples
# -d dirichlet prior
N = 10
if("-N" %in% args){
    N = args[which(args == "-N") + 1]
    N = as.numeric(N)
} 
dirich <- "-d" %in% args

# Read list of files ending with .RDS
rds_names <- list.files("data/unc_covid/", pattern = "*.RDS")
# The "-1" indicates that it's a copy. Remove it.
rds_names <- rds_names[!grepl("-1", rds_names)]

# Avoid re-tracing my steps
rds_done <- sapply(list.files("data/pangolineages/", pattern=ifelse(dirich, "*_d.csv", "*.csv")), 
    function(x) strsplit(x, "\\_")[[1]][1])
rds_todo <- sapply(strsplit(rds_names, "-"),
    function(x){
        strsplit(rev(x)[1], "\\.")[[1]][1]
    })

# Record accession names
asc_names <- c()
for(i in 1:length(rds_names)){
    asc <- strsplit(rev(strsplit(rds_names[i], "-")[[1]])[1], "\\.")[[1]][1]
    asc_names <- c(asc_names, sub(".RDS", "", asc))
}
print(asc_names)

if(!"--overwrite" %in% args) rds_names <- rds_names[which(!rds_todo %in% rds_done)]

# Prepare empty lists
S_list <- vector(mode = "list", length = length(rds_names))
header_list <- S_list
#print(rds_names)

# Read in uncertainty matrices 
for(i in 1:length(rds_names)){
    S_list[[i]] <- readRDS(paste0("data/unc_covid/", rds_names[i]))
    names(S_list)[i] <- asc_names[i]
}

# Check coverage of each read
#hist(apply(S_list[[1]], 1, sum), breaks = 30)

# Prepare empty list
sampled_files <- list()

# Set up timing system
t0 <- Sys.time() # Outer timer
nloops <- length(S_list)
collapsed <- estlapseds <- double(nloops)
for (i in 1:2){#nloops) {
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

    conseq <- apply(S, 1, function(x) {
        if(any(is.na(x)) | sum(x) < 10) {
            return("N")
        } else {
            return(alph[which.max(x)])
        }
    })

    # Convert sample letters to single string
    conseq <- paste0(conseq, collapse = "")
    name <- paste0("> ", asc_names[i], ".", 0, 
        collapse = "")

    filename <- paste0("data/sampled_covid/", 
        asc_names[i], "_sampled", 
        ifelse(dirich, "_d", ""), ".fasta", 
        collapse = "")
    system(paste0("echo '", name, "' > ", filename, 
        collapse = ""))
    system(paste0("echo '", conseq, "' >> ", filename, 
        collapse = ""))

    for (ii in 1:N) {
        sampleseq <- apply(S, 1, function(x) {
            if(any(is.na(x)) | (sum(x) < 10)){
                return("N")
            } else {
                if (dirich){
                    newx <- rdirichlet(1, x + 
                        c(rep(1/4, 4), 
                            rep(0, length(x) - 4)))
                    return(sample(alph, size = 1,
                        prob = newx, replace = TRUE))
                } else {
                    return(sample(alph, size = 1, 
                        prob = x, replace = TRUE))
                }
            }
        })

        sampleseq <- paste0(sampleseq, collapse = "")
        name <- paste0("> ", asc_names[i], ".", 
            ii, collapse = "")

        system(paste0("echo '", name, "' >> ", filename, 
            collapse = ""))
        system(paste0("echo '", sampleseq, "' >> ",
            filename,  collapse = ""))

    }

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













