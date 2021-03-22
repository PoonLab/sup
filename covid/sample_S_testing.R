# test sampleseq_mat
library(ape)
library(gtools) # rdirichlet
set.seed(2112) # \m/


N = 1000
dirich <- TRUE

# Read list of files ending with .RDS
rds_names <- list.files("data/unc_covid/", pattern = "*.RDS")
# The "-1" indicates that it's a copy. Remove it.
rds_names <- rds_names[!grepl("-1", rds_names)]


# Record accession names
asc_names <- c()
for(i in 1:length(rds_names)){
    asc <- strsplit(rev(strsplit(rds_names[i], "-")[[1]])[1], "\\.")[[1]][1]
    asc_names <- c(asc_names, sub(".RDS", "", asc))
}
print(asc_names)


# Prepare empty lists
S_list <- vector(mode = "list", length = length(rds_names))
header_list <- S_list
#print(rds_names)

# Read in uncertainty matrices 
for(i in 1:length(rds_names)){
    S_list[[i]] <- readRDS(paste0("data/unc_covid/", 
        rds_names[i]))
    names(S_list)[i] <- asc_names[i]
}

# Check coverage of each read
#hist(apply(S_list[[1]], 1, sum), breaks = 30)


i <- 1
is.null(dim(S_list[[i]]))

S <- S_list[[i]]
if(ncol(S) == 6) {
    S[,5] <- S[,5] + S[,6]
    S <- S[,1:5]
}
alph <- toupper(colnames(S))


N <- 1000
dirich <- !FALSE
sampleseq_mat <- apply(S, 1, function(x) {
    if(any(is.na(x)) | (sum(x) < 10)){
        return(rep("N", N))
    } else {
        if (dirich){
            newx <- rdirichlet(1, x + rep(1/4, length(x)))
            return(sample(alph, size = N, prob = newx,
                replace = TRUE))
        } else {
            return(sample(alph, size = N, prob = x, 
                replace = TRUE))
        }

    }
})

dim(sampleseq_mat)

# Manually re-ran, iterating the name
# It was quick and dirty and I refuse to be judged.
mytab_d4 <- table(apply(sampleseq_mat, 2, 
    function(x) length(unique(x))))

dmat <- rbind(mytab_r1, mytab_r2, mytab_r3, mytab_r4)
rmat <- rbind(mytab_d1, mytab_d2, mytab_d3, mytab_d4)
dmat <- cbind(dmat, rep(0, 4))

df <- data.frame(rbind(dmat, rmat))
names(df) <- 1:5
df$dirich <- c("rsample", "dirich")[grepl(
    row.names(df), pattern = "d") + 1]

library(ggplot2)
library(tidyr)
ggplot(pivot_longer(df, cols = -dirich), 
        aes(x = factor(name), y = value, fill = dirich)) + 
    geom_col(position = "dodge")
# So yes, this increases the variance in the base calls, as expected.

apply(as.data.frame(df[1:5]), 1, function(x) x[1] / sum(x[2:4]))
# Random sampling: ~40% == 1
# Dirichlet: ~25% == 1
# Dirichlet successfully adds variance to the genomes.
