
library(gtools) # rdirichlet

sample_S <- function(S, N = 10, dirich = TRUE){
    alph <- toupper(colnames(S))
    apply(S, 1, function(x) {
        if(any(is.na(x)) | (sum(x) < 10)){
            return(rep("N", N))
        } else {
            if (dirich){
                newx <- rdirichlet(1, x + c(rep(1/4, 4), rep(0, length(x) - 4)))
                return(sample(alph, size = N, prob = newx, replace = TRUE))
            } else {
                return(sample(alph, size = N, prob = x, replace = TRUE))
            }

        }
    })
}

conseq_S <- function(S){
    alph <- toupper(colnames(S))
    apply(S, 1, function(x) {
            if(any(is.na(x)) | sum(x) < 10) {
                return("N")
            } else {
                return(alph[which.max(x)])
        }
    })
}

calc_genlLik <- function(S, seq){
    if(length(seq) == 1) seq = strsplit(seq, "")
    if(length(seq) != nrow(S)) stop("Invalid Sequence")

    sum(sapply(1:nrow(S), function(x) {
        if(sum(S[x,]) > 10){
            return(log10(S[x, seq[x]] / sum(S[x, ])))
        } else {
            return(log10(1)) # TODO: Is there a better way?
        }
    }))
}

S <- readRDS("data/unc_covid/S-ERR4892386.RDS")
if(ncol(S) == 6) {
    S[,5] <- S[,5] + S[,6]
    S <- S[,1:5]
}

seqd1 <- sample_S(S, N=1000)
seqr1 <- sample_S(S, N=1000)
conseq1 <- conseq_S(S)

ld1 <- apply(seqd1, 1, calc_genlLik, S = S)
lr1 <- apply(seqr1, 1, calc_genlLik, S = S)
lcon <- calc_genlLik(S, conseq1)

mybreaks <- seq(min(ld1, lr1), max(ld1, lr1), length.out = 20)

hist(ld1, breaks = mybreaks, col = rgb(0,1,0,0.5))
hist(lr1, breaks = mybreaks, add = TRUE, 
    col = rgb(1,0,0,0.5))

ld1d <- density(ld1)
lr1d <- density(lr1)
ylime <- c(0, max(ld1d$y, lr1d$y))
xlime <- range(ld1d$x, lr1d$x)

plot(ld1d, ylim = ylime, xlim = xlime, col = "red")
lines(lr1d)














