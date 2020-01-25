library(dplyr)
library(ggplot2); theme_set(theme_bw())

set.seed(1234)


gen_seq_random <- function(seq.length, alphabet) {
  idx = sample(1:length(alphabet), size = seq.length, replace = TRUE)
  seq <- alphabet[idx]
  return(seq)
}


seq.length = 30
alphabet = c('A','C','G','T')
seq <- gen_seq_random(seq.length, alphabet)


gen_error_beta <- function(alpha, beta, seq.length, min.err = 1e-10) {
  err <- rbeta(n=seq.length, shape1 = alpha, shape2 = beta)
  err[err < min.err] <- min.err
  return(err)
}

alpha <- 0.1
beta  <- 20


err <- gen_error_beta(alpha, beta, seq.length)

err
hist(err)
summary(err)

translate_to_idx <- function(seq, alphabet) {
  idx <- numeric(seq.length)
  for(i in 1:seq.length){
    idx[i] <- which(alphabet == seq[i])
  }
  return(idx)
}

idx <- translate_to_idx(seq, alphabet)
idx
seq


fill_proba_unif <- function(v) {  # v = M[,1]
  x <- v[!is.na(v)]
  stopifnot(length(v) > length(x))
  y <- (1-sum(x)) / (length(v) - length(x))
  v[is.na(v)] <- y
  stopifnot(round(sum(v),2) == 1)
  return(v)
}

create_PS <- function(seq, err, alphabet) {
  
  M   <- matrix(nrow = length(alphabet), ncol = length(seq))
  idx <- translate_to_idx(seq, alphabet)
  
  for(j in 1:length(seq)){
    M[idx[j], j] <- 1 - err[j]
  }
  PS <- apply(M,MARGIN = 2, fill_proba_unif)
  rownames(PS) <- alphabet
  return(PS)
}

ps <- create_PS(seq, err, alphabet)


rmultinom2 <- function(prob, size, n) {
  return(rmultinom(n,size,prob))
}

z = rmultinom(n = 7, size = 1, prob = c(0.1,0.5,0.2,0.2))
w = apply(z, 2, function(x) { which(x>0)})
w

drawn.seq <- alphabet[w]
drawn.seq
# draw_from_PS

a = apply(ps, MARGIN = 2, FUN=rmultinom2, n=5, size = 1)


