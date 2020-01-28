library(dplyr)
library(ggplot2); theme_set(theme_bw())

library(seqinr)

set.seed(1234)

#' Generate a sequnce randomly (for tests & debugging)
#' @param seq.length Integer. Sequence length.
#' @param alphabet Character vector. Alphabet from which the sequence is built.
#' @return Character vector representing the molecular sequence. 
#' 
gen_seq_random <- function(seq.length, alphabet) {
  idx = sample(1:length(alphabet), size = seq.length, replace = TRUE)
  seq <- alphabet[idx]
  return(seq)
}


#' Generate observation-error probabilities 
#' for each position of a sequence.
#' Probabilities are drawn from a Beta distribution. 
#' @param alpha Numeric. First shape parameter of the Beta distribution.
#' @param beta Numeric. Second shape parameter of the Beta distribution.
#' @param seq.length Integer. Sequence length.
#' @param min.err Numeric. Floor for probability values. Default = 1e-10.
#' @return Numeric vector of observation-error probabilities of length `seq.length`.
#' 
gen_error_beta <- function(alpha, beta, seq.length, min.err = 1e-10) {
  err <- rbeta(n=seq.length, shape1 = alpha, shape2 = beta)
  err[err < min.err] <- min.err
  return(err)
}


#' Translate an alphabet element into its position in the alphabet.
#' For example, if alphabet = c('A','C','G','T') then 
#' A --> 1, C --> 2, etc. 
#' @param seq Character vector. Sequence. 
#' @param alphabet Character vector. Alphabet definition.
#' @return Integer vector of indices.
#' 
translate_to_idx <- function(seq, alphabet) {
  idx <- numeric(seq.length)
  for(i in 1:seq.length){
    idx[i] <- which(alphabet == seq[i])
  }
  return(idx)
}


#' Uniform filling. 
#' For a given vector with at least one of its element
#' a probability and the rest of its elements NA, 
#' fills the NA elements with equal probability values.
#' For example, v = c(NA, NA, 0.4, NA) --> c(0.2, 0.2, 0.4, 0.2)
#' @param v Numeric vector of NAs and probability values.
#' @return Numeric vector of probability values.
#' 
fill_proba_unif <- function(v) {  # v = M[,1]
  x <- v[!is.na(v)]
  stopifnot(length(v) > length(x))
  y <- (1-sum(x)) / (length(v) - length(x))
  v[is.na(v)] <- y
  stopifnot(round(sum(v),2) == 1)
  return(v)
}


#' Create a probability sequence. Uniform filling.
#' @param seq Character vector. Observed sequence. 
#' @param err Numeric vector. Observation-error probabilities for every position.
#' @param alphabet Character vector. Alphabet definition. 
#' @return Matrix of observation probabilities. Number of rows = alphabet length, number of columns = sequence length.
create_PS_unif <- function(seq, err, alphabet) {
  M   <- matrix(nrow = length(alphabet), ncol = length(seq))
  idx <- translate_to_idx(seq, alphabet)
  # Observed nucleotide:
  for(j in 1:length(seq)){
    M[idx[j], j] <- 1 - err[j]
  }
  # Other nucleotides:
  PS <- apply(M,MARGIN = 2, fill_proba_unif)
  rownames(PS) <- alphabet
  return(PS)
}


#' Generate `n` letters based on an alphabet probabilities.
#' @param prob Numeric vector. Probabilities of observation for every letter of the alphabet.
#' @param n Integer. Number of replicates.
#' @return Character vector of the letters generated.
#' 
gen_letters <- function(prob, n) { # n = 9
  # Only one nucleotide can be drawn 
  # from the probability sequence:
  x <- rmultinom(n, size=1, prob)
  
  # which nucleotide was drawn:
  y <- apply(x, 2, function(x) {which(x>0)})
  
  # Translate back into the alphabet:
  z <- alphabet[y]
  return(z)
}


#' Draw sequence replicates from an observed sequence
#' and a Beta/unif sequence uncertainty model.  
#' @param obs.seq Character vector. Observed sequence.
#' @param prm.beta List. Parameter (alpha,beta) of the Beta distribution that generates the observation-error probabilities.
#' @param n.repl Integer. Number of replicates drawn. 
#' @param alphabet Character vector. Alphabet definition.
#' @return List : `seq.drawn`:  `n.repl` sequences drawn from the sequence-uncertainty model.
#' `obs.err`: observation error for every position of the observed sequence. 
#' 
draw_sequences_beta_unif <- function(obs.seq, 
                                prm.beta, 
                                n.repl, 
                                alphabet) {
  # Draw errors:
  alpha <- prm.beta[["alpha"]]
  beta  <- prm.beta[["beta"]]
  obs.err <- gen_error_beta(alpha, beta, seq.length = length(obs.seq))
  
  # Calculate the Probability Sequence:
  ps <- create_PS_unif(obs.seq, obs.err, alphabet)
  
  # Draw `n.repl` replicates from the Probability Sequence:
  seq.drawn  <- apply(ps, MARGIN = 2, FUN=gen_letters, n=n.repl)
  
  # Convert matrix to list:
  seq.drawn.list <- lapply(seq_len(nrow(seq.drawn)), 
                           function(i) seq.drawn[i,])
  
  return(list(seq.drawn = seq.drawn.list, 
              obs.err   = obs.err))
}

# ---- RUN ----

seq.length = 1000
alphabet = c('A','C','G','T')
obs.seq <- gen_seq_random(seq.length, alphabet)
paste(obs.seq[1:300], collapse = '')

prm.beta <- list(alpha = 0.1, 
                 beta  = 30)
n.repl <- 100

system.time({
res <- draw_sequences_beta_unif(obs.seq, 
                           prm.beta, 
                           n.repl, 
                           alphabet)
})

seq.drawn <- res$seq.drawn
obs.err   <- res$obs.err

h <- hist(obs.err, breaks = 30, plot = F)
plot(h$mids, h$density, log='y', typ='o', las=1, pch=16, cex=1, lwd=2,
     main = paste("Observation-Error Probability Distribution\n",
                  "Beta Unif Model:", paste(prm.beta,collapse = ' ; ')),
     xlab = 'Obs. Err. Probability', ylab = 'Density')
grid()

#M <- matrix(unlist(seq.drawn), nrow = n.repl, byrow = TRUE)

