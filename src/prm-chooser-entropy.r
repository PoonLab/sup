# Choose parameters based on equally spaced entropy
# NOTE: Entropy from integrating beta distribution
# is different from calculating for data

beta_entropy_aux <- function(a, b, xdelta = 0.001) {
    # Riemann sum, midpoint rule to avoid inf
    xseq <- seq(xdelta / 2, 1 - xdelta / 2, xdelta)
    beta_dist <- dbeta(xseq, a, b)
    sum(beta_dist * log2(beta_dist) * xdelta)
}

mom <- function(xbar, s) {
    tmp <- xbar * (1 - xbar) / s - 1
    alpha <- xbar * tmp
    beta <- (1 - xbar) * tmp
    return(c(xbar = xbar, s = s, alpha = alpha, beta = beta))
}

thresh_crosser <- function(x, thresh) {
    if (all(x[!is.na(x)] < thresh)) return(length(x) + 1)
    start <- max(which(x <= thresh))
    start + min(which(x[start:length(x)] >= thresh))
}


beta_entropy <- function(beta_sd = 0.001) {
    beta_xbar <- seq(0.0001, 0.06, 0.0001)
    prmset <- t(sapply(beta_xbar, mom, s = beta_sd))
    prmset <- as.data.frame(prmset)
    prmset$entropy <- NA

    for (i in seq_len(nrow(prmset))) {
        a <- prmset$alpha[i]
        b <- prmset$beta[i]
        # skip invalid values
        if (a <= 0 | b <= 0) next

        prmset$entropy[i] <- beta_entropy_aux(a, b)
    }
    prmset
}

par(mfrow = c(3, 1))
for (beta_sd in c(0.0001, 0.001, 0.01)) {
    prmset <- beta_entropy(beta_sd)
    prmset$oom <- sapply(prmset$xbar, function(x) {
        min(which(10^(-(0:7)) <= x)) - 1
    })
    table(prmset$oom)
    plot(entropy ~ xbar, data = prmset)
    rug(prmset$xbar[is.na(prmset$entropy)])
    abline(h = 0:7)
    for (ent in 0:6) {
        print(prmset[thresh_crosser(prmset$entropy,
            ent), ])
    }
}