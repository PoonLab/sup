library(dplyr)
library(tidyr)
library(here)

in_file <- "open_seSAMe-S-ERR4363387.RDS"
in_path <- "data/unc_covid/"
sam1 <- readRDS(here(in_path, in_file))

dim(sam1)

if (ncol(sam1) == 6) {
    sam1[, 5] <- sam1[, 5] + sam1[, 6]
    sam1 <- sam1[, 1:5]
}

sam2 <- apply(sam1, 1, function(x) {
    if (sum(x) > 0 & !all(is.na(x))) {
        x_norm <- x / sum(x, na.rm = TRUE)
        if (any(x_norm > 0.99)) {
            newx <- rep(0, length(x_norm))
            newx[which(x_norm > 0.99)] <- 1
            return(newx)
        } else {
            return(x_norm)
        }
    } else {
        return(rep(NA, ncol(sam1)))
    }
})

sam2 <- t(sam2)
dim(sam1) == dim(sam2)

to_sample <- apply(sam2, 1, function(x) {
    (max(x) < 1) & (!is.na(max(x)))
})
table(to_sample)

non0 <- sam2[to_sample, ] %>%
    `colnames<-`(1:5) %>%
    as.data.frame() %>%
    tibble::rowid_to_column(var = "row_id") %>%
    tidyr::pivot_longer(cols = -row_id, names_to = "col_id") %>%
    mutate(order = order(value))

# Next steps:
    # 1. Expand grid to include labels of permutations
    # 2. Use this grid to lookup M and m
    # 3. Multiply!

bottom_n <- function(S, n) {
    M_all <- apply(S, 1, max)
    bottom <- which(rank(M_all) <= n)

    M <- M_all[bottom]
    m <- apply(S[bottom, ], 1, function(x) {
        sort(x, decreasing = TRUE)[2]
    })

    dummy_list <- list()
    for (i in seq_len(n)) {
        dummy_list[[i]] <- c(1, 2)
    }

    subs <- expand.grid(dummy_list)
    colnames(subs) <- bottom
    subs$diff_lik <- NA

    for (i in seq_len(nrow(subs))) {
        new_subs <- which(subs[i, 1:n] == 2)
        new_lik <- M
        new_lik[new_subs] <- m[new_subs]
        subs$diff_lik[i] <- sum(log(new_lik))
    }
    subs
}

# Testing
lik_mat <- bottom_n(sam2, 12)
lik_mat <- lik_mat[order(-lik_mat$diff_lik), ]
head(tmp)
nrow(tmp)

# now:
    # Calculate conseq
    # make the top 1000 substitutions
    # Run through pangolin
alph <- colnames(sam1)
conseq <- apply(S, 1, function(x) {
    if (any(is.na(x)) | sum(x) < 1) {
        return("N")
    } else {
        return(alph[which.max(x)])
    }
})

runner_up <- apply(S, 1, function(x) {
    if (any(is.na(x)) | sum(x) < 1) {
        return("N")
    } else {
        return(alph[order(-x)[2]])
    }
})

ordered_seq <- lapply(1:1000, function(x) conseq)
for (i in seq_along(ordered_seq)) {
    switch <- lik_mat[i, 1:(ncol(lik_mat) - 1)]
    to_switch <- colnames(switch)[as.numeric(switch) == 2]
    to_switch <- as.numeric(to_switch)
    ordered_seq[[i]][to_switch] <- runner_up[to_switch]
}

conseq <- paste(conseq, collapse = "")
ordered_seq <- sapply(ordered_seq, paste, collapse = "")

acc <- rev(strsplit(strsplit(in_file, "\\.")[[1]][1], "-")[[1]])[1]
fasta_labels <- paste("> ", acc,
        c(0, round(lik_mat$diff_lik, 6)),
    sep = "")

to_save <- paste(fasta_labels, c(conseq, ordered_seq),
    collapse = "\n", sep = "\n")

writeLines(to_save, here("data/ord_covid", "test1.fasta"))
