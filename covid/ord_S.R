library(dplyr)
library(tidyr)
library(here)

args <- commandArgs(trailingOnly = TRUE)

if ("-N" %in% args) {
    N <- as.numeric(args[which(args == "-N") + 1])
} else {
    N <- 1000
}

in_path <- "data/unc_covid/"
out_path <- "data/ord_covid"
in_files <- list.files(here(in_path), pattern = "RDS")

for (i in seq_along(in_files)) {
    n <- min(which(2 ^ (1:15) >= N))
    in_file <- in_files[i]
    print(in_file)

    sam1 <- readRDS(here(in_path, in_file))
    if (is.null(dim(sam1))) {
        print("Empty file")
        next
    }
    if ("X" %in% colnames(sam1)) {
        sam1 <- sam1[, -which(colnames(sam1) == "X")]
    }
    if (ncol(sam1) == 6) {
        sam1[, 5] <- sam1[, 5] + sam1[, 6]
        S <- sam1[, 1:5]
    }
    if (any(sam1[!is.na(sam1)] < 0 | sam1[!is.na(sam1)] > 10e8)) {
        print(paste0("Values too small or too large - ", in_file))
        next
    }

    alph <- toupper(colnames(sam1))

    acc <- rev(strsplit(strsplit(in_file, "\\.")[[1]][1], "-")[[1]])[1]

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

    bottom_n <- function(S, n) {
        M_all <- apply(S, 1, max)
        bottom <- which(rank(M_all) <= n)
        # If there aren't n uncertainties, just find the number of uncertainties
        if (length(bottom) < n) {
            print(paste0("Warning: ", acc, " did not have ",
                n, " uncertain bases"))
            bottom <- bottom[!is.na(bottom)]
            n <- length(bottom)
        }

        M <- M_all[bottom]
        m <- apply(S[bottom, ], 1, function(x) {
            sort(x, decreasing = TRUE)[2]
        })

        dummy_list <- list()
        for (i in seq_along(bottom[!is.na(bottom)])) {
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

    lik_mat <- bottom_n(sam2, n)
    lik_mat <- lik_mat[order(-lik_mat$diff_lik), ]
    # now:
    # Calculate conseq
    # make the top 1000 substitutions
    # Run through pangolin
    alph <- colnames(sam1)
    conseq <- apply(sam1, 1, function(x) {
        if (any(is.na(x)) | sum(x) < 1) {
            return("N")
        } else {
            return(alph[which.max(x)])
        }
    })

    runner_up <- apply(sam1, 1, function(x) {
        if (any(is.na(x)) | sum(x) < 1) {
            return("N")
        } else {
            return(alph[order(-x)[2]])
        }
    })

    ordered_seq <- lapply(1:min(1000, nrow(lik_mat)), function(x) conseq)
    for (j in seq_along(ordered_seq)) {
        switch <- lik_mat[j, 1:(ncol(lik_mat) - 1)]
        to_switch <- colnames(switch)[as.numeric(switch) == 2]
        to_switch <- as.numeric(to_switch)
        ordered_seq[[j]][to_switch] <- runner_up[to_switch]
    }

    conseq <- paste(conseq, collapse = "")
    ordered_seq <- sapply(ordered_seq, paste, collapse = "")

    fasta_labels <- paste("> ", acc,
        c(0, round(lik_mat$diff_lik[seq_along(ordered_seq)], 6)),
        sep = ""
    )

    to_save <- paste(fasta_labels, c(conseq, ordered_seq),
        collapse = "\n", sep = "\n"
    )

    writeLines(to_save, here(out_path, paste0(acc, "_ord.fasta")))
}
