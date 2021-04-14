library(dplyr)
library(tidyr)
library(here)

# Arg Parsing -----------------------------------
args <- commandArgs(trailingOnly = TRUE)
print(args)

# Finds the value after -N, otherwise N = 1000
N <- 1000
if ("-N" %in% args) {
    N <- as.numeric(args[which(args == "-N") + 1])
}

overwrite <- FALSE
if ("--overwrite" %in% args) {
    overwrite <- TRUE
}

# Read in files (using here::here()) ------------
in_path <- "data/unc_covid/"
out_path <- "data/ord_covid"
in_files <- list.files(here(in_path), pattern = "RDS")


# Loop though files
for (i in seq_along(in_files)) {
    # Smallest number of sequences to ensure an output of N
    n <- min(which(2 ^ (1:15) >= N))

    # Set up filenames, skip if --overwrite flag is present
    in_file <- in_files[i]

    # Ugly regex/string splitting to get accession number
    # Files are "blah-blah-blah-S-[accession number].RDS"
    acc <- rev(
        strsplit(
            strsplit(
                in_file, "\\."
                )[[1]][1], "-"
            )[[1]]
        )[1]
    out_file <- here(out_path, paste0(acc, "_ord.fasta"))
    if ((!overwrite) & 
            paste0(acc, "_ord.fasta") %in% list.files(out_path)) {
        print(paste0(in_file, " exists, skipping."))
        next
    }
    print(in_file)

    # File processing (TODO: make this a function)
    unc_mat <- readRDS(here(in_path, in_file))
    if (is.null(dim(unc_mat))) {
        print("Empty file")
        next
    }
    if ("X" %in% colnames(unc_mat)) {
        unc_mat <- unc_mat[, -which(colnames(unc_mat) == "X")]
    }
    if (ncol(unc_mat) == 6) {
        unc_mat[, 5] <- unc_mat[, 5] + unc_mat[, 6]
        unc_mat <- unc_mat[, 1:5]
    }
    if (any(unc_mat[!is.na(unc_mat)] < 0 | 
            unc_mat[!is.na(unc_mat)] > 10e8)) {
        print(paste0("Values too small or too large - ", in_file))
        next
    }



    bottom_n <- function(S, n) {
        # Find the n substitutions with smallest differences from the conseq
        # S: uncertainty matrix
        # n: number of substitution sites to check
        M_all <- apply(S, 1, max)
        M_all[M_all < 10] <- Inf
        bottom <- which(rank(M_all) <= n)
        # If there aren't n uncertainties, find the number of uncertainties
        if (length(bottom) < n | any(is.na(bottom))) {
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
        for (i in seq_along(bottom)) {
            dummy_list[[i]] <- c(1, 2)
        }

        subs <- expand.grid(dummy_list)
        colnames(subs) <- bottom
        subs$diff_lik <- NA

        for (j in seq_len(nrow(subs))) {
            new_subs <- which(subs[j, 1:n] == 2)
            new_lik <- M
            new_lik[new_subs] <- m[new_subs]
            new_lik <- new_lik/M
            subs$diff_lik[j] <- sum(log(new_lik))
        }
        subs
    }

    # Calculate likelihoods ---------------------
    lik_mat <- bottom_n(unc_mat, n)
    # put the sequences in order of the likelihood
    lik_mat <- lik_mat[order(-lik_mat$diff_lik), ]
    if (nrow(lik_mat) < N) {
        print(paste0(acc, " failed - not enough substitutions"))
        next
    }
    lik_mat <- lik_mat[1:N, ]

    # now:
    # Calculate conseq
    # make the top 1000 substitutions
    # Run through pangolin
    alph <- colnames(unc_mat)
    conseq <- apply(unc_mat, 1, function(x) {
        if (any(is.na(x)) | sum(x) < 1) {
            return("N")
        } else {
            return(alph[which.max(x)])
        }
    })

    # Second most likely base call at each locus.
    # Necessary to make the substitutions.
    runner_up <- apply(unc_mat, 1, function(x) {
        if (any(is.na(x)) | sum(x) < 1) {
            return("N")
        } else {
            return(alph[order(-x)[2]])
        }
    })

    # Prep empty list of sequences
    ordered_seq <- lapply(1:min(N, nrow(lik_mat)), 
        function(x) conseq)

    # Switch the conseq call with the relevant substitution
    for (j in seq_len(nrow(lik_mat))) {
        switch <- lik_mat[j, 1:(ncol(lik_mat) - 1)]
        to_switch <- colnames(switch)[as.numeric(switch) == 2]
        to_switch <- as.numeric(to_switch)
        ordered_seq[[j]][to_switch] <- runner_up[to_switch]
    }

    # Begin Fasta prep --------------------------
    conseq <- paste(conseq, collapse = "")
    ordered_seq <- sapply(ordered_seq, paste, collapse = "")

    accs <- rep(acc, length(ordered_seq))
    weights <- round(lik_mat$diff_lik[seq_along(ordered_seq)], 6)
    subs <- apply(lik_mat[, 1:n], 1, paste0, collapse = "")

    fasta_labels <- paste(">",
        paste(accs, weights, subs, sep = "_")
    )

    to_save <- paste(fasta_labels, c(conseq, ordered_seq),
        collapse = "\n", sep = "\n"
    )

    writeLines(to_save, out_file)
}
