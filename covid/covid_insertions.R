# ---------------------------------------------------------
# Investigating / implementing insertion code -------------
# ---------------------------------------------------------

library(gtools) # rdirichlet
library(dplyr)

N <- 10
dirich <- TRUE

# ---------------------------------------------------------
# Get and exampled sample ---------------------------------
# ---------------------------------------------------------

S <- readRDS("samDONE/S-SRR12349131.RDS")

S <- fix_unc(S) # aux_funk.R
if ("character" %in% class(S)) {
    print(paste("S", asc_names[i], sep = " - "))
    next
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

sampleseq_mat[1, ]




# ---------------------------------------------------------
# Parsing the insertions output ---------------------------
# ---------------------------------------------------------
insertions <- readRDS("samDONE/S-SRR12349131_insertions.RDS")
str(insertions)

# Add a line that puts in their group membership and group size
insertions <- insertions %>%
    mutate(Position2 = Position) %>%
    arrange(Line.Number) %>%
    group_by(Line.Number) %>%
    arrange(Position) %>%
    mutate(group = cumsum(c(1, diff(Position) != 1))) %>%
    group_by(Line.Number, group) %>%
    mutate(Position2 = rep(min(Position), n()), size = n(), ins_id = 1:n()) %>%
    ungroup()

# Sanity Checks
table(insertions$size)
table(insertions$size) / as.numeric(names(table(insertions$size)))
insertions$Position2 %>% table() %>% sort(decreasing = TRUE)

posies <- unique(insertions$Position2)

out_list <- list()
for(i in seq_along(posies)) {
    this_pos <- as.data.frame(insertions[insertions$Position2 == posies[i],])
    coverage <- sum(S[this_pos$Position2[1], ])
    unique_lens <- sapply(unique(this_pos$Line.Number), function(x) sum(this_pos$Line.Number == x))
    ins_length_counts <- table(unique_lens)
    ins_length_coverage <- ins_length_counts / coverage
    ins_sample_probs <- c("0" = max(0, (coverage - sum(ins_length_counts)) / coverage), 
        ins_length_coverage)
    ins_sample_probs <- ins_sample_probs / sum(ins_sample_probs)

    S_list <- list()
    sizes <- unique(this_pos$size)
    if(length(sizes) > 1) print(i)
    for(ii in seq_along(sizes)) {
        this_size <- sizes[ii]
        ins_mat <- matrix(0, nrow = this_size, ncol = length(alph))
        colnames(ins_mat) <- alph
        this_pos_size <- this_pos[this_pos$size == this_size,]
        for(iii in seq_len(nrow(this_pos_size))) {
            rownum <- this_pos_size$ins_id[iii]
            colnum <- which(alph == trimws(this_pos_size$Base[iii]))
            qual <- this_pos_size$Error.Probability..1.p.[iii]
            if(is.na(qual)) {
                print(this_pos_size$Error.Probability..1.p.[iii])
                print(qual)
                qual <- 0
            }
            paired <- as.numeric(as.logical(trimws(this_pos_size$Paired[iii]))) + 1
            if(is.na(paired)) {
                print(this_pos_size$Paired[iii])
                paired <- 1
            }
            ins_mat[rownum, colnum] <- ins_mat[rownum, colnum] + qual / paired
        }
        S_list[[ii]] <- ins_mat
        names(S_list)[ii] <- this_size
    }

    out_list[[i]] <- list(Position2 = posies[i], S_list = S_list, 
        coverage = coverage, ins_sampler = ins_sample_probs)
    names(out_list)[i] <- posies[i]
}





# ---------------------------------------------------------
# Sample from the insertions code -------------------------
# ---------------------------------------------------------
# The sequences will be different lengths; list is better than matrix
sampleseq <- lapply(seq_len(nrow(sampleseq_mat)), function(x) sampleseq_mat[x, ])

for (i in seq_len(length(sampleseq))) { # One seq at a time
    insertions_to_add <- c()
    for (ii in seq_along(out_list)) {
        ins_sampler <- out_list[[ii]]$ins_sampler
        if (any(ins_sampler < 0)) {
            ins_sampler[ins_sampler < 0] <- 0
        }
        insertion_size <- sample(names(ins_sampler),
            size = 1, prob = ins_sampler)
        if (insertion_size == 0) {
            next
        } else {
            mini_s <- out_list[[ii]]$S_list[[insertion_size]]
            sampled_bases <- double(nrow(mini_s))
            for (iii in seq_len(nrow(mini_s))) {
                sampled_bases[iii] <- sample(colnames(mini_s), size = 1, prob = mini_s[iii,]) 
            }
            sampled_bases <- paste(sampled_bases, collapse = "", sep = "")
            names(sampled_bases) <- out_list[[ii]]$Position2
            insertions_to_add <- c(insertions_to_add, sampled_bases)
        }
    }

    insertions_reversed <- insertions_to_add[order(as.numeric(names(insertions_to_add)), decreasing = TRUE)]
    for (ii in seq_along(insertions_reversed)) {
        sampleseq[[i]] <- append(sampleseq[[i]], insertions_reversed[ii], after = as.numeric(names(insertions_reversed)[ii]))
        names(sampleseq[[i]]) <- NULL
    }
}





