# open S matrix and sample from it
library(ape)
library(gtools) # rdirichlet
library(here)
#set.seed(2112) # \m/
source(here("covid", "aux_funk.R"))
library(dplyr)

args <- commandArgs(TRUE)
# --overwrite: Resample ALL
# -N number of resamples
# -d dirichlet prior
N <- 10
if ("-N" %in% args) {
    N <- args[which(args == "-N") + 1]
    N <- as.numeric(N)
}

append <- FALSE
if (!"--overwrite" %in% args) {
    append <- TRUE
}

dirich <- "-d" %in% args

print(args)

# Read list of files ending with .RDS
rds_names <- list.files("data/unc_covid/", pattern = "*RDS")
rds_names <- rds_names[!grepl(rds_names, pattern = "_insertions")]

# Record accession names
asc_names <- unique(parse_accession(rds_names)) # aux_funk.R
print(asc_names)

# Prepare empty lists
S_list <- vector(mode = "list", length = length(rds_names))
header_list <- S_list

# Read in uncertainty matrices
for (i in seq_along(rds_names)) {
    S_list[[i]] <- readRDS(paste0("data/unc_covid/", rds_names[i]))
    names(S_list)[i] <- asc_names[i]
}


# Set up timing system --------------------------
t0 <- Sys.time() # Outer timer
nloops <- length(S_list)
collapsed <- estlapseds <- double(nloops)
for (i in 1:nloops) {
    t1 <- Sys.time() # Inner timer
    print(asc_names[i])

    # Normalize matrix
    S <- S_list[[i]]
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
                    as.numeric(x)
                return(sample(alph, size = N, prob = newx, replace = TRUE))
            } else {
                return(sample(alph, size = N, prob = x, replace = TRUE))
            }

        }
    })

    # ---------------------------------------------------------
    # Parsing the insertions output ---------------------------
    # ---------------------------------------------------------
    insertion_name <- paste0("S-", asc_names[i], "_insertions.RDS")
    if(insertion_name %in% list.files(here("data/unc_covid"))) {
        insertions <- readRDS(here("data/unc_covid", insertion_name))

        # Add a line that puts in their group membership and group size
        insertions <- insertions %>%
            filter(!is.na(Position)) %>%
            mutate(Position2 = Position) %>%
            arrange(Line.Number) %>%
            group_by(Line.Number) %>%
            arrange(Position) %>%
            mutate(group = cumsum(c(1, diff(Position) != 1))) %>%
            group_by(Line.Number, group) %>%
            mutate(Position2 = rep(min(Position), n()), size = n(), ins_id = 1:n()) %>%
            ungroup() %>%
            filter(!is.na(Line.Number), !is.na(Position2)) %>%
            mutate(Cigar = trimws(Cigar), Base = trimws(Base), Paired = trimws(Paired))

        posies <- unique(insertions$Position2)

        out_list <- list()
        for(ii in seq_along(posies)) {
            this_pos <- as.data.frame(insertions[insertions$Position2 == posies[ii],])
            coverage <- sum(S[this_pos$Position2[1], ])
            unique_lens <- sapply(unique(this_pos$Line.Number), function(x) sum(this_pos$Line.Number == x))
            ins_length_counts <- table(unique_lens)
            ins_length_coverage <- ins_length_counts / coverage
            ins_sample_probs <- c("0" = max(0, (coverage - sum(ins_length_counts)) / coverage), 
                ins_length_coverage)
            ins_sample_probs <- ins_sample_probs / sum(ins_sample_probs)

            ins_list <- list()
            sizes <- unique(this_pos$size)
            #if(length(sizes) > 1) print(ii)
            for(iii in seq_along(sizes)) {
                this_size <- sizes[iii]
                ins_mat <- matrix(0, nrow = this_size, ncol = length(alph))
                colnames(ins_mat) <- alph
                this_pos_size <- this_pos[this_pos$size == this_size,]
                for(iiii in seq_len(nrow(this_pos_size))) {
                    rownum <- this_pos_size$ins_id[iiii]
                    colnum <- which(alph == trimws(this_pos_size$Base[iiii]))
                    qual <- this_pos_size$Error.Probability..1.p.[iiii]
                    if(is.na(qual)) {
                        print(this_pos_size$Error.Probability..1.p.[iiii])
                        print(qual)
                        qual <- 0
                    }
                    paired <- as.numeric(as.logical(trimws(this_pos_size$Paired[iiii]))) + 1
                    if(is.na(paired)) {
                        print(this_pos_size$Paired[iiii])
                        paired <- 1
                    }
                    ins_mat[rownum, colnum] <- ins_mat[rownum, colnum] + qual / paired
                }
                ins_list[[iii]] <- ins_mat
                names(ins_list)[iii] <- this_size
            }

            out_list[[ii]] <- list(Position2 = posies[ii], S_list = ins_list, 
                coverage = coverage, ins_sampler = ins_sample_probs)
            names(out_list)[ii] <- posies[ii]
        }


        # ---------------------------------------------------------
        # Sample from the insertions code -------------------------
        # ---------------------------------------------------------
        # The sequences will be different lengths; list is better than matrix
        sampleseq <- lapply(seq_len(nrow(sampleseq_mat)), function(x) sampleseq_mat[x, ])

        for (ii in seq_len(length(sampleseq))) { # One seq at a time
            insertions_to_add <- c()
            for (iii in seq_along(out_list)) {
                ins_sampler <- out_list[[iii]]$ins_sampler
                if (any(is.na(ins_sampler))) {
                    next
                }
                if (any(ins_sampler < 0)) {
                    ins_sampler[ins_sampler < 0] <- 0
                }
                insertion_size <- sample(names(ins_sampler),
                    size = 1, prob = ins_sampler)
                if (insertion_size == 0) {
                    next
                } else {
                    mini_s <- out_list[[iii]]$S_list[[insertion_size]]
		    if (length(as.numeric(mini_s) == 0)) { next }
                    if (length(is.na(mini_s)) == 0) { next }
                    if (is.na(mini_s)) { next }
	            if (is.na(nrow(mini_s))) { next } 
                    if (nrow(mini_s) == 0) { next }
                    sampled_bases <- double(nrow(mini_s))
                    for (iiii in seq_len(nrow(mini_s))) {
                        sampled_bases[iiii] <- sample(colnames(mini_s), size = 1, prob = mini_s[iiii,]) 
                    }
                    sampled_bases <- paste(sampled_bases, collapse = "", sep = "")
                    names(sampled_bases) <- out_list[[iii]]$Position2
                    insertions_to_add <- c(insertions_to_add, sampled_bases)
                }
            }

            insertions_reversed <- insertions_to_add[order(as.numeric(names(insertions_to_add)), decreasing = TRUE)]
            for (iii in seq_along(insertions_reversed)) {
                sampleseq[[ii]] <- append(sampleseq[[ii]], insertions_reversed[ii], after = as.numeric(names(insertions_reversed)[iii]))
                names(sampleseq[[ii]]) <- NULL
            }
        }
    } else {
        sampleseq <- lapply(seq_len(nrow(sampleseq_mat)), function(x) sampleseq_mat[x, ])
    }


    conseq <- apply(S, 1, function(x) {
        if (any(is.na(x)) | sum(x) < 10) {
            return("N")
        } else {
            return(alph[which.max(x)])
        }
    })

    # Convert sample letters to single string
    if (!append) {
        conseq <- paste(conseq, collapse = "", sep = "")
        sampleseq <- sapply(sampleseq, paste, collapse = "")
        sampleseq <- c(conseq, sampleseq)

        # Create well-formatted fasta file
        name <- paste0("> ", asc_names[i], ".", 0:(length(sampleseq) - 1))
        fasta <- paste(name, sampleseq, sep = "\n", collapse = "\n")

    } else {
        # When appending, conseq is not needed.
        sampleseq <- sapply(sampleseq, paste, collapse = "")
        name <- paste0("> ", asc_names[i], ".", 1:length(sampleseq))
        fasta <- paste(name, sampleseq, sep = "\n", collapse = "\n")
    }
    cat(fasta, append = append,
        file = paste0("data/sampled_covid/", asc_names[i],
            "_sampled", ifelse(dirich, "_d", ""), ".fasta"))

    # Loop Timing
    elapsed <- difftime(Sys.time(), t1, units = "mins")
    collapsed[i] <- elapsed
    estlapsed <- round((nloops - i) * mean(collapsed[collapsed > 0]), 3)
    estlapseds[i] <- estlapsed
    print(paste0("loop ", i, " of ", nloops, ", ", round(elapsed, 3),
        " mins | Approx. ", estlapsed, " mins (", round(estlapsed / 60, 2),
        " hours) remaining"))
}

#print(asc_names)
writeLines(asc_names, con = "data/sampled_covid/ascnames.txt")
