require(parallel)

# adapted from http://github.com/PoonLab/MiCall-Lite

#' re.findall
#'
#' Emulate behaviour of Python's re.findall() function
#'
#' @param pat:  regex pattern
#' @param s:  character, a single string
#' @return character, vector of all matching substrings
re.findall <- function(pat, s) {
    if (!is.character(s)) {
        stop("re.findall() requires a character object for input 's'")
    }
    matches <- gregexpr(pat, s)
    index <- as.integer(matches[[1]])
    match.length <- attr(matches[[1]], 'match.length')
    
    sapply(1:length(index), function(i) {
        start <- index[i]
        stop <- start + match.length[i] - 1
        substr(s, start, stop)
    })
}


#' Use CIGAR (Compact Idiosyncratic Gapped Alignment Report) string
#' to apply soft clips, insertions and deletions to the read sequence.
#' Any insertions relative to the reference sequence are removed to
#' enforce a strict pairwise alignment.
#'
#' @param cigar: character, string in CIGAR format
#' @param seq: character, read sequence
#' @param qual: character, quality string
#' @param pos: integer, first position of the read
#' @return character vector, c(sequence, quality)
apply.cigar <- function(cigar, seq, qual, pos) {
    # prepare outputs
    new.seq <- paste0(rep('-', pos), collapse = '')
    new.qual <- paste0(rep('!', pos), collapse = '')
    insertions <- list()
    
    # validate CIGAR
    is.valid <- grepl("^((\\d+)([MIDNSHPX=]))*$", cigar)
    if (!is.valid) {
        stop("Error: Invalid CIGAR string: \"", cigar, "\"")
    }
    
    pat <- "\\d+[MIDNSHPX=]"
    tokens <- re.findall(pat, cigar)
    left <- 1
    
    for (token in tokens) {
        len <- as.integer(gsub("^(\\d+).+$", "\\1", token))
        operand <- gsub("^\\d+([MIDNSHPX=])$", "\\1", token)
        
        if (operand == 'M') {
            # append matching substring
            new.seq <- paste0(new.seq,
                substr(seq, left, left + len - 1), collapse = '')
            new.qual <- paste0(new.qual,
                substr(qual, left, left + len - 1), collapse = '')
            left <- left + len
        }
        else if (operand == 'D') {
            # deletion relative to reference
            new.seq <- paste0(new.seq,
                paste0(rep('-', len), collapse = ""), collapse = "")
            new.qual <- paste0(new.qual,
                paste0(rep('!', len), collapse = ""), collapse = "")
        }
        else if (operand == 'I') {
            # insertion relative to reference
            insertions[[length(new.seq) + 1]] <- c(
                substr(seq, left, left + len - 1),
                substr(qual, left, left + len - 1)
            )
            left <- left + len
        }
        else if (operand == 'S') {
            # soft clip, omit
            left <- left + len
        }
        else if (operand == 'H') {
            # hard clip, do not advance along query
            next
        }
        else {
            stop("Error in apply.cigar: unsupported token ", token)
        }
    }
    
    # append quality string to sequence as an attribute
    attr(new.seq, 'qual') <- new.qual
    attr(new.seq, 'insertions') <- insertions
    return(new.seq)
}


# functions from https://datadebrief.blogspot.com/2011/03/ascii-code-table-in-r.html
ord <- function(x) { strtoi(charToRaw(x),16L) }
chr <- function(n) { rawToChar(as.raw(n)) }


atomize <- function(seq, qcutoff) {
    vec <- strsplit(seq, "")[[1]]
    qual <- strsplit(attr(seq, 'qual'), "")[[1]]
    to.censor <- (vec != '-') & sapply(qual, function(q) ord(q) - 33 < qcutoff)
    vec[to.censor] <- 'N'
    return(vec)
}

mixtures <- list(
    'AC' = 'M', 'AG' = 'R', 'AT' = 'W', 'CG' = 'S', 'CT' = 'Y', 'GT' = 'K',
    'AN' = 'A', 'CN' = 'C', 'GN' = 'G', 'NT' = 'T',
    '-A' = 'A', '-C' = 'C', '-G' = 'G', '-T' = 'T', '-N' = 'N'
)
resolve.mixture <- list(
    'M' = c('A', 'C'), 'R' = c('A', 'G'), 'W' = c('A', 'T'),
    'S' = c('C', 'G'), 'Y' = c('C', 'T'), 'K' = c('G', 'T')
)

#' Combine paired-end reads into a single sequence.  Manage discordant
#' base calls on the basis of quality scores.
#'
#' @param seq1: character, the first read of a pair, already processed with
#'              apply.cigar
#' @param seq2: character, the second read of a pair, already processed with
#'              apply.cigar
#' @param qcutoff: integer, bases with a quality score below this cutoff
#'                 will be replaced with "N"
#' @param min.q.delta: integer, if two reads disagree on a base call, the
#'                     base with higher quality will be assigned if it
#'                     exceeds this difference.
#' @return character, merged sequence
merge.pairs <- function(seq1, seq2, qcutoff = 10, min.q.delta = 5) {
    if (is.null(attr(seq1, "qual")) || is.null(attr(seq1, "qual"))) {
        stop("Error: merge.pairs() requires seq1 and seq2 to be processed by apply.cigar()")
    }
    
    v1 <- atomize(seq1, qcutoff)
    v2 <- atomize(seq2, qcutoff)
    ldiff <- length(v1) - length(v2)
    if (ldiff > 0) {
        v2 <- c(v2, rep("N", ldiff))
    } else if (ldiff < 0) {
        v1 <- c(v1, rep("N", -ldiff))
    }
    pair <- cbind(v1, v2)
    mseq <- apply(pair, 1, function(x) {
        if (x[1] == x[2]) { x[1] }
        else {mixtures[[paste0(sort(x), collapse = "")]]}
    })
    return(paste0(mseq, collapse = ""))
}


#' calculate length of terminal gap prefix/suffix
#' @param seq:  character, sequence to examine
#' @param prefix:  if FALSE, check gap suffix
len.terminal.gap <- function(seq, prefix=TRUE) {
    hits <- ifelse(prefix,
        re.findall('^-+', seq),
        re.findall('-+$', seq))
    ifelse(length(hits) > 0, nchar(hits[1]), 0)
}

is.first.read <- function(flag) {
    # actually flag>64 would do the job just as well..
    bitwAnd(flag, 0x40) != 0
}


parse.sam.line <- function(line) {
    tokens <- strsplit(line, '\t')[[1]]
    list(qname = tokens[1], flag = tokens[2], rname = tokens[3],
        pos = as.integer(tokens[4]), mapq = as.integer(tokens[5]),
        cigar = tokens[6], seq = tokens[10], qual = tokens[11])
}

parse.sam <- function(infile, verbose = TRUE){
    t0 <- Sys.time()
    
    timings <- c("First Read" = NA, "Cigar" = NA, "Paired" = NA, "PosRange" = NA,
        "Phred" = NA, "Update" = NA)
    # All lines as Tibble
    if (verbose) {
        cat("First Read...\n")
    }
    t1 <- Sys.time()
    con <- file(infile, open = 'r')
    s <- readLines(con)
    s <- s[!grepl("^@", s)]
    s <- lapply(s, function(x){
        as.data.frame(parse.sam.line(x))
    })
    s <- as.data.frame(do.call(rbind, s))
    timings["First Read"] <- difftime(Sys.time(), t1, units = "mins")
    
    
    # Run all cigar values
    if (verbose) {
        cat("Have a cigar, you're gonna go far...\n")
    }
    t1 <- Sys.time()
    mseqs <- lapply(which(s$cigar != '*'), function(i) {
        x <- s[i,]
        apply.cigar(cigar = x$cigar, seq = x$seq,
            qual = x$qual, pos = x$pos)
    })
    timings["Cigar"] <- difftime(Sys.time(), t1, units = "mins")
    
    
    
    # Check for paired neighbours
    ###FOR TEST FILE, THIS ACTUALLY == 0?
    ###THIS MERGED/PAIRED STUFF IS UNTESTED
    if (verbose) {
        cat("Checking paired...\n")
    }
    t1 <- Sys.time()
    iPaired <- sapply(1:length(s$qname), function(i){
        sum(s$qname == s$qname[i]) > 1
    })
    timings["Paired"] <- difftime(Sys.time(), t1, units = "mins")
    
    
    
    
    # Calculates the longest an mseq value could be
    # Used to create matrix and prevent dynamic growth
    maxLen <- max(sapply(mseqs, function(mseq){nchar(mseq)}))
    m <- matrix(0, nrow = maxLen, ncol = 4)
    colnames(m) <- c('A', 'C', 'G', 'T')
    
    
    
    # For Targetting
    if (verbose) {
        cat("Finding position ranges...\n")
    }
    t1 <- Sys.time()
    posRanges <- lapply(1:length(mseqs), function(i) {
        # Range of positions in df that are effected by mseq values
        return(len.terminal.gap(mseqs[[i]]):nchar(mseqs[[i]]))
    })
    timings["PosRange"] <- difftime(Sys.time(), t1, units = "mins")
    
    
    
    # For increasing
    if (verbose) {
        cat("Calculating/Translating Phred scores...\n")
    }
    t1 <- Sys.time()
    alphabet <- c('A', 'C', 'G', 'T')
    posVals <- lapply(1:length(mseqs), function(i){
        
        mseq <- mseqs[[i]]
        posRange <- posRanges[[i]]
        
        # Calculates a matrix.
        # By adding the values of this matrix to df, we can update it
        res <- sapply(posRange, function(pos){
            nt <- substr(mseq, pos, pos)
            qc <- substr(attr(mseq, 'qual'), pos, pos)
            
            if (nt == '-') { #If a gap, then nothing is updated
                return(rep(0,4))
            } else if (nt == 'N') { # If an 'N', then everything up by 0.25
                return(rep(0.25,4))
            } else {# If some base than that base up by 1-p and other bases up by p/3
                p <- 10^-((ord(qc) - 30)/10)
                # if paired, divide probabilities by 2
                # a pair counts as a single obs, each read counts half.
                # ipaired[i] == TRUE means ipaired[i] + 1 == 2.
                temp <- rep((p/3),4) / (iPaired[i] + 1)
                temp[which(alphabet %in% nt)] <- (1 - p) / (iPaired[i] + 1)
                return(temp)
            }
        })
        
        return(res)
    })
    timings["Phred"] <- difftime(Sys.time(), t1, units = "mins")
    
    
    
    # Increase targeted areas
    if (verbose) {
        cat("Updating matrix...\n")
    }
    t1 <- Sys.time()
    for (i in 1:length(mseqs)) {
        # Add the Transposed values to the
        m[posRanges[[i]],] <- m[posRanges[[i]],] + t(posVals[[i]])
    }
    timings["Update"] <- difftime(Sys.time(), t1, units = "mins")
    
    #print(t0)
    totaltimings <- difftime(Sys.time(), t0, units = "mins")
    print(paste0("Total time: ", as.numeric(totaltimings)))
    close(con)
    
    if(verbose){
        cat("\nRelative time (percent):\n")
        print(round(timings/as.numeric(totaltimings), 4)*100)
    }
    
    m
}

parse.sam_deprecated <- function(infile, paired=FALSE, chunk.size=1000,
    save.partial=TRUE, est.length=NULL, time.step=20) {
    df <- matrix(0, nrow=1000, ncol=6)
    colnames(df) <- c('A', 'C', 'G', 'T', 'N', 'gap')
    
    con <- file(infile, open='r')
    read.next <- TRUE
    i <- 0
    t0 <- Sys.time()
    while (TRUE) {
        # progress monitoring and timing
        i <- i + 1
        if (i %% time.step == 0) {
            dt <- difftime(Sys.time(), t0, units = "secs")
            cat(paste0(i, " ", round(dt/i, 5), " s/line\n"))
            
            if(!is.null(est.length)) {
                dt_hr <- dt
                units(dt_hr) <- "hours"
                est.remain <- (est.length - i)*dt_hr/i
                cat("Estimated ", est.remain, " hours remaining.\n")
            }
            if(save.partial) saveRDS(df, file = "sam_partial.RDS")
        }
        
        if (read.next) {
            line <- readLines(con, n=1, warn=FALSE)
            if (length(line) == 0) {
                # end of file
                break
            }
        }
        else {
            read.next <- TRUE
        }
        
        if (grepl("^@", line)) {
            # skip header line
            next
        }
        row1 <- parse.sam.line(line)
        
        # look ahead for second read of pair
        if (paired) {
            line <- readLines(con, n=1, warn=FALSE)
            if (length(line) == 0) {
                # end of file - process row1 and exit
                if (row1$cigar == '*') {
                    # unmapped
                    next
                }
                mseq <- apply.cigar(cigar=row1$cigar, seq=row1$seq,
                    qual=row1$qual, pos=row1$pos)
            }
            else {
                row2 <- parse.sam.line(line)
                if (row1$cigar == '*' || row2$cigar == '*') {
                    # unmapped
                    next
                }
                
                seq1 <- apply.cigar(cigar=row1$cigar, seq=row1$seq,
                    qual=row1$qual, pos=row1$pos)
                
                if (row2$qname == row1$qname) {
                    seq2 <- apply.cigar(cigar=row2$cigar, seq=row2$seq,
                        qual=row2$qual, pos=row2$pos)
                    mseq <- merge.pairs(seq1, seq2)
                }
                else {
                    read.next <- FALSE  # use current line as row1
                    if (row1$cigar == '*') {
                        # unmapped
                        next
                    }
                    # second read is not mate of first - process first as unpaired
                    mseq <- seq1
                }
            }
        }
        else {
            # unpaired SAM
            if (row1$cigar == '*') {
                # unmapped
                next
            }
            mseq <- apply.cigar(cigar=row1$cigar, seq=row1$seq,
                qual=row1$qual, pos=row1$pos)
        }
        
        # update data frame
        start <- len.terminal.gap(mseq)
        aligned <- apply.cigar(cigar=row1$cigar, seq=row1$seq,
            qual=row1$qual, pos=row1$pos)
        for (pos in seq(start, nchar(mseq))) {
            while (pos > nrow(df)) {
                chunk <- matrix(0, nrow=chunk.size, ncol=6)
                colnames(chunk) <- colnames(df)
                df <- rbind(df, chunk)
            }
            nt <- substr(mseq, pos, pos)
            qc <- substr(attr(aligned, 'qual'), pos, pos)
            if (nt == '-') {
                nt <- 'gap'
            } else if (nt == 'N') {
                # nanopore sequences have no N's
                for (nt2 in c('A', 'C', 'G', 'T')) {
                    df[pos, nt2] <- df[pos, nt2]+0.25
                }
            }  else {
                q <- ord(qc)-30
                p <- 10^-(q/10)  # Phred conversion
                df[pos, nt] <- df[pos, nt]+(1-p)
                for (nt2 in setdiff(c('A','C','G','T'), nt)) {
                    df[pos, nt2] <- df[pos, nt2]+(p/3)
                }
            }
        }
    }
    
    close(con)
    max.row <- which.max(cumsum(apply(df, 1, sum)))
    return(df[1:max.row, ])
}

parse.sam.dt <- function(inFile, nc = 1, cs = 5000){
    
    print("Reading File")
    #Reading data now uses this function to read by chunks
    con <- file(inFile, open="rb")
    
    f <- function(x, pos) {
        s <- x[!grepl("^@", x)]
        s <- bind_rows(mclapply(s, function(z){
            parse.sam.line(z)
        }, mc.cores=nc))
        return(as.data.table(s))
    }
    
    s <- read_lines_chunked(con, callback=DataFrameCallback$new(f), chunk_size = 100)
    
    #Data file has now been read
    print("File Read")
    close(con)
    
    #Run all cigar values 
    mseqs <- mclapply(which(s$cigar!='*'), function(i) {
        x <- s[i,]
        apply.cigar(cigar=x$cigar, seq=x$seq, 
            qual=x$qual, pos=x$pos)
    }, mc.cores=nc)
    print("Cigar Strings Processed")
    
    print("Preparing Position Data")
    #Check for paired neighbours
    ###FOR TEST FILE, THIS ACTUALLY == 0?
    ###THIS MERGED/PAIRED STUFF IS UNTESTED
    iPaired <- sapply(1:(length(s$qname)-1), function(i){
        s$qname[i] == s$qname[i+1]
    })
    
    #Merge those paired neighbours mseq values
    mseqsPaired <- sapply(which(iPaired), function(i){
        merge.pairs(mseqs[[i]], mseqs[[i+1]])
    })
    
    #Calculates the longest an mseq value could be
    #Used to create matrix and prevent dynamic growth
    maxLen <- max(sapply(mseqs, function(mseq){nchar(mseq)}))
    m <- matrix(0, nrow=maxLen, ncol=4)
    colnames(m) <- c('A', 'C', 'G', 'T')
    
    #For Targetting
    posRanges <- mclapply(1:length(mseqs), function(i) {
        
        #aligned normally == mseq. Only changes to merged value when paired
        ###UNTESTED
        if(i%in%iPaired) {
            mseq <- mseqsPaired[[which(iPaired==i)]]
        } else {
            mseq <- mseqs[[i]]
        }
        
        #Range of positions in df that are effected by mseq values
        return(len.terminal.gap(mseq):nchar(mseq))
    }, mc.cores=nc)
    
    #For increasing
    posVals <- mclapply(1:length(mseqs), function(i){
        
        #aligned normally == mseq. Only changes to merged value when paired
        ###UNTESTED
        if(i%in%iPaired) {
            mseq <- mseqsPaired[[which(iPaired==i)]]
        } else {
            mseq <- mseqs[[i]]
        }
        
        aligned <- mseqs[[i]]
        posRange <- posRanges[[i]]
        
        #Calculates a matrix. By adding the values of this matrix to df, we can update it 
        res <- sapply(posRange, function(pos){
            nt <- substr(mseq, pos, pos)
            qc <- substr(attr(aligned, 'qual'), pos, pos)
            
            if(nt == '-') { #If a gap, then nothing is updated
                return(rep(0,4))
                
            } else if(nt == 'N') { #If an 'N', then everything up by 0.25
                return(rep(0.25,4))
                
            } else { #If some base than that base up by 1-p and other bases up by p/3
                p <- 10^-((ord(qc)-30)/10)
                temp <- rep((p/3),4)
                temp[which(c('A', 'C', 'G', 'T') %in% nt)] <- 1-p
                return(temp)
            }
        })
        
        return(res)
    }, mc.cores=nc)
    
    print("Final Step: Applying position Values")
    
    #Increase targetted areas
    for(i in 1:length(mseqs)) {
        #Add the Transposed values to the 
        m[posRanges[[i]],] <- m[posRanges[[i]],] + t(posVals[[i]])
    }
    
    m
}


#' Parallel implementation specifically for unpaired
#' read files.
#' TODO: extend for paired read data
#' @param infile:  absolute or relative path to SAM file
#' @param n.cores:  number of cores to run in parallel
#' @return data frame of nucleotide counts
parse.sam.mp_deprecated <- function(infile, n.cores=2) {
    # parse header information to allocate data frame
    # FIXME: assumes only one reference!
    tot.len <- 0
    skip <- 0
    con <- file(infile, open='r')
    while (length(line <- readLines(con, n=1, warn=FALSE)) > 0) {
        if (grepl("^@", line)) {
            skip <- skip + 1
            if (grepl("@.+LN:[0-9]+", line)) {
                tot.len <- as.integer(gsub(".+LN:([0-9]+)$", "\\1", line))
            }
        } else {
            # assume first line not prefixed with "@" exits comment block
            break
        }
    }
    close(con)
    
    # load the entire file contents
    cat("Loading SAM file...\n")
    sam <- read.table(infile, sep='\t', skip=skip, quote="", fill=TRUE,
        stringsAsFactors = FALSE, comment.char="")
    sam <- sam[,1:11]  # strip optional fields
    names(sam) <- c('qname', 'flag', 'rname', 'pos', 'mapq', 'cigar',
        'rnext', 'pnext', 'tlen', 'seq', 'qual')
    
    cat("Launching processes...\n")
    res <- mclapply(1:n.cores, function(i) {
        df <- data.frame(matrix(0, nrow=tot.len, ncol=6))
        names(df) <- c('A', 'C', 'G', 'T', 'N', 'gap')
        for (j in seq(i, nrow(sam), by=n.cores)) {
            row <- sam[j, ]
            if (row$cigar == '*' || any(is.na(row))) {
                next
            }
            aligned <- apply.cigar(cigar=row$cigar, seq=row$seq,
                qual=row$qual, pos=row$pos)
            start <- len.terminal.gap(aligned)
            for (pos in seq(start, nchar(aligned))) {
                nt <- substr(aligned, pos, pos)
                qc <- substr(attr(aligned, 'qual'), pos, pos)
                if (nt == '-') {
                    # '-' is illegal column name
                    df[pos, 'gap'] <- df[pos, 'gap'] + 1
                }
                else if (nt == 'N') {
                    # nanopore sequences have no N's
                    for (nt2 in c('A', 'C', 'G', 'T')) {
                        df[pos, nt2] <- df[pos, nt2]+0.25
                    }
                }
                else {
                    q <- ord(qc)-30
                    p <- 10^-(q/10)  # Phred conversion
                    df[pos, nt] <- df[pos, nt]+(1-p)
                    for (nt2 in setdiff(c('A','C','G','T'), nt)) {
                        df[pos, nt2] <- df[pos, nt2]+(p/3)
                    }
                }
            }
        }
        df  # return data frame
    }, mc.cores=n.cores)
    Reduce('+', res)
}
