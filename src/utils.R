

xlogx <- function(x) {
    x <- x[x>0]
    return(x * log2(x))
}


entropy <- function(p){
    return(-sum(xlogx(p)))
}

#' Converts a dataframe of nucleotide and associated 
#' proba of observation to a probabilistic sequence
#' filling the probas for the other nucleotides
#' using a uniform distribution.
#' @param test.mode Boolean. Use for testing only!
#' 
pmaxseq_to_probseq <- function(df, test.mode = FALSE) {
    # Create empty structure:
    N = nrow(df)
    m <- data.frame(matrix(ncol=4, nrow = N))
    names(m) <- c('A','C','G','T')
    
    # Fill the input probabilities:
    if(!test.mode){
        for(i in 1:N){ # i=1
            m[i,df$base[i]] <- df$p[i]
        }
    }
    if(test.mode){
        m[,1] <- df$p
    }
    
    # Fill the 3 other probas using uniform distribution:
    if(!test.mode){
        for(i in 1:N){ # i=1
            y  <- runif(n=3)
            yy <- y/sum(y) * (1 - df$p[i])
            m[i, is.na(m[i,])] <- yy
        }
    }
    if(test.mode){
        mu <- matrix(data = runif(n = N*3), ncol=3)
        smu <- apply(mu, 1, sum)
        mu2 <- mu / smu * (1-df$p)
        m[,2:4] <- mu2
    }
    
    return(m)
}


get_prm <- function(prm, prm.name) {
    return(prm$value[prm$name == prm.name])
}

#' Define a probabilistic sequence, randomly assigning probabilities.
#' Four rows for each nucleotide: A,C,G,T
#' Column = nucleotide position
#' @param seq.length Integer. Length of the sequence
#' @param p Numerical. Fixed probability for one nucleotide. 
#' For example, if p=0.9 a randomly chosen nucleotide among A,C,G,T will
#' be assigned probability = 0.9. The three other will be assign
#' random numbers (such that sum=1).
#' @return Matrix (nrow=4, ncol=seq.length). Probabilistic sequence. 
build_random_seq <- function(seq.length, p = NULL) {
    
    M <- matrix(ncol = seq.length, 
                nrow = 4)
    
    if(is.null(p)){
        # Fill with probabilities:
        for(j in 1:seq.length){
            x <- runif(4)
            M[,j] <- x / sum(x) # proba must sum to 1.0
        }    
    }
    
    if(!is.null(p)){
        for(j in 1:seq.length){
            y <- runif(3)
            y <- y / sum(y) * (1-p) # proba must sum to 1.0
            tmp <- c(p,y)
            x <- tmp[sample(1:length(tmp))]
            M[,j] <- x
        }  
    }
    return(M)
}

#' Randomly assign probabilities from an existing sequence.
#' Four rows for each nucleotide: A,C,G,T
#' Column = nucleotide position
#' @param seq String. Sequence nucleotides.
#' @param p Numerical vector. Fixed probability for each position for the nucleotide given as input. 
#' For example, if p=0.9, the nucleotide given at any position will
#' be assigned probability = 0.9. The three other will be assign
#' random numbers (such that sum=1).
#' @return Matrix (nrow=4, ncol=seq.length). Sequence probabilities. 
build_from_seq <- function(seq, p) {
    # seq = 'GCCATA'
    
    # if character vector format, convert to string:
    if(length(seq) > 1) 
        seq = paste(seq,collapse = '')
    
    # If only one value is given, 
    # then same proba for all positions:
    if(length(p)==1) p <- rep(p, length(seq))
    
    seq.length <- nchar(seq)
    M <- matrix(ncol = seq.length, 
                nrow = 4)
    
    v <- strsplit(x = seq, split = '')[[1]]
    idx <- sapply(v, tr_idx)
    
    for(j in 1:seq.length){   # j=1
        y <- runif(3)
        y <- y / sum(y) * (1-p[j]) # proba must sum to 1.0
        x <- numeric(4)
        x[idx[j]] <- p[j]
        x[-idx[j]] <- y
        M[,j] <- x
    }  
    return(M)
}


build_from_zanini <- function(pos.start, pos.end, df) {
    
}



#' Translate index to nucleotide character.
tr_nucl <- function(i){
    res <- NA
    if(i==1) res ='A'
    if(i==2) res ='C'
    if(i==3) res ='G'
    if(i==4) res ='T'
    return(res)
}

#' Translate nucleotide character to index.
tr_idx <- function(nucl){
    res <- NA
    if(nucl=='A') res = 1
    if(nucl=='C') res = 2
    if(nucl=='G') res = 3
    if(nucl=='T') res = 4
    return(res)
}

#' Draw the jth nucleotide from 
#' a sequence probability distribution.
draw_nucleo <- function(j, M) {
    x <- rmultinom(n = 1, size = 1, prob = M[,j])
    res <- tr_nucl(which(x==1))
    return(res)
}

#' Draw a sequence from a sequence probability distribution.
draw_seq <- function(M, vect.format = TRUE){
    s <- sapply(1:ncol(M), FUN = draw_nucleo, M=M)
    if(!vect.format) s <- paste(s, collapse = '')
    return(s)
}



#' Draw multiple sequence probabilities tables.
#' @param n.seq Integer. Number of sequences to draw.
#' @param seq.length Integer. Length for each sequences (they all have the same length).
#' @return A list of sequence probabilities.
build_multiple_seqproba <- function(n.seq, 
                                    seq.length, 
                                    p = NULL) {
    # Random sequence proba distributions:
    M.list <- list()
    for(i in 1:n.seq){
        M.list[[i]] <- build_random_seq(seq.length, p)
    }
    return(M.list)
}


#' Draw multiple sequences. Each sequence was drawn from
#' a given sequence probability distribution.
#' @param M.list List of sequence probabilities. 
#' @param filename String. Name of the FASTA file where all the sequences drawn are saved.
#' @return A list of string vectors, each list element representing the sequence
draw_multiple_seq <- function(M.list,
                              filename) {
    
    # Draw the sequences:
    s.list <- lapply(M.list, draw_seq, vect.format = T)
    
    # Save in FASTA file:
    nm <- paste('seq',1:length(s.list), sep = '_')
    write.fasta(s.list, 
                names = nm, 
                nbchar = 50,
                file.out = filename)
    
    return(s.list)
}
