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
  new.seq <- paste0(rep('-', pos), collapse='')
  new.qual <- paste0(rep('!', pos), collapse='')
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
      new.seq <- paste0(new.seq, substr(seq, left, left+len-1), collapse='')
      new.qual <- paste0(new.qual, substr(qual, left, left+len-1), collapse='')
      left <- left + len
    }
    else if (operand == 'D') {
      # deletion relative to reference
      new.seq <- paste0(new.seq, rep('-', len), collapse='')
      new.qual <- paste0(new.qual, rep('!', len), collapse='')
    }
    else if (operand == 'I') {
      # insertion relative to reference
      insertions[[length(new.seq)+1]] <- c(
        substr(seq, left, left+len-1),
        substr(qual, left, left+len-1)
        )
      left <- left + len
    }
    else if (operand == 'S') {
      # soft clip, omit
      left <- left + len
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
merge.pairs <- function(seq1, seq2, qcutoff=10, min.q.delta=5) {
  mseq <- ''
  
  if (is.null(attr(seq1, "qual")) || is.null(attr(seq1, "qual"))) {
    stop("Error: merge.pairs() requires seq1 and seq2 to be processed by apply.cigar()")
  }
  
  # force second read to be the longer of the two
  if (nchar(seq1) > nchar(seq2)) {
    temp <- seq1
    seq1 <- seq2
    seq2 <- seq1
  }
  
  # retrieve quality strings
  qual1 <- attr(seq1, "qual")
  qual2 <- attr(seq2, "qual")
  
  qcutoff.char <- chr(qcutoff+33)
  is.forward.started <- FALSE
  is.reverse.started <- FALSE
  
  for (i in 1:nchar(seq2)) {
    c2 <- substr(seq2, i, i)
    if (c2 != '-' && !is.reverse.started) {
      is.reverse.started <- TRUE
    }
    
    if (i <= nchar(seq1)) {
      c1 <- substr(seq1, i, i)
      
      if (!is.forward.started) {
        if (c1 == '-' && c2 == '-') { next }
        is.forward.started <- TRUE
        mseq <- substr(seq1, 1, i-1)
      }
      else {
        if (c1 == '-' && c2 == '-') {
          mseq <- paste0(mseq, '-')
          next
        }
      }
      
      q1 <- ord(substr(qual1, i, i))
      q2 <- ord(substr(qual2, i, i))
      
      if (c1 == c2) {
        # reads agree on base
        if (q1 > qcutoff || q2 > qcutoff) {
          # at least one has sufficient confidence
          mseq <- paste0(mseq, c1)
        } else {
          # neither base call is confident
          mseq <- paste0(mseq, "N")
        }
      }
      else {
        if (abs(q1-q2) >= min.q.delta) {
          if (q1 > max(q2, qcutoff)) {
            mseq <- paste0(mseq, c1)
          } else if (q2 > max(q1, qcutoff)) {
            mseq <- paste0(mseq, c2)
          } else {
            mseq <- paste0(mseq, 'N')
          }
        }
        else {
          mseq <- paste0(mseq, 'N')
        }
      }
    }
    else {
      # past end of read 1
      if (c2 == '-') {
        if (is.reverse.started) {
          mseq <- paste0(mseq, c2)
        } else {
          mseq <- paste0(mseq, 'N')
        }
      } 
      else if (q2 > qcutoff) {
        mseq <- paste0(mseq, c2)
      }
      else {
        mseq <- paste0(mseq, 'N')
      }
    }
  }
  
  # TODO: reconcile insertions
  return (mseq)
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
  list(qname=tokens[1], flag=tokens[2], rname=tokens[3],
       pos=as.integer(tokens[4]), mapq=as.integer(tokens[5]),
       cigar=tokens[6], seq=tokens[10], qual=tokens[11])
}


parse.sam <- function(infile, paired=TRUE, chunk.size=1000) {
  df <- data.frame(matrix(0, nrow=1000, ncol=6))
  names(df) <- c('A', 'C', 'G', 'T', 'N', 'gap')
  
  con <- file(infile, open='r')
  read.next <- TRUE
  i <- 0
  
  while (TRUE) {
    i <- i + 1
    if (i %% 100 == 0) {
      print(i)
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
    if (row1$cigar == '*') {
      
    }
    
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
    for (pos in seq(start, nchar(mseq))) {
      while (pos > nrow(df)) {
        chunk <- data.frame(matrix(0, nrow=chunk.size, ncol=6))
        names(chunk) <- names(df)
        df <- rbind(df, chunk)
      }
      nt <- substr(mseq, pos, pos)
      if (nt == '-') {
        nt <- 'gap'
      }
      df[pos, nt] <- df[pos, nt]+1
    }
  }
  
  close(con)
  max.row <- which.max(cumsum(apply(df, 1, sum)))
  return(df[1:max.row, ])
}








