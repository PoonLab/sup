# adapted from http://github.com/PoonLab/MiCall-Lite

#' apply.cigar
#' 
#' Use CIGAR (Compact Idiosyncratic Gapped Alignment Report) string
#' to applly soft clips, insertions and deletions to the read sequence.
#' Any insertions relative to the reference sequence are removed to 
#' enforce a strict pairwise alignment.
#' 
#' @param cigar: string in CIGAR format
#' @param seq: read sequence
#' @param qual: quality string
#' @param pos: first position of the read
apply.cigar <- function(cigar, seq, qual, pos=0) {
  pos <- as.integer(pos)
  
  # prepare outputs
  new.seq <- rep('-', pos)
  new.qual <- rep('!', pos)
  insertions <- list()
  is.valid <- grepl("^((\d+)([MIDNSHPX=]))*$", cigar)
  if (!is.valid) {
    stop("Error: Invalid CIGAR string", cigar)
  }
  
}

parse.sam.line <- function(line) {
  
}