# Collection of helpful functions


parse_accession <- function(filename) {
    # Get the accession number from filenames
    # Accessions are either SRR[0-9]{8} or ERR[0-9]{7}.
    # Making this a function is easier than my usual
    # workflow of looking up regex every time.
    serr <- regexpr("ERR|SRR", filename)
    num_length <- ifelse(grepl("ERR", filename), 7, 8)
    unlist(substr(filename, serr, serr + num_length + 2))
}

# Testing
if (FALSE) {
    filename <- c("ERR5082673", "SRR11433882", "SRR11433888.sam",
        "S-ERR5081077.sam", "open-sesame-S-SRR12345678.RDS")
    parse_accession(filename)
}




fix_unc <- function(unc_mat) {
    # The uncertainty matrices can have a myriad of errors
    # To use: after loading, fix_unc(), then
    # check if "character" is in the class of unc_mat.
    # If so, skip that accession.
    msg <- ""
    if (is.null(dim(unc_mat))) {
        msg <- "Empty File"
        return(msg)
    } else if (any(unc_mat[!is.na(unc_mat)] < 0 |
            unc_mat[!is.na(unc_mat)] > 10e8)) {
        msg <- paste0("Values too small or too large")
        return(msg)
    } else {
        if ("X" %in% colnames(unc_mat)) {
            unc_mat <- unc_mat[, -which(colnames(unc_mat) == "X")]
        }
        if (ncol(unc_mat) == 6) {
            unc_mat[, 5] <- unc_mat[, 5] + unc_mat[, 6]
            unc_mat <- unc_mat[, 1:5]
        }
        return(unc_mat)
    }
}
