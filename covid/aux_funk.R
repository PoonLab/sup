# Collection of helpful functions



parse_accession <- function(filename) {
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
