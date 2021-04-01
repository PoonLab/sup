library(dplyr)
library(tidyr)

sam1 <- readRDS("data/unc_covid/open_seSAMe-S-ERR4363387.RDS")


dim(sam1)

if (ncol(sam1) == 6) {
    sam1[, 5] <- sam1[, 5] + sam1[, 6]
    sam1 <- sam1[, 1:5]
}

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
dim(sam1) == dim(sam2)

to_sample <- apply(sam2, 1, function(x) {
    (max(x) < 1) & (!is.na(max(x)))
})
table(to_sample)

non0 <- sam2[to_sample, ] %>%
    `colnames<-`(1:5) %>%
    as.data.frame() %>% 
    tibble::rowid_to_column(var = "row_id") %>%
    tidyr::pivot_longer(cols = -row_id, names_to = "col_id") %>%
    mutate(order = order(value))

# Next steps:
    # 1. Fix everything
    # 2, Figure it out from there.

