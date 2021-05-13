library(ggplot2)
library(dplyr)
options(scipen=6)

pad <- function(x, pad = -3){
    x <- as.character(x)
    if (length(gregexpr("\\.", x)[[1]]) > 1) {
        stop("Invalid number.")
    }
    if (pad < 0) {
        if (grepl("\\.", x)) {
            # nchar of everything after the decimal
            n <- nchar(strsplit(x, "\\.")[[1]][2])
            if (n < abs(pad)) {
                x <- paste0(x, 
                    paste0(rep(0, abs(pad) - n),
                        collapse = ""),
                    collapse = "")
            }
        } else {
            x <- paste0(x, ".")
            x <- paste0(x,
                paste0(rep(0, abs(pad)), collapse = ""),
                collapse = "")
        }
    } else { # pad > 0
    # nchar of everything before the decimal
        n <- nchar(strsplit(x, "\\.")[[1]][1])
        x <- paste0(rep(0, max(0, pad - n)), x)
    }
    x
}

# Method of Moments for Beta distribution
mom <- function(xbar, s){
    tmp <- xbar * (1 - xbar) / s - 1
    alpha <- xbar * tmp
    beta <- (1 - xbar) * tmp
    return(c(xbar = xbar, s = s, alpha = alpha, beta = beta))
}

# Choose beta parameters based on mean and sd
prmset <- rbind(
    mom(0.0002, 0.0001),
    mom(0.0025, 0.001),
    mom(0.005, 0.001),
    mom(0.008, 0.001),
    mom(0.0125, 0.001),
    mom(0.05, 0.005)
)

xseq <- seq(0, 0.075, 0.0001)

bts <- lapply(seq_len(nrow(prmset)), function(x) {
    thislabel <- paste0(
        "x̄=", pad(prmset[x, 1], -4),
        ", s=", pad(prmset[x, 2], -4),
        " | α=", pad(round(prmset[x, 3], 3), -3),
        ", β=", pad(round(prmset[x, 4], 3), -3),
        collapse = ""
    )
    data.frame(x = xseq,
        y = dbeta(xseq,
            shape1 = prmset[x, 3],
            shape2 = prmset[x, 4]),
        prm = thislabel)
}) %>% bind_rows()

ggplot(bts) +
    aes(x = x, y = y, colour = prm) +
    geom_line(size = 1) +
    scale_colour_viridis_d(option = 1) +
    coord_cartesian(ylim = c(0, 40)) +
    theme_dark() +
    theme(legend.position = c(0.8, 0.8)) +
    labs(x = "x", y = "Probability Density",
        colour = "Parameters")

write.table(prmset[, 3:4], file = "prm-btshp.csv",
    row.names = FALSE, col.names = FALSE, sep = ", ")
