library(ggplot2)
library(dplyr)
options(scipen=6)

pad <- function(x, pad = -3){
    x <- as.character(x)
    if (pad < 0) {
        if (grepl("\\.", x)) {
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
    } else {
        n <- nchar(x)
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
bt <- rbind(
    mom(0.0001, 0.00001),
    mom(0.0002, 0.00001),
    mom(0.0005, 0.00001),
    mom(0.001, 0.0001),
    mom(0.005, 0.0001),
    mom(0.01, 0.001),
    mom(0.05, 0.01)
)

xseq <- seq(0, 0.075, 0.0001)

bts <- lapply(seq_len(nrow(bt)), function(x) {
    thislabel <- paste0(
        "x̄=", pad(bt[x, 1], -4),
        ", s=", pad(bt[x, 2], -5),
        " | α=", pad(round(bt[x, 3], 3), -3),
        ", β=", pad(round(bt[x, 4], 3), -3),
        collapse = ""
    )
    data.frame(x = xseq,
        y = dbeta(xseq,
            shape1 = bt[x, 3],
            shape2 = bt[x, 4]),
        prm = thislabel)
}) %>% bind_rows()

ggplot(bts) +
    aes(x = x, y = y, colour = prm) +
    geom_line() +
    scale_colour_brewer(type = "div", palette = 8) +
    coord_cartesian(ylim = c(0, 40)) +
    theme_dark() +
    labs(x = "x", y = "Probability Density",
        colour = "Parameters")

write.table(bt[, 3:4], file = "prm.btshp",
    row.names = FALSE, col.names = FALSE)
