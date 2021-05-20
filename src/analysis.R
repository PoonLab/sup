library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)
library(here)
library(ggridges)

source("dist-fcts.R")
source("utils.R")

message("Analysis.R: Starting analysis...")

system("figgy=`whereis figlet`; figgylen=${#figgy}; if [ $figgylen > 10 ]; then figlet Analysis.R; fi")

# ---- Data ----

rdatas <- system("ls *treedist*RData", intern = TRUE)

# ---- Distances ----

btwn_tmp  <- list()
certain_tmp <- list()
for (i in seq_along(rdatas)) {
    load(rdatas[[i]])
    # between each inferred trees
    btwn_tmp[[i]]  <- data.frame(d.rf   = dist.list$d.rf,
                            d.kf   = dist.list$d.kf,
                            d.sh   = dist.list$d.sh,
                            #d.kern = dist.list$d.kern,
                            prmset = as.numeric(dist.list$prmset))
    # benchmark (difference from "true" tree)
    certain_tmp[[i]] <- data.frame(d.rf   = dist.list$d.rf.star,
                            d.kf   = dist.list$d.kf.star,
                            d.sh   = dist.list$d.sh.star,
                            #d.kern = dist.list$d.kern.star,
                            prmset = as.numeric(dist.list$prmset))
    prmsimlabel <- dist.list[["prmsimlab"]]
}
between_df  <- do.call("rbind", btwn_tmp) # between inferred
certain_df <- do.call("rbind", certain_tmp) # inferred versus certain


# Retrieve the sequence entropy for each prm set:
get_entropy_prmset <- function(fname = "entropy-prmset.out") {
    x <- read.csv(fname, header = F)
    names(x) <- c("prmset", "s1", "s2", "entropy")
    return(x)
}

digest_distances <- function(df) {
    p_lo <- 0.025
    p_hi <- 0.975

    df_ent <- get_entropy_prmset()

    tmp <- df %>%
        pivot_longer(cols = -prmset,
                     names_to = "distance.type",
                     values_to = "value") %>%
        group_by(prmset, distance.type) %>%
        summarise(m = mean(value, na.rm = TRUE),
                  md = median(value, na.rm = TRUE),
                  s = sd(value, na.rm = TRUE),
                  q.lo = quantile(value, probs = p_lo, na.rm = TRUE),
                  q.hi = quantile(value, probs = p_hi, na.rm = TRUE),
                  minv = min(value, na.rm = TRUE),
                  maxv = max(value, na.rm = TRUE),
                  n = sum(!is.na(value)),
                  .groups = "drop") %>%
        mutate(cv = s / m)

    beta_var <- function(a, b) (a * b) / ((a + b)^2 * (a + b + 1))
    res <- left_join(tmp, df_ent, by = "prmset") %>%
        mutate(beta_mean = s1 / (s1 + s2),
            beta_var = beta_var(s1, s2))

    return(res)
}



between  <- digest_distances(between_df)
certain <- digest_distances(certain_df)

saveRDS(between,
    file = here("data", "output", "between-inferred-distances.RDS"))
saveRDS(certain,
    file = here("data", "output", "inferred-to-certain-distances.RDS"))

# Pre-plot wrangling

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

options(scipen = 6)
meanset <- select(certain, prmset, s1, s2,
        beta_mean, beta_var) %>%
    distinct() %>%
    filter(prmset != 11) %>%
    arrange(prmset) %>%
    mutate(beta_mean = sapply(beta_mean, pad, -4)) %>%
    mutate(beta_var = sapply(beta_var, pad, -4)) %>%
    mutate(meansd = factor(paste(beta_mean, beta_var,
            sep = ", "),
        ordered = TRUE))
    
xseq <- seq(0, 0.075, 0.0001)
beta_dists <- lapply(seq_len(nrow(meanset)), function(i) {
    data.frame(x = xseq,
        y = dbeta(xseq,
            shape1 = meanset$s1[i],
            shape2 = meanset$s2[i]),
        label = meanset$meansd[i]
    )
}) %>% bind_rows()

ggplot(beta_dists) +
    aes(x = x, y = y, colour = label) +
    geom_line() +
    geom_line(size = 1) +
    scale_colour_viridis_d() +
    coord_cartesian(ylim = c(0, 40)) +
    theme_dark() +
    theme(legend.position = c(0.9, 0.7)) +
    labs(x = "x", y = "Probability Density",
        colour = "Parameters")

certain_long <- pivot_longer(certain_df, -prmset,
    names_to = "distance.type", values_to = "value") %>%
    filter(prmset != 11) %>%
    dplyr::left_join(meanset, by = "prmset") %>%
    mutate(distance = case_when(
        distance.type == "d.kf" ~ "Kuhner-Felsenstein",
        distance.type == "d.rf" ~ "Robinson-Foulds",
        distance.type == "d.sh" ~ "Shared Ancestry",
        TRUE ~ as.character(distance.type)
    ))
rm(meanset)

ggplot(certain_long) +
    aes(x = value, y = meansd, fill = meansd) +
    geom_density_ridges(bandwidth = 0.03) +
    facet_wrap(~ distance, scales = "free_x") +
    labs(y = "Mean and SD of Beta Dist",
        x = 'Distance from "certain" tree',
        title = 'KDE for distances between "certain" and simulated trees',
        subtitle = paste("Bandwidth = 0.03",
            "1000 trees for each parameter combo",
            "ordered by mean error rate (but note that the sd changes)",
            sep = "; ")) +
    theme(legend.position = "none")

# ---- TN93 distances ----

# Warning, these are distances between sequences
# in one tree, not distances between trees!

# The raw TN93 distances:
df_tn93 <- dist.tn93()

# Number of clusters based on TN93 distance:
thresh <- c(0.02, 0.30)

df_tn93_clustr_1 <- df_tn93 %>%
    clstr_num(dist.thresh.mean = thresh[1])

df_tn93_clustr_2 <- df_tn93 %>%
    clstr_num(dist.thresh.mean = thresh[2])
