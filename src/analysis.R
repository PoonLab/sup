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

between <- readRDS(file = here("data", "output", 
    "between-inferred-distances.RDS"))
certain <- readRDS(file = here("data", "output", 
    "inferred-to-certain-distances.RDS"))

# Pre-plot wrangling
options(scipen = 10)
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

# TODO: Find the problem in the TN93 code
    # Currently, NO clusters are found.
    # Not sure why/how the code calculates clusters