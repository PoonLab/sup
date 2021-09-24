# Finding the substitution rates and sd
library(dplyr)
library(lubridate)
library(here)
setwd(here("RTT/sampled_trees"))
library(ggplot2)
library(ggrepel)

get_slope_sd <- function(dir) {
    temp1 <- readLines(here("RTT/sampled_trees", dir, "molecular_clock.txt"))
    temp1 <- strsplit(temp1[2], split = " ")[[1]]
    slope <- as.numeric(strsplit(temp1[2], split = "\t")[[1]][2])
    st_dev <- as.numeric(temp1[4])
    data.frame(slope = slope, sd = st_dev, dir = dir)
}

tree_dirs <- list.dirs()
tree_dirs <- tree_dirs[grep(x = tree_dirs, pattern = "treetime")]
tree_dirs <- gsub(pattern = "\\./", replacement = "", x = tree_dirs)
tree_dirs <- lapply(strsplit(tree_dirs, split = "_"), function(x) {
    datestring <- strsplit(x[1], split = "-")[[1]]
    date <- paste(datestring[1:3], collapse = "-")
    sim <- ifelse(length(datestring) == 4, datestring[4], "0000")
    data.frame(date = ymd(date), sim = sim, dir = paste(x, collapse = "_"))
}) %>% bind_rows()

tree_dirs <- tree_dirs[tree_dirs$date == max(tree_dirs$date), ]


tree_dirs <- lapply(tree_dirs$dir, get_slope_sd) %>%
    bind_rows %>%
    right_join(tree_dirs, by = "dir") %>%
    select(-date) %>% 
    mutate(rank = rank(slope))


raw <- get_slope_sd("raw_tree")
raw$sim <- "raw"

rbind(raw, select(tree_dirs, -rank))


lit_clock <- as.data.frame(rbind(
    c("Duchene et al. 2020 (n=122, 95% BCI)", 1.1e-3, 7.03e-4, 1.5e-3),
    c("Choudhary et al. 2021 (n=261, 95% HPD)", 6.77e-4, 5.91e-4, 7.66e-4),
    c("Song et al. 2021 (n=29, 95% BCI)", 9.25e-4, 6.75e-4, 1.28e-3),
    c("Nie et al. 2020 (n=112, 95% BCI)", 9.9e-4, 6.29e-4, 1.35e-3),
    c("Giedelberg et al. 2021 (n=77, 95% HPD)", 1.3e-3, 9.8e-4, 1.7e-3)
))
names(lit_clock) <- c("study", "clock", "lo95", "hi95")
lit_clock_x <- (max(as.numeric(tree_dirs$rank))+1):
                (max(as.numeric(tree_dirs$rank)) + 
                    nrow(lit_clock))

save(raw, tree_dirs, lit_clock, 
    file = here("RTT", "RTT_results.RData"))


png(file = here("RTT", "Results_Slope.png"), width = 600, height = 500)
ggplot() +
    theme_bw() +
    geom_hline(yintercept = 0, colour = "darkgrey") +
    geom_hline(yintercept = raw$slope, col = "red") +
    geom_rect(
        mapping = aes(xmin = -10, xmax = 100,
            ymin = raw$slope - 1.96*raw$sd, ymax = raw$slope + 1.96*raw$sd),
        fill = "red", alpha = 0.2) +
    geom_errorbar(
        mapping = aes(x = as.numeric(rank),
            ymin = slope - 1.96*sd,
            ymax = slope + 1.96*sd),
        data = tree_dirs) +
    geom_point(
        mapping = aes(x = as.numeric(rank), y = slope),
        data = tree_dirs) +
    geom_errorbar(
        mapping = aes(x = lit_clock_x,
            ymin = as.numeric(lo95),
            ymax = as.numeric(hi95)),
        data = lit_clock,
        colour = "darkorchid") +
    geom_point(
        data = lit_clock,
        mapping = aes(x = lit_clock_x, 
            y = as.numeric(clock)),
        colour = "darkorchid"
        ) +
    geom_text(data = lit_clock,
        mapping = aes(
            y = as.numeric(hi95), 
            x = lit_clock_x,
            label = study),
        colour = "black", 
        angle = 90,
        hjust = -0.05
        ) +
    scale_colour_viridis_d() +
    coord_cartesian(xlim = c(0, 52)) +
    labs(x = "Index (Ordered by Slope)", 
        y = "Slope +/- 1.96 SD",
        title = "Root-to-Tip Slopes from collections of re-sampled genomes",
        subtitle = "Compared to conseqs (red) and published estimates (purple)") +
    theme(legend.position = "none")
dev.off()




if (FALSE) {
    # What is the estimated date? -----------------------------

    dates <- read.csv(here("RTT/sampled_trees", tree_dirs$dir[1], "dates.tsv"),
        sep = "\t")

    get_node0_date <- function(dir) {
        dates <- read.csv(here("RTT/sampled_trees", dir, "dates.tsv"),
            sep = "\t")
        data.frame(ancestor_date = min(as.numeric(dates$numeric.date)), 
            ancestor_ymd = as.numeric(dates$date[which.min(dates$numeric.date)]),
            dir = dir, stringsAsFactors = FALSE)
    }

    tree_dirs <- lapply(tree_dirs$dir, get_node0_date) %>%
        bind_rows() %>%
        right_join(tree_dirs, by = "dir")

    tree_dirs

    raw <- get_node0_date("raw_tree")

    hist(tree_dirs$ancestor_date)
}