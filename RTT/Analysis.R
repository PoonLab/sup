# Finding the substitution rates and sd
library(dplyr)
library(lubridate)
library(here)
library(ggplot2)

get_slope_sd <- function(dir) {
	temp1 <- readLines(here(dir, "molecular_clock.txt"))
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

tree_dirs <- tree_dirs[tree_dirs$date == max(tree_dirs$date),]


tree_dirs <- lapply(tree_dirs$dir, get_slope_sd) %>% 
	bind_rows %>% 
	right_join(tree_dirs, by = "dir") %>% 
	select(-date)

raw <- get_slope_sd("raw_tree")
raw$sim = "raw"

rbind(raw, tree_dirs)

ggplot() +
	geom_hline(yintercept = raw$slope) + 
	geom_rect(
		mapping = aes(xmin = -10, xmax = 60, 
			ymin = raw$slope - raw$sd, ymax = raw$slope + raw$sd),
		fill = "red", alpha = 0.2) + 
	geom_errorbar(
		mapping = aes(x = as.numeric(sim), ymin = slope - sd, ymax = slope + sd),
		data = tree_dirs) +
	coord_cartesian(xlim = c(0,50))


# Next steps: re-do the original tree with just the ones that were actually sampled
