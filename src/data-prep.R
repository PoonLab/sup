library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)
library(here)
library(ggridges)

source("dist-fcts.R")
source("utils.R")

message("data-prep.R: Starting analysis...")

system("figgy=`whereis figlet`; figgylen=${#figgy}; if [ $figgylen > 10 ]; then figlet data-prep.R; fi")

# ---- Data ----

rdatas <- system("ls *treedist*RData", intern = TRUE)
do_tn93 <- FALSE

# ---- Distances ----

btwn_tmp  <- list()
certain_tmp <- list()
for (i in seq_along(rdatas)) {
    load(rdatas[[i]])
    # between each inferred trees
    btwn_tmp[[i]]  <- data.frame(d.rf   = dist.list$d.rf,
                            d.wrf   = dist.list$d.wrf,
                            d.kf   = dist.list$d.kf,
                            d.sh   = dist.list$d.sh,
                            #d.kern = dist.list$d.kern,
                            prmset = as.numeric(dist.list$prmset))
    # benchmark (difference from "true" tree)
    certain_tmp[[i]] <- data.frame(d.rf   = dist.list$d.rf.star,
                            d.wrf   = dist.list$d.wrf.star,
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


if(do_tn93){
	# The raw TN93 distances:
	df_tn93 <- dist.tn93()

	saveRDS(df_tn93, file = "df_tn93.RDS")

	# Number of clusters based on TN93 distance:
	df_tn93_sub <- df_tn93 %>%
		group_by(prmset) %>% 
		slice(x = sample(1:n(), size = n()/8, replace = TRUE)) %>%
		ungroup()

	# Find the correct order(s) of magnitude
	df_tmp <- list()
	clust_tmp <- list()
	for (i in 1:5) {
		threshi <- 10^-i
		cat("Testing threshold ", threshi, "\n")
		df_tmp[[i]] <- df_tn93_sub %>%
		    clstr_num(dist.thresh.mean = threshi)
		clust_tmp[[i]] <- max(df_tmp[[i]]$n.clusters)
		cat("Max cluster number: ", clust_tmp[[i]], "\n")
	}

	print(unlist(clust_tmp))

	thresh <- c(0.02, 0.30)

	df_tn93_clustr_1 <- df_tn93 %>%
		    clstr_num(dist.thresh.mean = thresh[1])

	df_tn93_clustr_2 <- df_tn93 %>%
	    clstr_num(dist.thresh.mean = thresh[2])

}