library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)

source("dist-fcts.R")
source("utils.R")

message("Analysis.R: Starting analysis...")

system("figgy=`whereis figlet`; figgylen=${#figgy}; if [ $figgylen > 10 ]; then figlet Analysis.R; fi")

# ---- Data ----

rdatas <- system("ls *treedist*RData", intern = TRUE)

# ---- Distances ----

tmp  <- list()
tmps <- list()
for (i in seq_along(rdatas)) {
    load(rdatas[[i]])
    tmp[[i]]  <- data.frame(d.rf   = dist.list$d.rf,
                            d.kf   = dist.list$d.kf,
                            d.sh   = dist.list$d.sh,
                            d.kern = dist.list$d.kern,
                            prmset = as.numeric(dist.list$prmset))
    tmps[[i]] <- data.frame(d.rf   = dist.list$d.rf.star,
                            d.kf   = dist.list$d.kf.star,
                            d.sh   = dist.list$d.sh.star,
                            d.kern = dist.list$d.kern.star,
                            prmset = as.numeric(dist.list$prmset))
    prmsimlabel <- dist.list[["prmsimlab"]]
}
df_d  <- do.call("rbind", tmp)
df_ds <- do.call("rbind", tmps)


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
                  n = sum(!is.na(value))) %>%
        mutate(cv = s / m)

    res <- left_join(tmp, df_ent, by = "prmset")

    return(res)
}



df_d_m  <- digest_distances(df_d)
df_d_ms <- digest_distances(df_ds)

# ---- TN93 distances ----

# Warning, these are distances between sequences
# in one tree, not distances between trees!

# The raw TN93 distances:
df_tn93 <- dist.tn93()

# Number of clusters based on TN93 distance:
thresh <- c(0.02, 0.30)

df_tn93_clustr_1 <- dist.tn93() %>%
    clstr_num(dist.thresh.mean = thresh[1])

df_tn93_clustr_2 <- dist.tn93() %>%
    clstr_num(dist.thresh.mean = thresh[2])

# ---- Plot Fcts ----

plot_analysis_dig <- function(dfm, subtitle = "") {

    # Mean, min, max
    g_mmm <- dfm %>%
        ggplot(aes(x = entropy, y = m)) +
        geom_point(# ymin=minv, ymax=maxv,
                   # colour=prmset,
                   aes(shape = distance.type),
                   size = 1, alpha = 0.7) +
        geom_line() +
        geom_ribbon(aes(ymin = minv, ymax = maxv), alpha = 0.2) +
        scale_x_log10() +
        facet_wrap(~distance.type, nrow = 1, scales = "fixed") +
        ggtitle("Tree distance (mean, min, max)",
                subtitle) +
        ylab("distance") +
        guides(colour = FALSE, shape = FALSE) +
        theme(panel.grid.major.x = element_blank())

    # Standard-deviation:
    g_sd <- dfm %>%
        ggplot(aes(x = entropy, y = s,
                   colour = distance.type)) +
        geom_line(size = 1) +
        geom_point(size = 3, alpha = 0.7) +
        scale_x_log10() +
        # facet_wrap(~distance.type, nrow=1)+
        ggtitle("Tree distance standard-deviation", subtitle) +
        ylab("SD distance") +
        xlab("entropy") +
        # guides(colour=FALSE)+
        theme(panel.grid.major.x = element_blank())

    # Coefficient of variation
    g_cv <- dfm %>%
        ggplot(aes(x = entropy, y = cv,
                   colour = distance.type)) +
        geom_line(size = 1) +
        geom_point(size = 3, alpha = 0.7) +
        scale_x_log10() +
        # facet_wrap(~distance.type, nrow=1)+
        ggtitle("Tree distance coeff. of variation", subtitle) +
        ylab("CV distance") +
        xlab("prmset") +
        # guides(colour=FALSE) +
        theme(panel.grid.major.x = element_blank())
    return(list(g.mmm = g_mmm,
                g.cv  = g_cv,
                g.sd  = g_sd))
}

plot_analysis <- function(df, subtitle="") {

    dfl <- pivot_longer(df, -prmset,
                        names_to = "distance.type",
                        values_to = "value")

    g_hist <- dfl %>%
        ggplot() +
        geom_histogram(aes(x = value), bins = 20) +
        facet_grid(prmset~distance.type,
                   scales = "free_y") +
        ggtitle("Tree distance", subtitle) +
        xlab("distance") +
        theme(panel.grid = element_blank())

    g_dens <- dfl %>%
        ggplot() +
        geom_density(aes(x = value, fill = prmset, color = prmset),
                     alpha = 0.5, adjust = 2) +
        facet_grid(prmset ~ distance.type,
                   scales = "free") +
        ggtitle("Tree distance", subtitle) +
        xlab("Parameter set") +
        theme(panel.grid = element_blank())

    return(list(g.hist = g_hist,
                g.dens = g_dens))
}

plot_analysis_join <- function(df_d_m, df_d_ms) {
    df_d_m$ID  <- paste(df_d_m$distance.type, df_d_m$prmset)
    df_d_ms$ID <- paste(df_d_ms$distance.type, df_d_ms$prmset)

    dfj <- left_join(df_d_m, df_d_ms, by = "ID")
    g2 <- ggplot(dfj) +
        geom_point(aes(x = m.x, y = m.y,
                       shape = distance.type.x,
                       colour = factor(prmset.x)),
                       size = 2) +
        geom_segment(aes(x = m.x, xend = m.x,
                         y = m.y - s.y, yend = m.y + s.y,
                         colour = factor(prmset.x))) +
        geom_segment(aes(x = m.x - s.x, xend = m.x + s.x,
                         y = m.y, yend = m.y,
                         colour = factor(prmset.x))) +
        geom_smooth(method = "lm",
                    aes(x = m.x, y = m.y,
                        colour = factor(prmset.x),
                        fill = factor(prmset.x)),
                    alpha = 0.1) +
        xlab("RF distance between inferred trees") +
        ylab("RF distance from Benchmark") +
        ggtitle("RF distances (mean +/- sd)") +
        geom_vline(xintercept = 0) +
        geom_hline(yintercept = 0)
    return(g2)
}


#' Plot TN93 distances within all trees
#' across all parameter sets and MC iterations
#' @param df Data frame as output by the function `dist.tn93()`
plot_tn93_distances <- function(df) {
    df_entropy <- get_entropy_prmset()

    dfs <- df %>%
        group_by(prmset) %>%
        summarise(m = mean(Distance),
                  d.min = min(Distance),
                  d.max = max(Distance),
                  d.lo = quantile(Distance, probs = 0.05),
                  d.hi = quantile(Distance, probs = 0.95)
                  ) %>%
        left_join(df_entropy, by = "prmset")


    q <- df %>%
        mutate(ps = paste0("prmset #", rmset)) %>%
        ggplot() +
        geom_histogram(aes(x = Distance),
                       binwidth = 0.02) +
        facet_wrap(~ ps, ncol = 1) +
        xlab("TN93 distance")

    g <- dfs %>%
        ggplot(aes(x = entropy, y = m)) +
        geom_line() +
        geom_pointrange(aes(ymin = d.lo,
                            ymax = d.hi),
                        size = 1) +
        scale_x_log10() +
        ylab("TN93 distance") +
        ggtitle("Raw TN93 distances")

    return(list(g.ptrng = g,
                g.hist  = q))
}


plot_clstr_num <- function(dfclst,
                           subtitle="") {
    df_entropy <- get_entropy_prmset()

    dfs <- dfclst %>%
        group_by(prmset) %>%
        summarise(m = mean(n.clusters),
                  d.min = min(n.clusters),
                  d.max = max(n.clusters),
                  d.lo = quantile(n.clusters, probs = 0.05),
                  d.hi = quantile(n.clusters, probs = 0.95)) %>%
        left_join(df_entropy, by = "prmset")

    g <- dfs %>%
        ggplot(aes(x = entropy, y = m)) +
        geom_line() +
        geom_pointrange(aes(ymin = d.min,
                            ymax = d.max)) +
        scale_x_log10() +
        ylab("Number of clusters") +
        ggtitle("Number of clusters (TN93-based)",
                subtitle)

    return(g)
}


# ---- RUN ----
# The effect of changing the first shape parameter on the
# base call distribution
pdf("plot-proba-basecall-beta.pdf")
plot_prmset_distrib("prm-btshp.csv")
dev.off()


g  <- plot_analysis(df_d, "b/w inferred trees")
g0 <- plot_analysis(df_ds, "from benchmark")

g_digest  <- plot_analysis_dig(df_d_m, "b/w inferred trees")
g_digest0 <- plot_analysis_dig(df_d_ms, "from benchmark")

g_j <- plot_analysis_join(df_d_m, df_d_ms)


fname <- paste0("plot-analysis-", prmsimlabel, ".pdf")
pdf(fname, width = 12, height = 10)
plot(g$g.hist)
#plot(g$g.dens)
plot(g_digest$g.mmm)
grid.arrange(g_digest$g.sd, g_digest$g.cv, nrow = 1)

plot(g0$g.hist)
#plot(g0$g.dens)
plot(g_digest0$g.mmm)
grid.arrange(g_digest0$g.sd, g_digest0$g.cv, nrow = 1)

plot(g_j)

g_tn93 <- plot_tn93_distances(df_tn93)

g_tn93_clustr_1 <- plot_clstr_num(df_tn93_clustr_1,
                                subtitle = paste("Threshold(mean) =",
                                                 thresh[1]))
g_tn93_clustr_2 <- plot_clstr_num(df_tn93_clustr_2,
                                subtitle = paste("Threshold(mean) =",
                                                 thresh[2]))
plot(g_tn93$g.hist)
grid.arrange(g_tn93$g.ptrng,
             g_tn93_clustr_1,
             g_tn93_clustr_2,
             ncol = 1)

dev.off()


message("Analysis completed.")
