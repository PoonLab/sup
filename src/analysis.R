library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)

source('dist-fcts.R')

message('Starting analysis...')

# ---- Data ----

rdatas <- system('ls *treedist*RData', intern = TRUE)

# ---- Distances ----

tmp  <- list()
tmps <- list()
for(i in seq_along(rdatas)){  # i=1
    load(rdatas[[i]])
    tmp[[i]]  <- data.frame(d.rf   = dist.list$d.rf,
                            d.sh   = dist.list$d.sh,
                            d.kern = dist.list$d.kern,
                            prmset = as.numeric(dist.list$prmset))
    tmps[[i]] <- data.frame(d.rf   = dist.list$d.rf.star, 
                            d.sh   = dist.list$d.sh.star, 
                            d.kern = dist.list$d.kern.star, 
                            prmset = as.numeric(dist.list$prmset))
    prmsimlabel <- dist.list[['prmsimlab']]
}
df.d  <- do.call('rbind',tmp)
df.ds <- do.call('rbind',tmps)


# Retrieve the sequence entropy for each prm set:
get_entropy_prmset <- function(fname = 'entropy-prmset.out') {
    x <- read.csv(fname, header = F)
    names(x) <- c('prmset', 's1','s2','entropy')
    return(x)
}

digest_distances <- function(df) {
    p.lo = 0.025
    p.hi = 0.975
    
    df.ent <- get_entropy_prmset()
    
    tmp <- df %>%
        pivot_longer(cols = -prmset, 
                     names_to = 'distance.type', 
                     values_to = 'value') %>%
        group_by(prmset, distance.type) %>%
        summarise(m = mean(value, na.rm = TRUE), 
                  md = median(value, na.rm = TRUE),
                  s = sd(value, na.rm = TRUE),
                  q.lo = quantile(value, probs = p.lo, na.rm = TRUE),
                  q.hi = quantile(value, probs = p.hi, na.rm = TRUE),
                  minv = min(value, na.rm = TRUE),
                  maxv = max(value, na.rm = TRUE),
                  n = sum(!is.na(value))) %>%
        mutate(cv = s / m)
    
    res <- left_join(tmp, df.ent, by='prmset')
    
    return(res)
}



df.d.m  <- digest_distances(df.d)
df.d.ms <- digest_distances(df.ds)

# ---- TN93 distances ----

# Warning, these are distances between sequences 
# in one tree, not distances between trees!

# The raw TN93 distances:
df.tn93 <- dist.tn93()

# Number of clusters based on TN93 distance:
thresh <- c(0.02, 0.30)

df.tn93.clustr.1 <- dist.tn93() %>%
    clstr_num(dist.thresh.mean = thresh[1]) 

df.tn93.clustr.2 <- dist.tn93() %>%
    clstr_num(dist.thresh.mean = thresh[2]) 

# ---- Plot Fcts ----

plot_analysis_dig <- function(dfm, subtitle='') {
    
    # Mean, min, max
    g.mmm <- dfm %>%
        ggplot(aes(x= entropy, y=m))+
        geom_point(    # ymin=minv, ymax=maxv,
                       # colour=prmset,
                       aes(shape = distance.type),
                   size=1, alpha=0.7)+
        geom_line()+
        geom_ribbon(aes(ymin=minv, ymax=maxv), alpha=0.2)+
        scale_x_log10()+
        facet_wrap(~distance.type, nrow=1, scales = 'fixed')+
        ggtitle('Tree distance (mean, min, max)',
                subtitle)+
        ylab('distance')+
        guides(colour=FALSE, shape=FALSE)+
        theme(panel.grid.major.x = element_blank())
    
    # Standard-deviation:
    g.sd <- dfm %>%
        ggplot(aes(x=entropy, y=s, 
                   colour = distance.type))+
        geom_line(size=1)+
        geom_point(size=3, alpha=0.7)+
        scale_x_log10()+
        # facet_wrap(~distance.type, nrow=1)+
        ggtitle('Tree distance standard-deviation', subtitle)+
        ylab('SD distance')+
        xlab('entropy')+
        # guides(colour=FALSE)+
        theme(panel.grid.major.x = element_blank())
    
    # Coefficient of variation
    g.cv <- dfm %>%
        ggplot(aes(x=entropy, y=cv, 
                   colour = distance.type))+
        geom_line(size=1)+
        geom_point(size=3, alpha=0.7)+
        scale_x_log10()+
        # facet_wrap(~distance.type, nrow=1)+
        ggtitle('Tree distance coeff. of variation', subtitle)+
        ylab('CV distance')+
        xlab('prmset')+
        # guides(colour=FALSE)+
        theme(panel.grid.major.x = element_blank())
    return(list(g.mmm = g.mmm, 
                g.cv  = g.cv,
                g.sd  = g.sd))
}

plot_analysis <- function(df, subtitle='') {
    
    dfl <- pivot_longer(df, -prmset, 
                        names_to = 'distance.type', 
                        values_to = 'value')
    
    g.hist <- dfl %>%
        ggplot()+
        geom_histogram(aes(x=value, fill=prmset), 
                       bins = 20)+
        facet_grid(prmset~distance.type, 
                   scales = 'free_y')+
        ggtitle('Tree distance', subtitle)+
        xlab('distance')+
        theme(panel.grid = element_blank())
    
    g.dens <- dfl %>%
        ggplot()+
        geom_density(aes(x=value, fill=prmset, color=prmset),
                     alpha=0.5, adjust=2)+
        facet_grid(prmset~distance.type, 
                   scales = 'free')+
        ggtitle('Tree distance', subtitle)+
        xlab('Parameter set')+
        theme(panel.grid = element_blank())
    
    return(list(g.hist=g.hist, 
                g.dens=g.dens))
}

plot_analysis_join <- function(df.d.m, df.d.ms) {
    df.d.m$ID  <- paste(df.d.m$distance.type, df.d.m$prmset)
    df.d.ms$ID <- paste(df.d.ms$distance.type, df.d.ms$prmset)
    
    dfj <- left_join(df.d.m, df.d.ms, by='ID')
    g2 <- ggplot(dfj)+
        geom_point(aes(x=m.x, y=m.y, 
                       shape = distance.type.x,
                       colour = factor(prmset.x)), 
                   size=2)+
        geom_segment(aes(x=m.x, xend=m.x, y=m.y-s.y, yend=m.y+s.y,
                         colour = factor(prmset.x)))+
        geom_segment(aes(x=m.x-s.x, xend=m.x+s.x, y=m.y, yend=m.y,
                         colour = factor(prmset.x)))+
        geom_smooth(method = 'lm', 
                    aes(x=m.x, y=m.y,
                        colour = factor(prmset.x),
                        fill = factor(prmset.x)), 
                    alpha = 0.1) +
        xlab('RF distance between inferred trees')+
        ylab("RF distance from Benchmark")+
        ggtitle("RF distances (mean +/- sd)")+
        geom_vline(xintercept = 0)+
        geom_hline(yintercept = 0)
    return(g2)
}


#' Plot TN93 distances within all trees
#' across all parameter sets and MC iterations
#' @param df Data frame as output by the function `dist.tn93()`
plot_tn93_distances <- function(df) { 
    # df = dist.tn93()
    df.entropy <- get_entropy_prmset()
    
    dfs <- df %>%
        group_by(prmset) %>%
        summarise(m = mean(Distance),
                  d.min = min(Distance),
                  d.max = max(Distance),
                  d.lo = quantile(Distance, probs=0.05),
                  d.hi = quantile(Distance, probs=0.95)
                  ) %>%
        left_join(df.entropy, by='prmset')

    
    q <- df %>%
        mutate(ps = paste0('prmset #',prmset)) %>%
        ggplot()+
        geom_histogram(aes(x=Distance),
                       binwidth = 0.02)+
        facet_wrap(~ps, ncol=1)+
        xlab('TN93 distance')
    
    g <- dfs %>%
        ggplot(aes(x=entropy, y=m))+
        geom_line()+
        geom_pointrange(aes(ymin=d.lo, 
                            ymax=d.hi),
                        size=1)+
        scale_x_log10()+
        ylab('TN93 distance')+
        ggtitle('Raw TN93 distances')
    
    return(list(g.ptrng = g,
                g.hist  = q))    
}


plot_clstr_num <- function(dfclst, 
                           subtitle='') {
    # dfclst = df.tn93.clustr
    
    df.entropy <- get_entropy_prmset()
    
    dfs <- dfclst %>%
        group_by(prmset) %>%
        summarise(m = mean(n.clusters),
                  d.min = min(n.clusters),
                  d.max = max(n.clusters),
                  d.lo = quantile(n.clusters, probs=0.05),
                  d.hi = quantile(n.clusters, probs=0.95)) %>%
        left_join(df.entropy, by='prmset')
    
    g <- dfs %>%
        ggplot(aes(x=entropy, y=m ))+
        geom_line()+
        geom_pointrange(aes(ymin=d.min, 
                            ymax=d.max))+
        scale_x_log10()+
        ylab('Number of clusters')+
        ggtitle('Number of clusters (TN93-based)',
                subtitle)
    
    return(g)
}


# ---- RUN ----

g  <- plot_analysis(df.d, 'b/w inferred trees')
g0 <- plot_analysis(df.ds, 'from benchmark')

g.digest  <- plot_analysis_dig(df.d.m, 'b/w inferred trees')
g.digest0 <- plot_analysis_dig(df.d.ms, 'from benchmark')

g.j <- plot_analysis_join(df.d.m,df.d.ms)


fname <- paste0('plot-analysis-',prmsimlabel,'.pdf')
pdf(fname, width = 12 , height = 10)
plot(g$g.hist)
#plot(g$g.dens)
plot(g.digest$g.mmm)
grid.arrange(g.digest$g.sd, g.digest$g.cv, nrow=1)

plot(g0$g.hist)
#plot(g0$g.dens)
plot(g.digest0$g.mmm)
grid.arrange(g.digest0$g.sd, g.digest0$g.cv, nrow=1)

plot(g.j)

g.tn93 <- plot_tn93_distances(df.tn93)

g.tn93.clustr.1 <- plot_clstr_num(df.tn93.clustr.1,
                                subtitle = paste('Threshold(mean) =',
                                                 thresh[1]))
g.tn93.clustr.2 <- plot_clstr_num(df.tn93.clustr.2,
                                subtitle = paste('Threshold(mean) =',
                                                 thresh[2]))
plot(g.tn93$g.hist)
grid.arrange(g.tn93$g.ptrng, 
             g.tn93.clustr.1,
             g.tn93.clustr.2,
             ncol=1)

dev.off()


message('Analysis completed.')
