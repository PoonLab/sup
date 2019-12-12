library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(gridExtra)

message('Starting analysis...')

# ---- Data ----

rdatas <- system('ls *treedist*RData', intern = TRUE)

# ---- Distances ----

tmp  <- list()
tmps <- list()
for(i in seq_along(rdatas)){  # i=1
    load(rdatas[[i]])
    tmp[[i]]  <- data.frame(d.rf = dist.list$d.rf,
                            d.sh = dist.list$d.sh,
                            prmset = dist.list$prmset)
    tmps[[i]] <- data.frame(d.rf = dist.list$d.rf.star, 
                            d.sh = dist.list$d.sh.star, 
                            prmset = dist.list$prmset)
    prmsimlabel <- dist.list[['prmsimlab']]
}
df.d  <- do.call('rbind',tmp)
df.ds <- do.call('rbind',tmps)



digest_distances <- function(df) {
    p.lo = 0.025
    p.hi = 0.975
    
    res <- df %>%
        pivot_longer(cols = -prmset, 
                     names_to = 'distance.type', 
                     values_to = 'value') %>%
        group_by(prmset, distance.type) %>%
        summarise(m = mean(value), 
                  md = median(value),
                  s = sd(value),
                  q.lo = quantile(value, probs = p.lo),
                  q.hi = quantile(value, probs = p.hi),
                  minv = min(value),
                  maxv = max(value)) %>%
        mutate(cv = s / m)
    return(res)
}


df.d.m  <- digest_distances(df.d)
df.d.ms <- digest_distances(df.ds)

# ---- Plot Fcts ----

plot_analysis_dig <- function(dfm, subtitle='') {
    
    # Mean, min, max
    g.mmm <- dfm %>%
        ggplot()+
        geom_pointrange(aes(x=prmset, 
                            y=m, ymin=minv, ymax=maxv,
                            colour=prmset,
                            shape = distance.type),
                        size=1, alpha=0.7)+
        facet_wrap(~distance.type, nrow=1)+
        ggtitle('Tree distance (mean, min, max)',
                subtitle)+
        ylab('distance')+
        guides(colour=FALSE, shape=FALSE)+
        theme(panel.grid.major.x = element_blank())
    
    # Standard-deviation:
    g.sd <- dfm %>%
        ggplot(aes(x=as.numeric(prmset), y=s, 
                   colour = distance.type))+
        geom_line(size=1)+
        geom_point(size=3, alpha=0.7)+
        # facet_wrap(~distance.type, nrow=1)+
        ggtitle('Tree distance standard-deviation', subtitle)+
        ylab('SD distance')+
        xlab('prmset')+
        # guides(colour=FALSE)+
        theme(panel.grid.major.x = element_blank())
    
    # Coefficient of variation
    g.cv <- dfm %>%
        ggplot(aes(x=as.numeric(prmset), y=cv, 
                   colour = distance.type))+
        geom_line(size=1)+
        geom_point(size=3, alpha=0.7)+
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
                       colour = prmset.x), 
                   size=2)+
        geom_segment(aes(x=m.x, xend=m.x, y=m.y-s.y, yend=m.y+s.y, colour=prmset.x))+
        geom_segment(aes(x=m.x-s.x, xend=m.x+s.x, y=m.y, yend=m.y, colour=prmset.x))+
        geom_smooth(method = 'lm', aes(x=m.x, y=m.y, 
                                       colour=distance.type.x), 
                    alpha = 0.2) +
        xlab('RF distance between inferred trees')+
        ylab("RF distance from Benchmark")+
        ggtitle("RF distances (mean +/- sd)")+
        geom_vline(xintercept = 0)+
        geom_hline(yintercept = 0)
    return(g2)
}

# ---- Plots ----

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
dev.off()


message('Analysis completed.')
