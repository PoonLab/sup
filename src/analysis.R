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
    tmp[[i]]  <- data.frame(d = dist.list$d, 
                            prmset = dist.list$prmset)
    tmps[[i]] <- data.frame(d.star = dist.list$d.star, 
                            prmset = dist.list$prmset)
}
df.d  <- do.call('rbind',tmp)
df.ds <- do.call('rbind',tmps)

df.d.ms <- df.d %>%
    group_by(prmset) %>%
    summarise(m = mean(d), s = sd(d)) 

df.ds.ms <- df.ds %>%
    group_by(prmset) %>%
    summarise(m = mean(d.star), s = sd(d.star)) 


# ---- Plots ----

# Between inferred trees

g.d.ms <- df.d.ms %>%
    ggplot()+
    geom_pointrange(aes(x=prmset, y=m, ymin=m-s, ymax=m+s,
                        colour=prmset),
                    size=1)+
    ggtitle('RF distance between inferred trees (mean +/- sd)')+
    ylab('RF distance')+
    guides(colour=FALSE)+
    theme(panel.grid.major.x = element_blank())

g.d.v <- df.d.ms %>%
    ggplot()+
    geom_point(aes(x=prmset, y=s), size=3)+
    geom_line(aes(x=as.numeric(prmset), y=s))+
    ggtitle('Variance (sd) RF distance between inferred trees') +
    ylab('RF distance SD') +
    theme(panel.grid.major.x = element_blank())

g.d <- df.d %>%
    ggplot()+
    geom_histogram(aes(x=d, fill=prmset), 
                  bins = 20)+
    facet_wrap(~prmset, ncol=1, scales = 'free_y')+
    ggtitle('RF distance between inferred tree ')+
    xlab('RF distance')+
    theme(panel.grid = element_blank())


# Distance from true tree:
gs <- df.ds %>%
    ggplot()+
    geom_histogram(aes(x=d.star, fill=prmset), 
                 bins = 20)+
    facet_wrap(~prmset, ncol=1, scales = 'free_y')+
    ggtitle('RF distance from true tree T*')+
    theme(panel.grid = element_blank())


gs2 <- df.ds %>%
    ggplot() +
    geom_boxplot(aes(x=prmset, y=d.star, fill=prmset))+
    theme(panel.grid.major.x = element_blank())+
    ggtitle('RF distance from true tree T*')


gs.ms <- df.ds.ms %>%
    ggplot(aes(x=prmset))+
    geom_pointrange(aes(y=m, ymin=m-s, ymax=m+s,
                        colour=prmset),
                    size=1)+
    theme(panel.grid.major.x = element_blank())+
    ggtitle('RF distance from true tree T* (mean +/- sd)')+
    guides(colour=FALSE)+
    ylab('RF distance')

g.ds.v <- df.ds.ms %>%
    ggplot()+
    geom_point(aes(x=prmset, y=s), size=3)+
    geom_line(aes(x=as.numeric(prmset), y=s))+
    ggtitle('Variance (sd) RF distance from true tree T*') +
    ylab('RF distance SD') +
    theme(panel.grid.major.x = element_blank())


dfj <- left_join(df.d.ms, df.ds.ms, by='prmset')
g2 <- ggplot(dfj)+
    geom_point(aes(x=m.x, y=m.y, colour = prmset), size=2)+
    geom_segment(aes(x=m.x, xend=m.x, y=m.y-s.y, yend=m.y+s.y, colour=prmset))+
    geom_segment(aes(x=m.x-s.x, xend=m.x+s.x, y=m.y, yend=m.y, colour=prmset))+
    geom_smooth(method = 'lm', aes(x=m.x, y=m.y), alpha = 0.2) +
    xlab('RF distance between Ts')+
    ylab("RF distance from T*")+
    ggtitle("RF distances (mean +/- sd)")+
    geom_vline(xintercept = 0)+
    geom_hline(yintercept = 0)
g2

pdf('plot-analysis.pdf', width = 8 , height = 6)
grid.arrange(g.d.ms, gs.ms, ncol=1)
grid.arrange(g.d.v, g.ds.v, ncol=1)
plot(g2)
plot(g.d)
plot(gs)
plot(gs2)
dev.off()

message('Analysis completed.')
