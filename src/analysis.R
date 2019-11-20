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
    geom_pointrange(aes(x=prmset, y=m, ymin=m-2, ymax=m+s,
                        colour=prmset),
                    size=1)+
    ggtitle('RF distance between inferred trees (mean +/- sd)')+
    ylab('RF distance')+
    guides(colour=FALSE)+
    theme(panel.grid.major.x = element_blank())

g.d <- df.d %>%
    ggplot()+
    geom_histogram(aes(x=d, fill=prmset), 
                  binwidth = 2)+
    facet_wrap(~prmset, ncol=1, scales = 'free_y')+
    ggtitle('RF distance between inferred tree ')+
    xlab('RF distance')+
    theme(panel.grid = element_blank())


# Distance from true tree:
gs <- df.ds %>%
    ggplot()+
    geom_histogram(aes(x=d.star, fill=prmset), 
                 binwidth = 2)+
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

pdf('plot-analysis.pdf', width = 8 , height = 6)
grid.arrange(g.d.ms, gs.ms, ncol=1)
plot(g.d)
plot(gs)
plot(gs2)
dev.off()

message('Analysis completed.')
