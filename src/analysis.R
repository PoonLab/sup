library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())


rdatas <- system('ls *treedist*RData', intern = TRUE)


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



g <- df.d %>%
    ggplot()+
    geom_density(aes(x=d, fill=prmset, colour=prmset), 
                 alpha=0.2)+
    facet_wrap(~prmset, ncol=1)
plot(g)

gs <- df.ds %>%
    ggplot()+
    geom_density(aes(x=d.star, fill=prmset, colour=prmset), 
                 alpha=0.2)+
    facet_wrap(~prmset, ncol=1)
plot(gs)

gs2 <- df.ds %>%
    ggplot() +
    geom_violin(aes(x=prmset, y=d.star),fill='azure2')
plot(gs2)
