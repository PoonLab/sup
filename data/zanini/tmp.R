library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())

# Download data at: https://hiv.biozentrum.unibas.ch/data/




read_data <- function(patient, timepoint) {
    flist <- system(paste0('ls act_p', patient), intern = TRUE) 
    flist <- flist[!grepl('README', flist)]
    pos   <- gregexpr('[0-9]+', flist)
    
    smpld <- numeric(length = length(pos))
    for(i in seq_along(pos)){
        smpld[i] <- substr(flist[i],
                           start = pos[[i]][1], 
                           stop = attr(pos[[i]],'match.length')) %>%
            as.numeric()
    }
    smpld.sort <- sort(smpld)
    flist.sort <- flist[order(smpld)]
    
    if(timepoint==0) timepoint <- length(flist)
    
    fname <- paste0('act_p',patient,'/',
                    flist.sort[timepoint])
    
    fname2 <- paste('patient', patient,
                    'day', smpld.sort[timepoint], 
                    sep='_')
    
    dat <- read.table(file = fname, 
                      sep = '\t', 
                      stringsAsFactors = F , 
                      header = F, 
                      comment.char = '#')
    
    headers <- c('A','C','G','T','d','N')
    names(dat) <- headers
    attr(dat,'fname2') <- fname2
    return(dat)
}



calc_entropy_one <- function(p){
    res <- NA
    if(p==0) res <- 0
    if(p>0) res <- -sum(p * log(p,base = 2))
    return(res)
}


calc_entropy <- function(dfp) {
    # Calculate "p*log(p)" for 
    # each position and base:
    tmp <- list()
    for(b in c('A','C','G','T','d','N')){
        tmp[[b]] <- sapply(dfp[[b]], calc_entropy_one)
    }
    
    # Entropy for each position:
    m <- do.call('cbind', tmp)
    dfp$entropy <- apply(m, MARGIN = 1, FUN=sum)    
    return(dfp)
}


digest <- function(dat) {
    
    df <- dat %>%
        mutate(total = rowSums(.[1:ncol(dat)]))
    
    dfp <- (df / df$total) %>%
        select(-total) %>%
        mutate(position = 1:nrow(df))
    
    dfpb <- dfp %>%
        gather(key='base', value='p', -position) 
    
    # Add after dfpb creation to avoid messing up (?)
    dfp <- calc_entropy(dfp)
    
    return(list(df=df, dfp=dfp, dfpb=dfpb))
}


proba_by_base <- function(dat) {
    x <- digest(dat)
    dfpb <- x$dfpb
    
    y <- dfpb %>%
        filter(base=='A')
    
    pr <- y$p
    mean(pr)
    pr[pr==0] <- 1e-10
    pr[pr==1] <- 1-1e-10
    
    
    nllk <- function(x, pr) {
        return(sum(dbeta(x = pr, 
                  shape1 = x[1], 
                  shape2 = x[2], 
                  log = TRUE)))
    }
    a <- optim(par = c(0.5,0.5), fn = nllk, pr=pr)
    a.fit <- a$par
    a.fit
    a.fit[1]/sum(a.fit)
    
    xx <- seq(0,1,length.out = 100)
    yy <- dbeta(xx, a.fit[1], a.fit[2])
    plot(xx,yy,typ='l')
    
}

plots <- function(dat) {
    z    <- digest(dat)
    df   <- z$df
    dfp  <- z$dfp
    dfpb <- z$dfpb
    
    fname2 <- attr(dat, 'fname2')
    
    # Coverage:
    dfc <- df %>%
        mutate(position = 1:nrow(df))
    g.cvg <- dfc %>%
        ggplot(aes(x=position, y=total)) + 
        geom_line() +
        scale_y_log10()+
        ggtitle(fname2)
    
    # Visualize probabilities for each position:
    g.pviz <- dfpb %>%
        filter(position < 50) %>%
        ggplot(aes(x=position)) +
        geom_bar(aes(y=p, fill=base), 
                 colour='grey',
                 stat = 'identity')+
        ggtitle(fname2)
    
    # Histogram of probabilities by base:
    g.histp <- dfpb %>%
        ggplot(aes(x=p)) +
        geom_histogram(aes(fill=base), bins=10)+
        facet_wrap(~base, scales = 'fixed') + 
        scale_y_log10()+
        ggtitle(fname2)
    
    # Histogram of entropy by position:
    g.histe <- dfp %>%
        ggplot(aes(x=entropy))+
        geom_histogram(binwidth = 0.05)+
        scale_y_log10()+
        ggtitle(fname2)
    g.histe
    
    # Entropy by position
    bucket <- 50
    g.pose <- dfp %>%
        mutate(roundpos = round(position/bucket)*bucket) %>%
        group_by(roundpos) %>%
        summarise(m.entropy = mean(entropy),
                  min.entropy = min(entropy),
                  max.entropy = max(entropy)) %>%
        ggplot(aes(x=roundpos)) + 
        geom_point(aes(y=m.entropy), size=0.3) + 
        geom_line(aes(y=m.entropy), size=0.3) + 
        geom_ribbon(aes(ymin=min.entropy, 
                        ymax=max.entropy),
                        alpha=0.4)+
        ggtitle(fname2, 
                paste('Bucket size =',bucket, 
                      '; Total entropy =', round(sum(dfp$entropy),3)))
    # g.pose
    
    # Entropy & coverage
    
    g.pose2 <- g.pose + 
        geom_line(data = dfc, 
                  aes(x=position, y=log(total)/max(log(total))),
                  colour = 'blue', alpha=0.8)
    
    dfj <- data.frame(entropy=dfp$entropy, 
                      cvg = dfc$total,
                      position = dfc$position)
    
    g.cvgent <- dfj %>% 
        ggplot()+
        geom_point(aes(y=entropy, 
                       x=cvg,
                       colour = position),
                   alpha = 0.25, size=3)+
        scale_x_log10()+
        scale_y_log10()
    
    
    fplot <- paste0('plot-',fname2,'.pdf')
    pdf(file = fplot, width=12, height = 10)
    plot(g.pose2)
    plot(g.cvgent)
    plot(g.histp)
    plot(g.cvg)
    # plot(g.pviz)
    dev.off()
}

# --- RUN ----

patient.vec   <- 1 # c(1,3,9,11)
timepoint.vec <- 0 # c(1,3,5,0)

for(patient in patient.vec){
    for(timepoint in timepoint.vec){
        print(paste('patient:',patient,'timepoint:',timepoint))
        dat <- read_data(patient, timepoint)
        plots(dat)        
    }
}


