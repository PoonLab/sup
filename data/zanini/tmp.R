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




digest <- function(dat) {
    df <- dat %>%
        mutate(total = A+C+T+G+d+N)
    
    dfp <- (df / df$total) %>%
        select(-total) %>%
        mutate(position=1:nrow(df))
    
    dfpb <- dfp %>%
        gather(key='base', value='p', -position) 
    
    return(list(df=df, dfp=dfp, dfpb=dfpb))
}


calc_entropy_one <- function(p){
    res <- NA
    if(p==0) res <- 0
    if(p>0) res <- -sum(p * log(p,base = 2))
    return(res)
}

calc_entropy <- function(dat) {
    z <- digest(dat)
    dfpb = z$dfpb
    
    a <- dfpb %>%
        group_by(position) %>%
        summarise(entropy = calc_entropy_one(p))
    
    # STOPPED HERE
    # GTG 
}


plots <- function(dat) {
    z <- digest(dat)
    df = z$df
    dfpb = z$dfpb
    fname2 <- attr(dat, 'fname2')
    
    # Coverage:
    g.cvg <- df %>%
        mutate(position = 1:nrow(df)) %>%
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
    
    
    g.histp <- dfpb %>%
        ggplot(aes(x=p)) +
        geom_histogram(aes(fill=base), bins=10)+
        facet_wrap(~base, scales = 'fixed') + 
        scale_y_log10()+
        ggtitle(fname2)
    
    fplot <- paste0('plot-',fname2,'.pdf')
    pdf(file = fplot, width=12, height = 10)
    plot(g.histp)
    plot(g.cvg)
    # plot(g.pviz)
    dev.off()
}

# --- RUN ----

patient.vec   <- c(1,3,9,11)
timepoint.vec <- c(1,5,0)

for(patient in patient.vec){
    for(timepoint in timepoint.vec){
        print(paste('patient:',patient,'timepoint:',timepoint))
        dat <- read_data(patient, timepoint)
        plots(dat)        
    }
}


