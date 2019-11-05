###
###   Estimate mean probabilities
###   from Zanini's data for a 
###   multinomial distribution 
###   (typically used in the probabilistic sequence)
###

source('utils-zanini.R')

# Select patients:
pv <- c(1:11)
# From paper: "patients p4 and p7 were excluded from analysis 
# because of sus- pected superinfection and failure 
# to amplify early samples with low virus levels, respectively"
pv <- pv[-c(4,7)]


timepoint = 2

ldat <- lapply(pv, FUN = read_data, timepoint=timepoint) %>%
    lapply(digest) %>%
    lapply('[[','dfpb')

for(i in seq_along(ldat)){
    ldat[[i]]$patient = pv[i]
    ldat[[i]]$timepoint = timepoint
}
df <- do.call('rbind',ldat)


df.probas <- df %>%
    filter(base %in% c('A','C','G','T')) %>%  # include "d" ???
    group_by(position, base) %>%
    summarize(m = mean(p, na.rm=TRUE), n=n())


prob.seq.zanini <- df.probas %>%
    select(-n) %>%
    pivot_wider(names_from = base, values_from = m)

save(list = 'prob.seq.zanini', 
     file = 'prob_seq_zanini.RData')
