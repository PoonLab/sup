source('utils-zanini.R')

patient.vec   <- c(1:11)
timepoint.vec <- c(1:5,0)

dat.list <- list() ; k=1
for(patient in patient.vec){
    for(timepoint in timepoint.vec){
        print(paste('Entropy :: patient:',patient,'timepoint:',timepoint))
        tmp <- read_data(patient, timepoint) %>%
            digest()
        dat.list[[k]] <- tmp$dfp %>%
            select(position, entropy)
        dat.list[[k]]$patient <- patient
        dat.list[[k]]$tp <- timepoint
        k = k+1
    }
}

zanini.entropy <- do.call('rbind', dat.list)
save(list = 'zanini.entropy', 
     file = 'zanini-entropy.RData')
