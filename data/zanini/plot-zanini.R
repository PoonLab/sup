source('utils-zanini.R')

patient.vec   <- c(1,3,9,11)
timepoint.vec <- c(1,3,5,0)

for(patient in patient.vec){
    for(timepoint in timepoint.vec){
        print(paste('patient:',patient,'timepoint:',timepoint))
        dat <- read_data(patient, timepoint)
        plots(dat)        
    }
}