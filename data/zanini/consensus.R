###
###   Consensus sequence for all patients
###



library(seqinr)
source('utils-zanini.R')


#' Extract the most frequent nucleotide.
most_frequent <- function(dat) {
    df <- select(dat, A,C,G,T) 
    mf <- apply(df,1,which.max)
    a <- c('A','C','G','T')
    seq <- a[mf]
    return(seq)
}

# ---- RUN ----

patient.vec   <- c(1,3,9,11)
timepoint.vec <- c(1,3,5,0)

seq.patient.time <- list() ; k = 1

for(patient in patient.vec){
    for(timepoint in timepoint.vec){
        print(paste('patient:',patient,'timepoint:',timepoint))
        dat <- read_data(patient, timepoint)
        
        # Consensus = most frequent:
        seq.k <- most_frequent(dat)
        
        # Save FASTA file:
        id <- paste('patient', patient, 'TP', timepoint, sep='_')
        write.fasta(sequences = seq.k, 
                    names     = id, 
                    file.out  = paste0(id,'_consensus.fasta'))
        
        # Store in list (not used yet)
        seq.patient.time[[k]] <- list(patient = patient, 
                                      timepoint = timepoint,
                                      seq = seq.k)
        k=k+1
    }
}




