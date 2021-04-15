## 
##  Simulate evolution of a 
##  "root" nucleotide sequence 
##

suppressPackageStartupMessages(
    {library(phylosim);
    library(dplyr);
    library(seqinr)}
)

source('utils.R')

set.seed(1234)
t1 <- as.numeric(Sys.time())

message("\nStarting phylogeny simulation...")

prm <- read.csv('prm.csv')

# Define the evolution process:
ev.proc    <- phylosim::JC69()  
summary(ev.proc)

# Define the "root" sequence from
# which evolution occurs. 
# Either generate randomly or
# read from an existing sequence.

random.root <- get_prm(prm,'phylosim.root.random')

message(paste("Random root sequence:",random.root))
message('Building root sequence... ')

if(random.root){
    root.seq.length <- get_prm(prm,'phylosim.root.seq.length')
    root.seq <- NucleotideSequence(length = root.seq.length, 
                                   processes = list(list(ev.proc)) ) %>%
        sampleStates()
}

if(!random.root){
    rs <- read.fasta(file = 'seqs/patient_1_TP_3_consensus.fasta',
                     as.string = TRUE, 
                     forceDNAtolower = FALSE) %>% 
        as.character()
    root.seq <- NucleotideSequence(string = rs, 
                                   processes = list(list(ev.proc)) )
    root.seq.length <- nchar(rs)
}

message("\nThe root sequence is (first 200 nucleotides):")
message(substr(root.seq,1,200))
message(paste("Root sequence length:", root.seq.length))

# Evolution rate (for all positions):
evol.rate   <- get_prm(prm,'phylosim.evol.rate')
setRateMultipliers(root.seq, ev.proc, evol.rate)
message(paste("Evolution rate:", evol.rate))

# Draw invariable positions 
invar.p   <- get_prm(prm,'phylosim.prop.invar')
invar.n   <- round(invar.p*root.seq.length)
invar.pos <- sample(1:root.seq.length, invar.n)
message(paste("Number of invariant positions:",invar.n, '(',root.seq.length,')'))
setRateMultipliers(root.seq,ev.proc,0,invar.pos)
# getRateMultipliers(root.seq,ev.proc)


# Define the evolution tree
# and how spreaded it is:
n.tips   <- get_prm(prm,'phylosim.n.tips')
tree.sim <- ape::rcoal(n.tips)
tree.sim$tip.label <- paste('seq', 1:n.tips, sep='_')
message(paste('Tree with',n.tips,'tips built.'))


# Simulate evolution of 
# the root seq on the tree:
message('Simulating phylogeny...')
sim <- PhyloSim(phy  = tree.sim,
                root = root.seq )
Simulate(sim)
message('Simulation done.')

# Save
message('Saving...', appendLF = FALSE)
fname.seqs <- 'seqs/sim.fasta'
fname.tree <- 'trees/sim.nwk'
ape::write.tree(tree.sim, file = fname.tree)
phylosim::saveAlignment(this = sim, 
                        file = fname.seqs, 
                        skip.internal = TRUE, 
                        paranoid = TRUE)

message(paste(" done.\nSimulated phylogeny saved in:",
              fname.seqs, fname.tree))

pdf('plot-sim-phylo.pdf', 
    width=15, height = 10)
plot(ev.proc)
plot(tree.sim, main='Benchmark tree simulated')
plot(tree.sim, main='Benchmark tree simulated\nwith branch lengths')
edgelabels(round(tree.sim$edge.length,2), 
           bg="black", col="white", font=2, cex=.6)
plot(sim, num.pages = 1)
dev.off()

t2 <- as.numeric(Sys.time())
dt <- round((t2-t1)/60, 1)

message(paste("Phylogeny simulation done in",dt,"minutes."))
