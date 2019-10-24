## 
##  Simulate evolution of a 
##  "root" nucleotide sequence 
##

library(phylosim, quietly = TRUE, warn.conflicts = F)
library(dplyr,quietly = TRUE, warn.conflicts = F)
source('utils.R')

set.seed(1234)

message("\nStarting phylogeny simulation...")

prm <- read.csv('prm.csv')

# Define the evolution process:
ev.proc    <- phylosim::JC69()  
summary(ev.proc)

# Define the "root" sequence from
# which evolution occurs:
root.seq.length <- get_prm(prm,'phylosim.root.seq.length')
root.seq <- NucleotideSequence(length = root.seq.length, 
                               processes = list(list(ev.proc)) ) %>%
    sampleStates()
print("The root sequence is:")
print(root.seq)
print(paste("Root sequence length:", root.seq.length))

# Draw invariable positions 
invar.n   <- get_prm(prm,'phylosim.n.invar')
invar.pos <- sample(1:root.seq.length, invar.n)
print("Invariant positions:")
print(sort(invar.pos))
setRateMultipliers(root.seq,ev.proc,0,invar.pos)
# getRateMultipliers(root.seq,ev.proc)


# Define the evolution tree
# and how spreaded it is:
n.tips   <-  get_prm(prm,'phylosim.n.tips')
tree.sim <- rcoal(n.tips)
message(paste('Tree with',n.tips,'tips built.'))

# Simulate evolution of 
# the root seq on the tree:
sim <- PhyloSim(phy  = tree.sim,
                root = root.seq )
Simulate(sim)

# Save
fname <- 'seqs/sim.fasta'
phylosim::saveAlignment(this = sim, 
                        file = fname, 
                        skip.internal = TRUE, 
                        paranoid = TRUE)
# print(sim$alignment)

message(paste("Simulated phylogeny saved in:",
              fname))

pdf('plot-sim-phylo.pdf', 
    width=15, height = 10)
plot(ev.proc)
plot(tree.sim, type = 'clado')
plot(sim, num.pages = 1)
dev.off()

message("Phylogeny simulation done.")
