#!/bin/sh

./clean.sh


# Generate the true phylogeny.
# Sequences saved in "seqs/sim.fasta"
# Tree saved in "trees/sim.nwk"
Rscript simul-seq.R

N_PRM_SET=$(wc -l < prm-btshp.csv)

for i in $(seq 1 $N_PRM_SET)  
do
./run-unit.sh $i &
done

# Analyze the distances:
# Rscript analysis.R

echo " * * run-multi.sh completed * * "

