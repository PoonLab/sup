#!/bin/sh

./clean.sh


# Generate the true phylogeny.
# Sequences saved in "seqs/sim.fasta"
# Tree saved in "trees/sim.nwk"
Rscript simul-seq.R

# Tree infered when considering
# the sequence as certain:
RAXML_OPT="-m GTRGAMMA -p 12345 --JC69 --silent" # <-- make sure it's the same as in `infr-tree-raxml.sh`
raxmlHPC $RAXML_OPT -s seqs/sim.fasta -n tree-raxml-certain.out > trees/certain.out
mv RAxML_*certain*.out* trees/

# Tree inference for 
# probabilistic sequences:
N_PRM_SET=$(wc -l < prm-btshp.csv)

pids=""
for i in $(seq 1 $N_PRM_SET)  
do
  ./run-unit.sh $i &
  pids="$pids $!"
done

# Run analysis once all is done:
wait $pids
Rscript analysis.R

echo " * * run-multi.sh completed * * "

