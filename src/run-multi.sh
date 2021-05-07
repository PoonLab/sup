#!/bin/sh

./clean.sh

# Install missing R packages where Rscript can find them.
#Rscript check-R-packages.R

# Generate the true phylogeny.
# Parameters taken from prm.csv
#   - These params affect all of the rest of the analysis
# Sequences saved in "seqs/sim.fasta"
# Tree saved in "trees/sim.nwk"
Rscript simul-seq.R

# Tree infered when considering the sequence as certain:
RAXML_OPT="-m GTRGAMMA -p 12345 --JC69 --silent" # <-- make sure it's the same as in `infr-tree-raxml.sh`
raxmlHPC $RAXML_OPT -s seqs/sim.fasta -n tree-raxml-certain.out > trees/certain.out
mv RAxML_*certain*.out* trees/

# Tree inference for probabilistic sequences:
N_PRM_SET=$(wc -l < prm-btshp.csv)

pids=""
for i in $(seq 1 $N_PRM_SET)
do
  # Calls gen_uncertain.R, then
    # infr-tree-raxml.sh, then 
    # dist-tn93.sh, then
    # dist-calc.R
  # This is run in parallel:
  ./run-unit.sh $i &
  pids="$pids $!"
done


# Run analysis once all is done:
wait $pids
Rscript analysis.R

echo " * * run-multi.sh completed * * "

