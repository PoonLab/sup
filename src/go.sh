#!/bin/sh

./clean.sh

Rscript simul-seq.R
Rscript draw-tip-seqs.R
./calc-tree.sh
# Rscript ss.R

