#!/bin/sh

./clean.sh
./sim-all-reads.sh
Rscript read-fastq-proba.R

echo DONE.
