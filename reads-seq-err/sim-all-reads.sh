#!/bin/sh

NREADS=5000

./sim-reads.sh hiseq $NREADS
./sim-reads.sh miseq $NREADS
./sim-reads.sh novaseq $NREADS

echo "Reads simulation completed."



