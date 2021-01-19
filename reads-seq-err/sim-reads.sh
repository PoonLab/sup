#!/bin/sh

# Use the software "InSilicoSeq" to generate
# short reads from a source DNA sample 
# as an Illumina instrument would produce.
# 
# Source:
# https://insilicoseq.readthedocs.io/en/latest/
# https://github.com/HadrienG/InSilicoSeq

# miseq , hiseq , novaseq <-- Type of Illumina instrument  
MODELSEQ=$1            
# Number of reads generated- default is 1000000
NREADS=$2               
# FASTA file that contains the "true" source sample
FASTAFILE=../data/seqs/hiv-db.fasta   
# output files generic name:
OUTFILES=$MODELSEQ\_reads
# Number of cores used:
NCPUS=4

iss generate --genomes $FASTAFILE --model $MODELSEQ --output $OUTFILES --n_reads $NREADS --cpus $NCPUS --gc_bias


# hiv-db.fasta: HIV sequences (full genome) randomly picked from https://www.hiv.lanl.gov/content/sequence/HIV/mainpage.html
