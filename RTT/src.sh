# Add worldwide covid case counts as columns
# Creates sequences_descr_wk.csv
Rscript make_descr_wk.R

# Sample uniformly from weeks
Rscript sample-seqs_unif.R -N 10

# Multiple sequence alignment (afa = Aligned FAsta)
# http://www.drive5.com/muscle/downloads.htm
# alias muscle=/usr/bin/muscle3.8.31_i86linux64
nohup muscle -in sampled_sequences.fasta -out sampled_sequences.afa > muscle.out &

# Construct Phylogeny
# FastTree 2.1.11
# http://www.microbesonline.org/fasttree/
./FastTree sampled_sequences.afa > raw_seq.nwk

# Sample uncertainty
# To save time, sample new sequences from ALIGNED sequences



