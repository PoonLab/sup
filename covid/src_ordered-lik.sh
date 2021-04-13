#!/bin/sh

# Create files
Rscript covid/ord_S.R -N 1000


source ~/anaconda3/etc/profile.d/conda.sh
conda init
conda activate pangolin 


accs=`ls data/ord_covid/*.fasta`
for acc in $accs
do
    sampled_fasta=$acc
    newacc1=${acc/"data/ord_covid/"/}
    newacc2=${newacc1/"_ord.fasta"/}
    echo $newacc2
    out_fasta="data/pangordlineages/"$newacc2"_pangolineages.csv"
    pangolin $sampled_fasta --outfile $out_fasta
done
