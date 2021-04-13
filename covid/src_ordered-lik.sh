#!/bin/sh

rm data/pangordlineages/*

# I'm bad at arg parsing
N=5
overwrite=false
num_args=$#
if [ $num_args -gt 0 ]; then
    for (( i=1; i <= $#; i++)); do
        arg=${!i}
        if [ $arg = -N ]; then
            arg2_pos=$((i+1))
            N=${!arg2_pos}
        fi

        if [ $arg = --overwrite ]; then
            overwrite=true
        fi
    done
fi



# Create files
if [ $overwrite = true ]; then
    Rscript covid/ord_S.R -N $N --overwrite
else 
    Rscript covid/ord_S.R -N $N
fi



source ~/anaconda3/etc/profile.d/conda.sh
conda init
conda activate pangolin 


accs=`ls data/ord_covid/*.fasta`
for acc in $accs
do
    sampled_fasta=$acc
    arg_chars=${#sampled_fasta}
    echo $arg_chars
    newacc1=${acc/"data/ord_covid/"/}
    newacc2=${newacc1/"_ord.fasta"/}
    echo $newacc2

    out_fasta="data/pangordlineages/"$newacc2"_pangolineages.csv"
    pangolin $sampled_fasta --outfile $out_fasta
done


rm data/ord_covid/*

Rscript -e "rmarkdown::render('figures/ord-results.Rmd')"
