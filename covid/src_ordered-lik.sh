#!/bin/sh

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

if [ $overwrite = true]; then
    rm data/pangordlineages/*
    rm data/ord_covid/*
fi


# Create files
if [ $overwrite = true ]; then
    Rscript covid/ordered-lik.R -N $N --overwrite
else 
    Rscript covid/ordered-lik.R -N $N
fi



source ~/miniconda3/etc/profile.d/conda.sh
conda init
conda activate pangolin


accs=(data/ord_covid/*.fasta)
for acc in ${accs[@]}
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


Rscript -e "rmarkdown::render('figures/ord-results.Rmd')"
