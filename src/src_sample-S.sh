#! bin/sh

source ~/.bashrc

# Suggested usage:
# bash covid/src_sample-pangolin-vis.sh -d -N 1000


# Arguments
    # -N 10
        # Number of samples to take
    # -d
        # Dirichlet prior for sampling
    # --overwrite
        # Re-sample existing files

# I'm bad at arg parsing
N=5
dirich=false
overwrite=false
num_args=$#
if [ $num_args -gt 0 ]; then
    for (( i=1; i <= $#; i++)); do
        arg=${!i}
        if [ $arg = -N ]; then
            arg2_pos=$((i+1))
            N=${!arg2_pos}
        fi

        if [ $arg = -d ]; then
            dirich=true
        fi

        if [ $arg = --overwrite ]; then
            overwrite=true
        fi
    done
fi



# Sample from Uncertainty Matrix
# Again, I'm bad at argparsing
# I have no idea why my for loop failed,
# but it only needed to do four loops.
if [ $dirich = true ]; then # if dirichlet sampling
    Rscript src/sample-S.R -d -N $N --overwrite
    Rscript src/sample-S.R -d -N $N
    Rscript src/sample-S.R -d -N $N
    Rscript src/sample-S.R -d -N $N
    Rscript src/sample-S.R -d -N $N
else 
    Rscript src/sample-S.R -N $N --overwrite
    Rscript src/sample-S.R -N $N
    Rscript src/sample-S.R -N $N
    Rscript src/sample-S.R -N $N
    Rscript src/sample-S.R -N $N
fi




# Initialize Conda
[[ -d "~/anaconda3" ]]; source ~/anaconda3/etc/profile.d/conda.sh || source ~/miniconda3/etc/profile.d/conda.sh
conda init
conda activate pangolin 

# Find all names
asc_file="data/sampled_covid/ascnames.txt"
accs=$(cat $asc_file)

# Assign lineages with pangolin
for acc in $accs
do
    sampled_fasta=("data/sampled_covid/"$acc"_sampled.fasta")
    echo $sampled_fasta

    if [ $dirich = true ]; then 
        sampled_fasta=("data/sampled_covid/"$acc"_sampled_d.fasta")
        echo $sampled_fasta
        out_fasta=("data/pangolineages/"$acc"_pangolineages_d.csv")
    else 
        sampled_fasta=("data/sampled_covid/"$acc"_sampled.fasta")
        echo $sampled_fasta
        out_fasta=("data/pangolineages/"$acc"_pangolineages.csv")
    fi
    
    pangolin $sampled_fasta --outfile $out_fasta
done

# Visualize!
if [ $dirich = true ]; then
    Rscript -e "rmarkdown::render('results/pangolin_results_report.Rmd', params = list(dirich=TRUE), output_file='pangolin_results_report_d.pdf')"
else 
    Rscript -e "rmarkdown::render('results/pangolin_results_report.Rmd', params = list(dirich=FALSE))"
fi

