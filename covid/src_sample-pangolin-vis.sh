#! bin/sh

# Suggested usage:
# /usr/bin/time -v --outfilebash covid/src_sample-pangolin-vis.sh -d --overwrite -N 5000; /usr/bin/time -v bash covid/src_sample-pangolin-vis.sh -d --overwrite -N 5000


# Arguments
    # -N 10
        # Number of samples to take
    # -d
        # Dirichlet prior for sampling
    # --overwrite
        # Re-sample existing files

# I'm bad at arg parsing
N=10
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
if [ $dirich = true ]; then # if dirichlet sampling
    if [ $overwrite = true ]; then # if overwrite
        Rscript covid/sample_S_dir.R -d -N $N --overwrite
    else 
        Rscript covid/sample_S_dir.R -d -N $N
    fi
else 
    if [ $overwrite = true ]; then # if overwrite
        Rscript covid/sample_S_dir.R -N $N --overwrite
    else 
        Rscript covid/sample_S_dir.R -N $N
    fi
fi


# Initialize Conda
source ~/anaconda3/etc/profile.d/conda.sh
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
    
    if [ $overwrite = true ]; then # if overwrite
        pangolin $sampled_fasta --outfile $out_fasta
    else
        if [ -f "$out_fasta" ]; then
            echo "$out_fasta already exists"
        else 
            pangolin $sampled_fasta --outfile $out_fasta
        fi
    fi
done

# Visualize!
if [ $dirich = true ]; then
    Rscript -e ".libPaths('/home/devan/R/x86_64-pc-linux-gnu-library/4.0'); rmarkdown::render('figures/pangolin_results_report.Rmd', params = list(dirich=TRUE), output_file='pangolin_results_report_d.pdf')"
else 
    Rscript -e ".libPaths('/home/devan/R/x86_64-pc-linux-gnu-library/4.0'); rmarkdown::render('figures/pangolin_results_report.Rmd', params = list(dirich=FALSE))"
fi

