figlet "grep>"
# Get the fasta description line (with line numbers)
grep -n ">" sequences.fasta > sequences_descr_raw.txt




figlet "make_descr"
# Add worldwide covid case counts as columns
	# Creates sequences_descr_mt.csv
Rscript make_descr_mt.R




figlet "sample-seqs"
# Sample uniformly from months
	# Creates sampled_seqs.txt and sampled_SRA.txt
Rscript sample-seqs_unif.R -N 100





figlet "download SRA"
# Checks if file has been downloaded, then proceeds if not
# Adds hours to the runtime if there are 
# a bunch of new things to download
# apt install sra-toolkit
Rscript SRA_downloader.R
Rscript quality-by-date.R # Moves very uncertain SAM files to samBAD/
Rscript make_sequences_descr_downloaded.R




figlet "seqtk subseq"
# Gather sampled sequences
seqtk subseq sequences.fasta sampled_seq_dl.txt > sampled_seqs.fasta




figlet "minimap2.py"
# Align sequences
#minimap2 -a NC_045512.fa sampled_seqs.fasta > sampled_seqs_aligned.sam
sed -i "s/\ //g" sampled_seqs.fasta
python minimap2.py sampled_seqs.fasta -o sampled_seqs_aligned.fasta -a --ref NC_045512.fa
#samtools fasta sampled_seqs_aligned.sam > sampled_seqs_aligned.fa




# Clean up
#rm sampled_seqs_aligned.sam sampled_seqs.fasta sequences_descr_mt.csv
sed -i "s/,/_/g" sampled_seqs_aligned.fasta
sed -i "s/:/_/g" sampled_seqs_aligned.fasta
grep ">" sampled_seqs_aligned.fasta > sampled_seqs_aligned_descr.txt
Rscript clean_names.R sampled_seqs_aligned_descr.txt sampled_metadata.csv




figlet "treetime"
# Use tree in treetime
# https://treetime.readthedocs.io/en/latest/
# pip install phylo-treetime
treetime --dates sampled_metadata.csv --aln sampled_seqs_aligned.fasta --outdir sampled_trees/raw_tree --covariation



figlet "resample nucleotides"
# Creates a "sampled_trees" folder with subfolders
# Each subfolder has a fasta resulting from sample
Rscript sample-S-collections.R -N 50
grep -F ">" sampled_trees/sampled_tree_1.fasta > sampled_trees/sample_descr.txt
Rscript clean_names.R sampled_trees/sample_descr.txt sampled_trees/sampled_metadata.csv




figlet "Trees4samples"
# For simplicity of bash script, filenames are *.fastaligned
cd sampled_trees
samples=`ls *fasta`
for sample in $samples; do
	python ../minimap2.py $sample -o "${sample}ligned" -a --ref ../NC_045512.fa
	treetime --dates sampled_metadata.csv --aln "${sample}ligned" --covariation
done
cd ..



