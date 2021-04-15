
The goal is to simulate sequencing reads with `InSilicoSeq` and visualize the error probabilities of the base calls.
The simulations mimic *Illumina* instruments: HiSeq, MiSeq and NovaSeq. 

To install `InSilicoSeq`, run the script `install-iss.sh`. 

Run `go.sh` for the full execution:

 1. simulate reads from a "true" source sequence
 2. translate FASTQ quality scores into error probability (for each position in each read)
 3. plot the mean error probability, by position (`plot-errorProba.pdf`)


