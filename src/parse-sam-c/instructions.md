# Dependencies

glib - GHashTable from https://developer.gnome.org/

# Running the program
`make` will compile the program and build output to `sam2aln`

`./sam2aln -f [INPUT FILENAME] -t [NUM_THREADS]` will run the program with the SAM filename and number of threads provided as a command line argument

`make clean` will remove `sam2aln` and the generated CSV file

# Notes
Program has been tested on Ubuntu 16.04.7
