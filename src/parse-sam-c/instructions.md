# Dependencies

glib - GHashTable from https://developer.gnome.org/

# Running the program
`make` will compile the program and build output to `sam2main`

`./sam2main [INPUT FILENAME]` will run the program with the SAM filename provided as a command line argument

`make clean` will remove `sam2main` and `matrix.csv` (stores the matrix as a table)

# Notes
Program has been tested on Ubuntu 16.04.7

# TODOs
- Checking for paired neighbours
- Multithreading to run the `apply_cigar()` function concurrently
