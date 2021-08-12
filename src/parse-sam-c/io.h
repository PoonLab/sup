#ifndef IO_H
#define IO_H
#define  _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include "types.h"

void read_input_file(char *file, in_file *in, int n);
int get_number_rows(char *file);
void write_matrix(double ** m, int maxLen, char *filename);
void write_insertions(insertion_info *insertions, char *filename);

#endif