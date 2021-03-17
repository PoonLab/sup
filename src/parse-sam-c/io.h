#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include "types.h"

void read_input_file(char *file, in_file *in, int n);
int get_number_rows(char *file);
void write_matrix(double m[][4], int maxLen);

#endif