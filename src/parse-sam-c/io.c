#include "io.h"

void read_input_file(char *file, in_file *in, int n) {
    FILE *fp = fopen(file, "r");
    if (fp == NULL)
    {
        perror(file);
        exit(1);
    }

    char *line = NULL;
    size_t line_buf_size = 0;
    ssize_t line_size;

    char *tok;
    const char delim[2] = "\t";
    int i = 0;
    while((line_size = getline(&line, &line_buf_size, fp)) != -1)
    {
        if (line[0] != '@') {
            tok = strtok(line, delim);

            // Assuming that each line in the file has the same number of tokens
            in[i].qname = malloc(strlen(tok) + 1);
            strcpy(in[i].qname, tok);
            
            tok = strtok(NULL, delim);
            in[i].flag = malloc(strlen(tok) + 1);
            strcpy(in[i].flag, tok);

            tok = strtok(NULL, delim);
            in[i].rname = malloc(strlen(tok) + 1);
            strcpy(in[i].rname, tok);

            tok = strtok(NULL, delim);
            in[i].pos = malloc(strlen(tok) + 1);
            strcpy(in[i].pos, tok);

            tok = strtok(NULL, delim);
            in[i].mapq = malloc(strlen(tok) + 1);
            strcpy(in[i].mapq, tok);

            tok = strtok(NULL, delim);
            in[i].cigar = malloc(strlen(tok) + 1);
            strcpy(in[i].cigar, tok);

            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);

            tok = strtok(NULL, delim);
            in[i].seq = malloc(strlen(tok) + 1);
            strcpy(in[i].seq, tok);

            tok = strtok(NULL, delim);
            in[i].qual = malloc(strlen(tok) + 1);
            strcpy(in[i].qual, tok);

            i++;
        }
    }
    fclose(fp);
    free(line);
}

int get_number_rows(char *file) {
    FILE *fp = fopen(file, "r");
    if (fp == NULL)
    {
        exit(1);
    }

    char *line = NULL;
    size_t line_buf_size = 0;
    ssize_t line_size;
    int num_lines = 0;

    while((line_size = getline(&line, &line_buf_size, fp)) != -1)
    {
        if (line[0] != '@') {
            num_lines++;
        }
    }
    fclose(fp);
    free(line);
    return num_lines;
}

void write_matrix(double** m, int maxLen, char *filename) {
    char *tok;
    const char delim[2] = ".";
    tok = strtok(basename(filename), delim);

    char outputFile[strlen(tok) + 5];
    snprintf(outputFile, sizeof(outputFile), "%s.csv", tok);

    FILE *fp = fopen(outputFile, "w+");
    fprintf(fp, " , A, C, G, T\n");
    for (int i = 0; i < maxLen; i++) {
        fprintf(fp, "%d, %f, %f, %f, %f\n", (i+1), m[i][0], m[i][1], m[i][2], m[i][3]);
    }
    fclose(fp);
}

