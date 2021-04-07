#include "io.h"
#include <time.h>
#include <regex.h>
#include <stdbool.h>
#include <math.h>
#include <glib.h>

int len_terminal_gap(char *seq) {
    bool prefix = true;
    int i;

    if (prefix) {
        // length of terminal gap prefix
        for (i = 0; i < strlen(seq); i++) {
            if (seq[i] != '-')
                break;
        }
        return (i > 0) ? i-1 : 0;
    }
    else {
        // length of terminal gap suffix
        for (i = strlen(seq) - 1; i >= 0; i--) {
            if (seq[i] != '-')
                break;
        }
        return ((i == strlen(seq) - 1) ? 0 : strlen(seq) - i);
    }
}

newSeq * apply_cigar(char *cigar, char *pos, char *seq, char *qual) {
    // TODO: Check if CIGAR is valid

    int num_valid = 0;
    
    /************************** Tokenize CIGAR ************************/
    // printf("Tokenize CIGAR\n");

    if (strcmp(cigar, "*") == 0) return NULL;

    tokens *first = NULL;
    tokens *prev_tok = first;

    // printf("Start Token\n");
    regex_t pexp;
    regmatch_t whole_match;

    char *sz = cigar;

    regcomp(&pexp, "[0-9]+[MIDNSHPX=]", REG_EXTENDED);

    int eflags = 0;
    int match = 0;
    size_t offset = 0;
    size_t length = strlen(sz);

    while (regexec(&pexp, sz + offset, 1, &whole_match, eflags) == 0) {
        eflags = REG_NOTBOL;
        match = 1;

        tokens *current_tok = malloc(sizeof(tokens));
        if (current_tok == NULL) printf("ERROR");

        current_tok->next = NULL;

        size_t leng = (offset + whole_match.rm_eo) - (offset + whole_match.rm_so);
        current_tok->value = malloc(leng + 1);
        if (current_tok->value == NULL) printf("ERROR");

        strncpy(current_tok->value, sz + offset + whole_match.rm_so, leng);
        (current_tok->value)[leng] = '\0';

        if (first != NULL) prev_tok->next = current_tok;
        else first = current_tok;
        
        prev_tok = current_tok;

        // increase the starting offset
        offset += whole_match.rm_eo;

        // Increase offset to prevent infinite loop (zero-length match)
        if (whole_match.rm_so == whole_match.rm_eo) {
            offset += 1;
        }

        // End of the string
        if (offset > length) {
            break;
        }
    }
    if (! match) {
        printf("\"%s\" does not contain a match\n", sz);
    }

    regfree(&pexp);
    num_valid++;

    /******************************** Edits *********************************/

    edits *edits_start = NULL;
    edits *prev_edits = edits_start;

    int seql = 1;
    int left = 1;
    int flag = 1;
    int curr_pos = atoi(pos);

    tokens *temp_tok = first;
    while(temp_tok!=NULL) {
        char *curr_tok = temp_tok->value;
        regex_t exp;
        regmatch_t num_match;

        regcomp(&exp, "[0-9]+", REG_EXTENDED);
        if (regexec(&exp, curr_tok, 1, &num_match, 0) == 0) 
        {
            char *num;
            size_t num_digits = num_match.rm_eo - num_match.rm_so;
            num = malloc(num_digits + 1);
            strncpy(num, curr_tok, num_digits);
            num[num_digits] = '\0';
            edits *current_edit = malloc(sizeof(edits));
            current_edit->len = atoi(curr_tok);
            current_edit->operand = curr_tok[strlen(curr_tok) - 1];
            current_edit->seql = seql;
            current_edit->left = left;

            // Free num 
            free(num);

            // First element have - and ! for seqPart and qualPart
            if (flag == 1) {
                current_edit->seqPart = malloc(curr_pos + 1);
                memset(current_edit->seqPart, '-', curr_pos);
                (current_edit->seqPart)[curr_pos] = '\0';

                current_edit->qualPart = malloc(curr_pos + 1);
                memset(current_edit->qualPart, '!', curr_pos);
                (current_edit->qualPart)[curr_pos] = '\0';

                flag = 0;
            }
            else if (current_edit->operand == 'D') {
                current_edit->seqPart = malloc(current_edit->len + 1);
                memset(current_edit->seqPart, '-', current_edit->len);
                (current_edit->seqPart)[current_edit->len] = '\0';

                current_edit->qualPart = malloc(current_edit->len + 1);
                memset(current_edit->qualPart, '!', current_edit->len);
                (current_edit->qualPart)[current_edit->len] = '\0';
            }
            else {
                current_edit->seqPart = NULL;
                current_edit->qualPart = NULL;
            }

            if (current_edit->operand == 'M' || current_edit->operand == 'I') {
                if (current_edit->seqPart != NULL) {
                    char *sequence = malloc(current_edit->len + 1);
                    strncpy(sequence, seq + (current_edit->left - 1), current_edit->len);
                    (sequence)[current_edit->len] = '\0';
                    int len_target = strlen(current_edit->seqPart) + current_edit->len + 1;
                    current_edit->seqPart = realloc(current_edit->seqPart,len_target);
                    strncat(current_edit->seqPart, sequence, current_edit->len + 1);
                    free(sequence);
                }
                else {
                    current_edit->seqPart = malloc(current_edit->len + 1);
                    strncpy(current_edit->seqPart, seq + (current_edit->left - 1), current_edit->len);
                    (current_edit->seqPart)[current_edit->len] = '\0';
                }

                if (current_edit->qualPart != NULL) {
                    char *qual_part = malloc(current_edit->len + 1);
                    strncpy(qual_part, qual + (current_edit->left - 1), current_edit->len);
                    (qual_part)[current_edit->len] = '\0';
                    int len_target = strlen(current_edit->qualPart) + current_edit->len + 1;
                    current_edit->qualPart = realloc(current_edit->qualPart ,len_target);
                    strncat(current_edit->qualPart, qual_part, current_edit->len + 1);
                    free(qual_part);
                }
                else {
                    current_edit->qualPart = malloc(current_edit->len + 1);
                    strncpy(current_edit->qualPart, qual + (current_edit->left - 1), current_edit->len);
                    (current_edit->qualPart)[current_edit->len] = '\0';
                }
            }

            // Update lenSeq
            curr_pos+=current_edit->len;
            current_edit->lenSeq = curr_pos;
            current_edit->next = NULL;

            // Update left
            if (current_edit->operand == 'M' || current_edit->operand == 'I' || current_edit->operand == 'S') {
                left += current_edit->len;
            }

            if (edits_start != NULL) prev_edits->next = current_edit;
            else edits_start = current_edit;

            prev_edits = current_edit;
        }
        temp_tok=temp_tok->next;
        regfree(&exp);
    }


    /****************************** Pull sequences by seqPart row ************************/

    edits *start_edits = edits_start;
    newSeq *sequences_list = NULL;
    newSeq *prev_seq = sequences_list;
    int currSeql = edits_start->seql;
    while (start_edits != NULL) {
        newSeq *curr_sequence = malloc(sizeof(newSeq));
        curr_sequence->seq = NULL;
        curr_sequence->qualseq = NULL;
        curr_sequence->next=NULL;

        while (start_edits != NULL && start_edits->seql == currSeql) {
            if (start_edits->operand != 'I' && start_edits->seqPart != NULL) {
                if (curr_sequence->seq == NULL){
                    curr_sequence->seq = malloc(strlen(start_edits->seqPart) + 1);
                    strncpy(curr_sequence->seq, start_edits->seqPart, strlen(start_edits->seqPart));
                    (curr_sequence->seq)[strlen(start_edits->seqPart)] = '\0';

                    curr_sequence->qualseq = malloc(strlen(start_edits->qualPart) + 1);
                    strncpy(curr_sequence->qualseq, start_edits->qualPart, strlen(start_edits->qualPart));
                    (curr_sequence->qualseq)[strlen(start_edits->qualPart)] = '\0';
                }
                else {
                    int curr_len = strlen(curr_sequence->seq) + strlen(start_edits->seqPart) + 1;
                    curr_sequence->seq = realloc(curr_sequence->seq, curr_len);
                    strncat(curr_sequence->seq, start_edits->seqPart, strlen(start_edits->seqPart) + 1);

                    int qual_len = strlen(curr_sequence->qualseq) + strlen(start_edits->qualPart) + 1;
                    curr_sequence->qualseq = realloc(curr_sequence->qualseq, qual_len);
                    strncat(curr_sequence->qualseq, start_edits->qualPart, strlen(start_edits->qualPart) + 1);
                }
            }

            start_edits = start_edits->next;
        }

        if (start_edits != NULL) currSeql = start_edits->seql;
        if (sequences_list != NULL) prev_seq->next = curr_sequence;
        else sequences_list = curr_sequence;

        prev_seq = curr_sequence;
    }

    /****************************** Free lists ******************************/

    tokens *tmp_tok;
    while(first!=NULL) {
        free(first->value);
        tmp_tok=first;
        first = first->next;
        free(tmp_tok);
    }

    edits *temp_edits;
    while (edits_start!=NULL) {
        free(edits_start->seqPart);
        free(edits_start->qualPart);
        temp_edits = edits_start;
        edits_start = edits_start->next;
        free(temp_edits);
    }

    return sequences_list;
}


// Merge two matched reads into a single aligned read
int main(int argc, char **argv)
{
    if (argc != 2)
    {
        printf("Error: No input file provided\n");
        exit(-1);
    }

    bool paired = true;

    clock_t t;
    t = clock();

    FILE *fp = fopen(argv[1], "r");
    if (fp == NULL)
    {
        perror(argv[1]);
        exit(-1);
    }

    // Matrix
    int maxLength = 100;
    double **m = (double **)malloc(sizeof(double *) * maxLength);
    for (int i = 0; i < maxLength; i++) {
        m[i] = (double *)malloc(sizeof(double) * 4);
        for (int j = 0; j < 4; j++) {
            m[i][j] = 0;
        }
    }

    char *line = NULL;
    size_t line_buf_size = 0;
    ssize_t line_size;

    char *tok;
    const char delim[2] = "\t";

    // Hash table to keep track of repeated qname
    GHashTable* hash = g_hash_table_new_full(g_str_hash, g_str_equal, free, NULL);
    
    if (paired) {
        while((line_size = getline(&line, &line_buf_size, fp)) != -1)
        {
            if (line[0] == '@') continue;
            char *qname;

            tok = strtok(line, delim);

            // Assuming that each line in the file has the same number of tokens
            qname = malloc(strlen(tok) + 1);
            strcpy(qname, tok);
            
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);
            tok = strtok(NULL, delim);


            gpointer ret = g_hash_table_lookup(hash, qname);
            if (ret != NULL) 
                g_hash_table_insert(hash, qname, GINT_TO_POINTER(-1));
            else 
                g_hash_table_insert(hash, qname, GINT_TO_POINTER(1));
            
        }

        // Moves file pointer to the beginning
        fseek(fp, 0, SEEK_SET);
    }
    

    while((line_size = getline(&line, &line_buf_size, fp)) != -1)
    {
        // printf("Iter #: %d\n", iter++);
        if (line[0] == '@') continue;
        
        char *qname, *pos, *cigar, *seq, *qual;

        tok = strtok(line, delim);

        // Assuming that each line in the file has the same number of tokens
        qname = malloc(strlen(tok) + 1);
        strcpy(qname, tok);
        
        tok = strtok(NULL, delim);
        tok = strtok(NULL, delim);

        tok = strtok(NULL, delim);
        pos = malloc(strlen(tok) + 1);
        strcpy(pos, tok);

        tok = strtok(NULL, delim);

        tok = strtok(NULL, delim);
        cigar = malloc(strlen(tok) + 1);
        strcpy(cigar, tok);

        tok = strtok(NULL, delim);
        tok = strtok(NULL, delim);
        tok = strtok(NULL, delim);

        tok = strtok(NULL, delim);
        seq = malloc(strlen(tok) + 1);
        strcpy(seq, tok);

        tok = strtok(NULL, delim);
        qual = malloc(strlen(tok) + 1);
        strcpy(qual, tok);

        newSeq *mseqs = apply_cigar(cigar, pos, seq, qual);

        bool isRepeated = false;
        if (paired) {
            gpointer ret = g_hash_table_lookup(hash, qname);
            if (ret != NULL) {
                if (GPOINTER_TO_INT(ret) < 0) isRepeated = true;
            }
        }

        if (mseqs == NULL) continue;

        int posRanges[2];
        posRanges[0] = len_terminal_gap(mseqs->seq);
        posRanges[1] = strlen(mseqs->seq);

        char *aligned = mseqs->seq;
        int posRange = posRanges[0];
        int range = strlen(mseqs->seq);
        double **posVals=(double**)malloc(sizeof(double*) * 4);

        for (int j = 0; j < 4; j++) {
            posVals[j] = (double *)malloc((range - posRange) * sizeof(double));
        }
        
        double p;

        for (int k = posRange; k < range; k++) {
            char nt = aligned[k];
            char qc = (mseqs->qualseq)[k];

            if (nt == '-') {
                for (int j = 0; j < 4; j++) {
                    posVals[j][k - posRange] = 0;
                }
            }
            else if (nt == 'N') {
                for (int j = 0; j < 4; j++) {
                    posVals[j][k - posRange] = 0.25;
                }
            }
            else {
                p = pow(10, -1 * (((int)qc - 30)/ (double)10));

                for (int j = 0; j < 4; j++) {
                    posVals[j][k - posRange] = p/(double)3;
                }
                switch (nt) {
                    case 'A':
                        posVals[0][k - posRange] = 1-p;
                        break;
                    case 'C':
                        posVals[1][k - posRange] = 1-p;
                        break; 
                    case 'G':
                        posVals[2][k - posRange] = 1-p;
                        break; 
                    case 'T':
                        posVals[3][k - posRange] = 1-p;
                        break; 
                    default:
                        break;
                }
            }

            if (isRepeated) {
                for (int l = 0; l < 4; l++) {
                    posVals[l][k - posRange]/=2;
                }
            }   
        }     

        // int posRange = posRanges[0];
        // int range = posRanges[1];

        if (maxLength < range) {
            m = realloc(m, (sizeof(double *) * range));
            for (int t = maxLength; t < range; t++) {
                m[t] = (double *)malloc(sizeof(double) * 4);
                for (int y = 0; y < 4; y++) {
                    m[t][y] = 0;
                }
            }
            maxLength = range;
        }

        for (int k = posRange; k < range; k++) {
            for (int j = 0; j < 4; j++) {
                m[k][j] = posVals[j][k - posRange] + m[k][j];
            }
        }

        // Deallocate memory
        newSeq *temp_newSeq = mseqs;
        while (mseqs != NULL) {
            temp_newSeq = mseqs;
            free(temp_newSeq->seq);
            free(temp_newSeq->qualseq);
            mseqs = mseqs->next;
            free(temp_newSeq);
        }

        for (int f = 0; f < 4; f++) {
            // for (int k = 0; k < 4; k++) {
            //     free(posVals[f][k]);
            // }
            free(posVals[f]);
        }
        free(posVals);

        free(qname);
        free(pos);
        free(cigar);
        free(seq);
        free(qual);
    }
    fclose(fp);
    free(line);

    g_hash_table_destroy(hash);

    printf("Writing Matrix to CSV\n");
    write_matrix(m, maxLength, argv[1]);

    for (int i = 0; i < maxLength; i++) {
        free(m[i]);
    }
    free(m);

    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("Total Time: %f s\n", time_taken);

    return 0;
}