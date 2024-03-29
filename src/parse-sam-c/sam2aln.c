#include "io.h"
#include <time.h>
#include <regex.h>
#include <pthread.h>
#include <unistd.h>

static const bool paired = true;
int maxLength;
pthread_mutex_t lock;
double **m;
insertion_info *insertions_list;
insertion_info *prev_insertion;

void print_usage(char *pname) {
    printf("Usage: %s [OPTIONS]\n", pname);
    printf("Options:\n");
    printf("  -h            Display the help message\n");
    printf("  -f FILE       Input SAM filename\n");
    printf("  -t NUM        Number of threads\n");
}

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
    if (strcmp(cigar, "*") == 0) return NULL;

    tokens *first = NULL;
    tokens *prev_tok = first;

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

    // POS is  1-indexed, so shift to the left
    int curr_pos = atoi(pos) - 1;

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
                memset(current_edit->seqPart, 'x', current_edit->len);
                (current_edit->seqPart)[current_edit->len] = '\0';

                current_edit->qualPart = malloc(current_edit->len + 1);
                memset(current_edit->qualPart, '!', current_edit->len);
                (current_edit->qualPart)[current_edit->len] = '\0';
            }
            else {
                current_edit->seqPart = NULL;
                current_edit->qualPart = NULL;
            }

            // Appends sequence together
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
                else if (current_edit->left - 1 < strlen(qual)) {
                    current_edit->qualPart = malloc(current_edit->len + 1);
                    strncpy(current_edit->qualPart, qual + (current_edit->left - 1), current_edit->len);
                    (current_edit->qualPart)[current_edit->len] = '\0';
                }
            }

            // Update lenSeq
            if (current_edit->operand == 'I') {
                // Store the position where the insertion starts
                current_edit->lenSeq = curr_pos;
                curr_pos+=current_edit->len;
            }
            else {
                curr_pos+=current_edit->len;
                current_edit->lenSeq = curr_pos;
            }

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

        insertions *head_insertions = NULL;
        insertions *prev_insert = head_insertions;

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

		            if (strcmp(qual, "*") !=0) {
                      int qual_len = strlen(curr_sequence->qualseq) + strlen(start_edits->qualPart) + 1;
                      curr_sequence->qualseq = realloc(curr_sequence->qualseq, qual_len);
                      strncat(curr_sequence->qualseq, start_edits->qualPart, strlen(start_edits->qualPart) + 1);
                    }
                    else {
                        curr_sequence->qualseq = realloc(curr_sequence->qualseq, curr_len);
                        memset(curr_sequence->qualseq, '!', curr_len);
                    }
                }
            }

            if (start_edits->operand == 'I') {
                insertions *curr_insert = malloc(sizeof(insertions));
                curr_insert->lenSeq = start_edits->lenSeq + 1 - (atoi(pos) - 1);
                curr_insert->seq = malloc(strlen(start_edits->seqPart) + 1);
                curr_insert->qual = malloc(strlen(start_edits->qualPart) + 1);
                curr_insert->next = NULL;

                strncpy(curr_insert->seq, start_edits->seqPart, strlen(start_edits->seqPart));
                (curr_insert->seq)[strlen(start_edits->seqPart)] = '\0';
                
                strncpy(curr_insert->qual, start_edits->qualPart, strlen(start_edits->qualPart));
                (curr_insert->qual)[strlen(start_edits->qualPart)] = '\0';

                if (head_insertions != NULL) prev_insert->next = curr_insert;
                else head_insertions = curr_insert;
                prev_insert = curr_insert;
            }

            start_edits = start_edits->next;
        }

        curr_sequence->insertions = head_insertions;


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

void *processLines(void *args)
{
    threadArg *arg = (threadArg *)args;
    char *filename = arg->filename;
    GHashTable *hash = arg->hash;
    int start = arg->start_line;
    int end = arg->end_line;
    int numLines = 0;
    FILE *fp = fopen(filename, "r");
    if (fp == NULL)
    {
        perror(filename);
        exit(-1);
    }
    char *line = NULL;
    size_t line_buf_size = 0;
    ssize_t line_size;

    char *tok;
    char *context = NULL;
    const char delim[2] = "\t";
    while((line_size = getline(&line, &line_buf_size, fp)) != -1)
    {
        numLines++;
        if (numLines > start && numLines <= end) {
            if (line[0] == '@') continue;
            
            char *qname, *pos, *cigar, *seq, *qual;

            tok = strtok_r(line, delim, &context);

            // Assuming that each line in the file has the same number of tokens
            qname = malloc(strlen(tok) + 1);
            strcpy(qname, tok);
            
            tok = strtok_r(NULL, delim, &context);
            tok = strtok_r(NULL, delim, &context);

            tok = strtok_r(NULL, delim, &context);
            pos = malloc(strlen(tok) + 1);
            strcpy(pos, tok);

            tok = strtok_r(NULL, delim, &context);

            tok = strtok_r(NULL, delim, &context);
            cigar = malloc(strlen(tok) + 1);
            strcpy(cigar, tok);

            tok = strtok_r(NULL, delim, &context);
            tok = strtok_r(NULL, delim, &context);
            tok = strtok_r(NULL, delim, &context);

            tok = strtok_r(NULL, delim, &context);
            seq = malloc(strlen(tok) + 1);
            strcpy(seq, tok);

            tok = strtok_r(NULL, delim, &context);
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

            if (mseqs->insertions != NULL) {
                // Add insertions to list
                pthread_mutex_lock(&lock);
                insertions *mseq_insertions = mseqs->insertions;
                while (mseq_insertions != NULL) {
                    insertion_info *insert = malloc(sizeof(insertion_info));
                    insert->lineNum = numLines;
                    insert->cigar = malloc(strlen(cigar) + 1);
                    strcpy(insert->cigar, cigar);
                    insert->seqPosition = mseq_insertions->lenSeq;
                    insert->isRepeated = isRepeated;
                    insert->seq = malloc(strlen(mseq_insertions->seq) + 1);
                    insert->qual = malloc(strlen(mseq_insertions->qual) + 1);
                    strcpy(insert->seq, mseq_insertions->seq);
                    strcpy(insert->qual, mseq_insertions->qual);
                    insert->next = NULL;
                    mseq_insertions = mseq_insertions->next;

                    if (insertions_list != NULL) prev_insertion->next = insert;
                    else insertions_list = insert;

                    prev_insertion = insert;
                }
                pthread_mutex_unlock(&lock);
            }

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
                else if (nt == 'x') {
                    for (int j = 0; j < 4; j++) {
                        posVals[j][k - posRange] = 1;
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

            pthread_mutex_lock(&lock);
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
            pthread_mutex_unlock(&lock);

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
                insertions *temp_insert;
                while (temp_newSeq->insertions != NULL) {
                    temp_insert = temp_newSeq->insertions;
                    free(temp_insert->seq);
                    free(temp_insert->qual);
                    temp_newSeq->insertions = temp_newSeq->insertions->next;
                    free(temp_insert);
                }
                mseqs = mseqs->next;
                free(temp_newSeq);
            }

            for (int f = 0; f < 4; f++) {
                free(posVals[f]);
            }
            free(posVals);

            free(qname);
            free(pos);
            free(cigar);
            free(seq);
            free(qual);
        }
    }
    fclose(fp);
    free(line);

    free(arg);
    return 0;
}



// Merge two matched reads into a single aligned read
int main(int argc, char **argv)
{
    int c;
    char *input_file = NULL;
    char *nthreads = NULL;

    while ((c = getopt(argc, argv, "hf:t:")) != -1) {
        switch (c) {
            case 'h':
                print_usage(argv[0]);
                return 0;
            case 'f':
                input_file = optarg;
                break;
            case 't':
                nthreads = optarg;
                break;
            case '?':
                fprintf(stderr, "Invalid option: %c\n", optopt);
                break;
        }
    }

    if (input_file == NULL) {
        fprintf(stderr, "Error: input file is required\n");
        print_usage(argv[0]);
        return 1;
    }

    if (nthreads == NULL) {
        printf("Note: 2 threads will be used by default\n");
        nthreads = "2";
    }

    if (pthread_mutex_init(&lock, NULL) != 0)
    {
        printf("\n Mutex Initialization Failed\n");
        exit(1);
    }

    clock_t t;
    t = clock();

    FILE *fp = fopen(input_file, "r");
    if (fp == NULL)
    {
        perror(input_file);
        exit(-1);
    }

    // Matrix
    maxLength = 100;
    m = (double **)malloc(sizeof(double *) * maxLength);
    for (int i = 0; i < maxLength; i++) {
        m[i] = (double *)malloc(sizeof(double) * 4);
        for (int j = 0; j < 4; j++) {
            m[i][j] = 0;
        }
    }

    char *tok, *line = NULL;
    size_t line_buf_size = 0;
    ssize_t line_size;
    const char delim[2] = "\t";
    int numLines = 0;

    // Hash table to keep track of repeated qname
    GHashTable* hash = g_hash_table_new_full(g_str_hash, g_str_equal, free, NULL);
    
    if (paired) {
        while((line_size = getline(&line, &line_buf_size, fp)) != -1)
        {
            numLines++;
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
    }
    else {
        numLines = get_number_rows(input_file);
    }
    fclose(fp);
    free(line);
    
    int num_read_lines = numLines / atoi(nthreads);
    int num_threads = atoi(nthreads);
    int start_line = 0, end_line = num_read_lines;

    if (numLines % num_threads > 0) {
        num_threads++;
    }  

    pthread_t tid[num_threads];


    for (int k = 0; k < num_threads; k++) {
        if (end_line > numLines) end_line = numLines;

        threadArg *arg = malloc(sizeof(threadArg));
        arg->start_line = start_line;
        arg->end_line = end_line;
        arg->filename = input_file;
        arg->hash = hash;
        pthread_create(&tid[k], NULL, processLines, (void *)arg);


        start_line += num_read_lines;
        end_line += num_read_lines;
    }
  
    for (int k = 0; k < num_threads; k++) {
        pthread_join(tid[k], NULL);
    }

    g_hash_table_destroy(hash);
    pthread_mutex_destroy(&lock);

    printf("Writing Matrix to CSV\n");
    write_matrix(m, maxLength, input_file);

    if (insertions_list!=NULL) {
        printf("Writing Insertions to CSV\n");
        write_insertions(insertions_list, input_file);
    }

    for (int i = 0; i < maxLength; i++) {
        free(m[i]);
    }
    free(m);

    // Free List of insertions
    insertion_info *tempInsertions;
    while (insertions_list!=NULL) {
        free(insertions_list->cigar);
        free(insertions_list->seq);
        free(insertions_list->qual);
        tempInsertions = insertions_list;
        insertions_list = insertions_list->next;
        free(tempInsertions);
    }

    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    printf("Total Time: %f s\n", time_taken);

    return 0;
}