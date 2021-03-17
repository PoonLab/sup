#include "io.h"
#include <time.h>
#include <regex.h>
#include <stdbool.h>
#include <math.h>

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

newSeq * apply_cigar(in_file *table, int n) {
    // TODO: Check if CIGAR is valid

    int num_valid = 0;
    
    /************************** Tokenize CIGAR ************************/
    list_tokens *head = NULL;
    list_tokens *prev = head;

    for (int i = 0; i < n; i++) {
        if (strcmp(table[i].cigar, "*") == 0) continue;

        tokens *first = NULL;
        tokens *prev_tok = first;

        regex_t pexp;
        regmatch_t whole_match;

        char *sz = table[i].cigar;

        regcomp(&pexp, "[0-9]+[MIDNSHPX=]", REG_EXTENDED);

        int eflags = 0;
        int match = 0;
        // int leng;
        size_t offset = 0;
        size_t length = strlen(sz);

        while (regexec(&pexp, sz + offset, 1, &whole_match, eflags) == 0) {
            eflags = REG_NOTBOL;
            match = 1;

            tokens *current_tok = malloc(sizeof(tokens));
            current_tok->next = NULL;

            size_t leng = (offset + whole_match.rm_eo) - (offset + whole_match.rm_so);
            current_tok->value = malloc(leng + 1);
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

        list_tokens *current_tok_list = malloc(sizeof(list_tokens));
        current_tok_list->toks = first;
        current_tok_list->next = NULL;
        char *pos_table = table[i].pos;
        current_tok_list->pos = atoi(pos_table);
       
        current_tok_list->seq = malloc(strlen(table[i].seq) + 1);
        strncpy(current_tok_list->seq, table[i].seq, strlen(table[i].seq));
        (current_tok_list->seq)[strlen(table[i].seq)] = '\0';

        current_tok_list->qual = malloc(strlen(table[i].qual) + 1);
        strncpy(current_tok_list->qual, table[i].qual, strlen(table[i].qual));
        (current_tok_list->qual)[strlen(table[i].qual)] = '\0';

        if (head != NULL) prev->next = current_tok_list;
        else head=current_tok_list;

        prev = current_tok_list;

        regfree(&pexp);
        num_valid++;
    }

    /******************************** Edits *********************************/

    edits *edits_start = NULL;
    edits *prev_edits = edits_start;

    list_tokens *temp_list = head;
    int seql = 0;
    int left;
    int flag;
    int curr_pos;
    while(temp_list!=NULL) {
        tokens *temp_tok = temp_list->toks;
        left = 1;
        flag = 1;
        seql++;
        curr_pos = temp_list->pos;
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
                current_edit->len = atoi(num);
                current_edit->operand = curr_tok[strlen(curr_tok) - 1];
                current_edit->seql = seql;
                current_edit->left = left;

                // Free num as it is no longer needed
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
                        strncpy(sequence, temp_list->seq + (current_edit->left - 1), current_edit->len);
                        (sequence)[current_edit->len] = '\0';
                        int len_target = strlen(current_edit->seqPart) + current_edit->len + 1;
                        current_edit->seqPart = realloc(current_edit->seqPart,len_target);
                        strncat(current_edit->seqPart, sequence, current_edit->len + 1);
                    }
                    else {
                        current_edit->seqPart = malloc(current_edit->len + 1);
                        strncpy(current_edit->seqPart, temp_list->seq + (current_edit->left - 1), current_edit->len);
                        (current_edit->seqPart)[current_edit->len] = '\0';
                    }

                    if (current_edit->qualPart != NULL) {
                        char *qual_part = malloc(current_edit->len + 1);
                        strncpy(qual_part, temp_list->qual + (current_edit->left - 1), current_edit->len);
                        (qual_part)[current_edit->len] = '\0';
                        int len_target = strlen(current_edit->qualPart) + current_edit->len + 1;
                        current_edit->qualPart = realloc(current_edit->qualPart ,len_target);
                        strncat(current_edit->qualPart, qual_part, current_edit->len + 1);
                    }
                    else {
                        current_edit->qualPart = malloc(current_edit->len + 1);
                        strncpy(current_edit->qualPart, temp_list->qual + (current_edit->left - 1), current_edit->len);
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
        }
        temp_list = temp_list->next;
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

                    int qual_len = strlen(curr_sequence->qualseq) + strlen(start_edits->qualPart) + 1;
                    curr_sequence->qualseq = realloc(curr_sequence->qualseq, qual_len);
                    strncat(curr_sequence->qualseq, start_edits->qualPart, strlen(start_edits->qualPart) + 1);
                }
            }

            if (start_edits->operand == 'I') {
                insertions *curr_insert = malloc(sizeof(insertions));
                curr_insert->lenSeq = start_edits->lenSeq;
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

    list_tokens *tmp_list;
    while(head!=NULL) {
        tokens *tmp_tok;
        while(head->toks!=NULL) {
            free(head->toks->value);
            tmp_tok=head->toks;
            head->toks = head->toks->next;
            free(tmp_tok);
        }
        tmp_list = head;
        head = head->next;
        free(tmp_list);
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
        exit(1);
    }

    time_t t0, t1;
    float total_timing = 0;
    time(&t0);

    printf("Reading File\n");

    // Creating the Data Frame
    int n = get_number_rows(argv[1]);
    in_file *table = malloc(n * sizeof(in_file));
    read_input_file(argv[1], &table[0], n);

    time(&t1);
    total_timing += difftime(t1, t0);
    printf("File Read\n");
    printf("First Read: %f s\n", difftime(t1, t0));

    // Apply Cigar values
    time(&t0);
    newSeq *mseqs = apply_cigar(&table[0], n);
    time(&t1);
    total_timing += difftime(t1, t0);
    printf("Cigar Strings Processed\n");
    printf("Cigar: %f s\n", difftime(t1, t0));


    printf("Preparing Position Data\n");

    // TODO:  Merge those paired neighbours mseq values

    // Calculates the longest an mseq value could be
    int maxLen = 0;
    int numSeqs = 0;
    int currlength;
    newSeq *head_mseqs = mseqs;
    while (head_mseqs!=NULL) {
        currlength = strlen(head_mseqs->seq);
        if (currlength > maxLen) maxLen = currlength;
        head_mseqs = head_mseqs->next;
        numSeqs++;
    }

    // printf("%d\n", maxLen);

    // 0 -> A, 1 -> C, 2 -> G, 3 -> T
    double m[maxLen][4];

    // Initialize the m matrix
    for (int r = 0; r < maxLen; r++) {
        for (int c = 0; c < 4; c++) {
            m[r][c] = 0;
        }
    }

    // For Targetting, 0->number of '-' at the start/end of sequence, 1->length of the sequence
    int posRanges[numSeqs][2];
    int i = 0;
    newSeq *seqTem = mseqs;
    while(seqTem != NULL) {
        posRanges[i][0] = len_terminal_gap(seqTem->seq);
        posRanges[i][1] = strlen(seqTem->seq);
        seqTem = seqTem->next;
        i++;
    }

    //TODO: Only changes to merged value when paired

    // double posVals[numSeqs][4];
    double ***posVals = (double***)malloc(sizeof(double**) * numSeqs);
    seqTem = mseqs;
    i = 0;
    double p;
    while (seqTem != NULL) {
        char *aligned = seqTem->seq;
        int posRange = posRanges[i][0];
        int range = strlen(seqTem->seq);
        posVals[i] = (double **)malloc(sizeof(double *) * 4);
        
        for (int j = 0; j < 4; j++) {
            posVals[i][j] = (double*) malloc((range - posRange) * sizeof(double));
        }

        for (int k = posRange; k < range; k++) {
            char nt = aligned[k];
            char qc = (seqTem->qualseq)[k];

            if (nt == '-') {
                for (int j = 0; j < 4; j++) {
                    posVals[i][j][k - posRange] = 0;
                }
            }
            else if (nt == 'N') {
                for (int j = 0; j < 4; j++) {
                    posVals[i][j][k - posRange] = 0.25;
                }
            }
            else {
                p = pow(10, -1 * (((int)qc - 30)/ (double)10));

                for (int j = 0; j < 4; j++) {
                    posVals[i][j][k - posRange] = p/(double)3;
                }
                switch (nt) {
                    case 'A':
                        posVals[i][0][k - posRange] = 1-p;
                        break;
                    case 'C':
                        posVals[i][1][k - posRange] = 1-p;
                        break; 
                    case 'G':
                        posVals[i][2][k - posRange] = 1-p;
                        break; 
                    case 'T':
                        posVals[i][3][k - posRange] = 1-p;
                        break; 
                    default:
                        break;
                }
            }
        }
        seqTem = seqTem->next;
        i++;
    }

    printf("Final Step: Applying position Values\n");

    // Increase targeted areas
    for (int i = 0; i < numSeqs; i++) {
        int posRange = posRanges[i][0];
        int range = posRanges[i][1];
        for (int k = posRange; k < range; k++) {
            for (int j = 0; j < 4; j++) {
                m[k][j] = posVals[i][j][k - posRange] + m[k][j];
            }
        }
    }

    printf("Writing Matrix to CSV\n");
    write_matrix(m, maxLen);

    // Free input struct
    for (int i=0; i < n; i++) 
    {
        free(table[i].qname);
        free(table[i].flag);
        free(table[i].rname);
        free(table[i].pos);
        free(table[i].mapq);
        free(table[i].cigar);
        free(table[i].seq);
        free(table[i].qual);
    }
    free(table);

    newSeq *temp_newSeq = mseqs;
    while (mseqs != NULL) {
        insertions *temp_insert;
        while (mseqs->insertions != NULL) {
            temp_insert = mseqs->insertions;
            free(temp_insert->seq);
            free(temp_insert->qual);
            mseqs->insertions = mseqs->insertions->next;
            free(temp_insert);
        }
        temp_newSeq = mseqs;
        free(temp_newSeq->seq);
        free(temp_newSeq->qualseq);
        mseqs = mseqs->next;
        free(temp_newSeq);
    }

    for (int f = 0; f < numSeqs; f++) {
        for (int k = 0; k < 4; k++) {
            free(posVals[f][k]);
        }
        free(posVals[f]);
    }
    free(posVals);

    return 0;
}