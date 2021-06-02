#ifndef TYPES_H
#define TYPES_H
#include <glib.h>
#include <stdbool.h>

typedef struct in_file
{
    char *qname;
    char *flag;
    char *rname;
    char *pos;
    char *mapq;
    char *cigar;
    char *seq;
    char *qual;
} in_file;


typedef struct tokens
{
	char *value;
	struct tokens *next;
} tokens;

typedef struct edits {
    int len;
    char operand;
    int seql;
    int left;
    int lenSeq;
    char *seqPart;
    char *qualPart;
    struct edits *next;
} edits;

typedef struct insertions {
    int lenSeq; // Position of insertion
    char *seq; // Base
    char *qual;
    struct insertions *next;
} insertions;

typedef struct insertion_info {
    int lineNum;
    char *cigar;
    int seqPosition; // Position of insertion
    char *seq; // Base
    char *qual;
    bool isRepeated;
    struct insertion_info *next;
} insertion_info;

typedef struct newSeq {
    char *seq;
    char *qualseq;
    struct insertions *insertions;
    struct newSeq *next;
} newSeq;

typedef struct threadArg {
    int start_line;
    int end_line;
    char *filename;
    GHashTable *hash;
} threadArg;

#endif