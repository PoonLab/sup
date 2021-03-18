#ifndef TYPES_H
#define TYPES_H

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


typedef struct timings
{
    double first_read;
    double cigar;
    double paired;
    double posrange;    
    double phred;
    double update;
} timings;

typedef struct tokens
{
	char *value;
	struct tokens *next;
} tokens;

typedef struct list_tokens
{
    struct tokens *toks;
    int pos;
    char *seq;
    char *qual;
    struct list_tokens *next;
} list_tokens;

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
    int lenSeq;
    char *seq;
    char *qual;
    struct insertions *next;
} insertions;

typedef struct newSeq {
    char *seq;
    char *qualseq;
    struct insertions *insertions;
    struct newSeq *next;
} newSeq;


#endif