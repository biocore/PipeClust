#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <stdio.h>
#include <stdbool.h>

struct sequence_str{
    char* sequence __attribute__ ((aligned (16)));
    char* label;
    int seq_length;
    int label_length;
};
typedef struct sequence_str sequence;

/*
    Allocates the memory for a new sequence structure

    Returns a pointer to the new sequence structure
*/
sequence* new_sequence();

/*
    Reads the next sequence present in the file pointed by fd

    Returns a pointer to the read sequence structure or NULL if no
    more sequences are present on the file
*/
sequence* read_sequence(FILE *fd);

/*
    Reads the sequence number `idx` present in the file pointed by fd

    Returns a pointer to the read sequence structure or NULL if no more
    sequences are present on the file or the sequence has been already read
*/
sequence* read_sequence_by_idx(FILE *fd, int idx);

/*
    Writes the sequence pointed by seq in FASTA format to the
    file pointed by fd
*/
void write_sequence(sequence* seq, FILE* fd);

/*
    De-allocates all the memory used by sequence structure seq
*/
void free_sequence(sequence* seq);

/* Hackish - will revisit later */
void restore_counter();

#endif