#include "sequence.h"
#include "util.h"
#include <stdlib.h>
#include <string.h>

#define BUFFER_SIZE 2000 * sizeof(char)

// Global variable used to read by idx
int CURR_SEQ = 0;

/*
    Allocates the memory for a new sequence structure

    Returns a pointer to the new sequence structure
*/
sequence* new_sequence(){
    // Allocate memory for the sequence structure
    sequence* seq = (sequence*) malloc(sizeof(sequence));
    // Initialize all pointers to NULL
    seq->label = NULL;
    seq->sequence = NULL;
    // Initialize all lengths to 0
    seq->label_length = 0;
    seq->seq_length = 0;
    return seq;
}

/*
    Reads the next sequence present in the file pointed by fd

    Returns a pointer to the read sequence structure or NULL if no
    more sequences are present on the file
*/
sequence* read_sequence(FILE *fd){
    char* ret;
    // Allocate memory for reading buffer
    char* buffer = (char*) malloc(BUFFER_SIZE);
    
    // Read sequence label
    ret = fgets(buffer, BUFFER_SIZE, fd);
    if(ret == NULL){
        // Free allocated memory
        free(buffer);
        // Check if there has been an error while reading the FASTA
        // file or we simply have reach the end of the file
        if(feof(fd) == 0)
            // An error occurred, terminate execution
            error_handler(FATAL_ERROR, "Error reading the FASTA file");
        return NULL;
    }
    
    // Allocate memory for sequence
    sequence* seq = (sequence*) malloc(sizeof(sequence));
    // Allocate memory for label scanning
    char* label_buffer = (char*) malloc(BUFFER_SIZE);
    // Scan the label
    if(sscanf(buffer, ">%s%*s", label_buffer) != 1){
        // Free allocated memory
        free(buffer);
        free(label_buffer);
        free(seq);
        // An error occurred, terminate execution
        error_handler(FATAL_ERROR, "Error parsing sequence %d label", CURR_SEQ);
    }
    // Set the label info in the sequence struct
    seq->label_length = strlen(label_buffer);
    seq->label = (char*) malloc(sizeof(char) * seq->label_length + 1);
    memcpy(seq->label, label_buffer, seq->label_length);
    seq->label[seq->label_length] = '\0';
    // The label buffer is no longer needed - free its memory
    free(label_buffer);

    // Read sequence
    ret = fgets(buffer, BUFFER_SIZE, fd);
    if(ret == NULL){
        // Free allocated memory
        free(buffer);
        free(seq->label);
        free(seq);
        // Since we have already read the label, we can safely throw an error
        // because either it has been an error reading the FASTA file or 
        // the FASTA file is not correct -> it ends with a label...
        error_handler(FATAL_ERROR, "Error reading sequence %d from the FASTA file", CURR_SEQ);
    }
    // Set the sequence info in the sequence struct
    // fgets puts the \n character also in the buffer,
    // so the real sequence length is strlen(buffer) - 1
    seq->seq_length = strlen(buffer) - 1;
    // seq->sequence = (char*) malloc(sizeof(char) * (seq->seq_length+1));
    posix_memalign((void **) &seq->sequence, 16, sizeof(char) * (seq->seq_length+1));
    memcpy(seq->sequence, buffer, seq->seq_length);
    seq->sequence[seq->seq_length] = '\0';

    // Free reading buffer memory
    free(buffer);
    // Update CURR_SEQ, as we have read a sequence
    ++CURR_SEQ;
    // Return the sequence
    return seq;
}

/*
    Reads the sequence number `idx` present in the file pointed by fd

    Returns a pointer to the read sequence structure or NULL if no more
    sequences are present on the file or the sequence has been already read
*/
sequence* read_sequence_by_idx(FILE *fd, int idx){
    // Check that we don't have already read the sequence
    if(CURR_SEQ > idx)
        return NULL;
    // Skip sequences until the one we have to read
    sequence* seq;
    do{
        seq = read_sequence(fd);
    } while(CURR_SEQ <= idx && seq != NULL);
    // Return the sequence to read
    return seq;
}

/*
    Writes the sequence pointed by seq in FASTA format to the
    file pointed by fd
*/
void write_sequence(sequence* seq, FILE* fd){
    fprintf(fd, "%s\n%s\n", seq->label, seq->sequence);
}

/*
    De-allocates all the memory used by sequence structure seq
*/
void free_sequence(sequence* seq){
    free(seq->label);
    free(seq->sequence);
    free(seq);
}

/* Hackish - will revisit later */
void restore_counter(){
    CURR_SEQ = 0;
}