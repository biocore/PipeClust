#include <stdlib.h>
#include <string.h>
#include "derep_db.h"

/************************************
 *   Replica structure functions    *
************************************/

/*
    Creates a new seq_replicas structure with the 'seq'

    Inputs:
        seq: pointer to the sequence structure

    Returns a pointer to the new seq_replicas structure
*/
seq_replicas* create_seq_replica(sequence* seq){
    // Allocate memory for the new replica structure
    seq_replicas* r = (seq_replicas*) malloc(sizeof(seq_replicas));
    // Initialize replica structure with the new sequence
    r->count = 1;
    // Allocate memory for the sequence string
    posix_memalign((void **) &r->sequence, 16, sizeof(char) * (seq->seq_length+1));
    // Copy the sequence string
    memcpy(r->sequence, seq->sequence, seq->seq_length+1);
    // Initialize labels array
    utarray_new(r->labels, &ut_str_icd);
    // Insert the sequence's label to the labels array
    utarray_push_back(r->labels, &seq->label);
    return r;
}

/*
    Adds the sequence seq as a replica of r

    Inputs:
        r: pointer to the seq_replicas structure
        seq: pointer to the sequence structure
*/
void add_replica(seq_replicas* r, sequence* seq){
    // Add the sequence's label to the labels array
    utarray_push_back(r->labels, &seq->label);
    // Update the counter
    ++r->count;
}

/*
    Destroy the replica structure r

    Inputs:
        r: pointer to the seq_replicas structure
*/
void destroy_seq_replica(seq_replicas* r){
    // Free up the labels array
    utarray_free(r->labels);
    // Free the sequence memory
    free(r->sequence);
    // Free the whole structure memory
    free(r);
}

/*
    Auxiliary function that compares two seq_replica structures by abundance
    For descendant sorting purposes

    Inputs:
        a, b: pointers to the seq_replicas structures

    Returns:
        -# if a > b
        0 if a == b
        +# if a < b
*/
int compare_abundances(seq_replicas* a, seq_replicas* b){
    return (b->count - a->count);
}

/************************************
* De-replication database functions *
************************************/

/*
    Creates a new de-replication database

    Returns:
        the new derep_db structure
*/
derep_db* create_derep_db(void){
    // Allocate memory for new structure
    derep_db* db = (derep_db*) malloc(sizeof(derep_db));
    // There is no sequence still, so count and unique are 0
    db->count = 0;
    db->unique = 0;
    // Initialize the hash table to NULL
    db->seqs = NULL;
    // Return the new database
    return db;
}

/*
    Destroys the de-replication database db

    Inputs:
        db: pointer to the derep_db structure to destroy
*/
void destroy_derep_db(derep_db* db){
    // Free up the memory of the seq_replicas structure
    seq_replicas* current;
    seq_replicas* tmp;
    HASH_ITER(hh, db->seqs, current, tmp){
        HASH_DEL(db->seqs, current);
        destroy_seq_replica(current);
    }
    // Free all memory
    free(db);
}

/*
    De-replicates the sequence seq. If the sequence is already present,
    the label is recorded and the sequence struct freed up. Otherwise,
    it is added as a new sequence in the de-replication database and 
    the structure is kept.

    Inputs:
        db: pointer to the derep_db structure
        seq: pointer to the sequence structure to be de-replicated
*/
void dereplicate_db(derep_db* db, sequence* seq){
    // Check if the sequence already exists on the DB
    seq_replicas* r;
    HASH_FIND_STR(db->seqs, seq->sequence, r);
    if(r){
        // The sequence was already present on the DB
        add_replica(r, seq);
    }
    else{
        // The sequence didn't exist, add as new sequence
        r = create_seq_replica(seq);
        HASH_ADD_KEYPTR(hh, db->seqs, r->sequence, strlen(r->sequence), r);
        // Update unique counter
        ++db->unique;
    }
    // Update counter
    ++db->count;
}

/*
    Sorts the de-replication database by sequence abundance

    Inputs:
        db: pointer to the derep_db structure to sort
*/
void sort_db(derep_db* db){
    HASH_SORT(db->seqs, compare_abundances);
}

/*
    Writes the de-replication database db in FASTA format to the `fasta` file
    and in an OTU map format in the `map` file

    Inputs:
        db: pointer to the derep_db struct
        fasta: string with the output fasta filename
        map: string with the output OTU map filename
*/
void write_output(derep_db* db, char* fasta, char* map){
    int i;
    // Open fasta and OTU map files
    FILE* fasta_fd = fopen(fasta, "w");
    FILE* map_fd = fopen(map, "w");

    // Loop through all the sequences
    seq_replicas* current;
    seq_replicas* tmp;
    char **l;
    i = 0;
    HASH_ITER(hh, db->seqs, current, tmp){
        // Write sequence into the fasta file
        fprintf(fasta_fd, ">Seq_%d count=%d\n%s\n", i, current->count,current->sequence);
        // Write OTU id in the OTU map
        fprintf(map_fd, "Seq_%d", i);
        // Write all the labels
        l = NULL;
        while((l=(char**)utarray_next(current->labels, l))){
            fprintf(map_fd, "\t%s", *l);
        }
        // Current OTU done - write new line character in the OTU map
        fprintf(map_fd, "\n");
        // and increment OTU counter
        i++;
    }
    // Close files
    fclose(fasta_fd);
    fclose(map_fd);
}
