#ifndef __DEREP_DB_H__
#define __DEREP_DB_H__

#include "sequence.h"
#include "uthash.h"
#include "utarray.h"

typedef struct seq_replicas_str {
    char* sequence __attribute__ ((aligned (16)));
    int count;
    UT_array *labels;
    UT_hash_handle hh;
} seq_replicas;

typedef struct derep_db_str {
    int count;
    int unique;
    seq_replicas* seqs;
} derep_db;

/*
    Creates a new de-replication database

    Returns:
        the new derep_db structure
*/
derep_db* create_derep_db(void);

/*
    Destroys the de-replication database db

    Inputs:
        db: pointer to the derep_db structure to destroy
*/
void destroy_derep_db(derep_db* db);

/*
    De-replicates the sequence seq against the de-replication database db.

    If the sequence is already present, the label is recorded. Otherwise,
    it is added as a new sequence in the de-replication database.

    Inputs:
        db: pointer to the derep_db structure
        seq: pointer to the sequence structure to be de-replicated
*/
void dereplicate_db(derep_db* db, sequence* seq);

/*
    Sorts the de-replication database by sequence abundance

    Inputs:
        db: pointer to the derep_db structure to sort
*/
void sort_db(derep_db* db);

/*
    Writes the de-replication database db in FASTA format to the `fasta` file
    and in an OTU map format in the `map` file

    Inputs:
        db: pointer to the derep_db struct
        fasta: string with the output fasta filename
        map: string with the output OTU map filename
*/
void write_output(derep_db* db, char* fasta, char* map);

/*
    Collects all the information about the derep_db spread across multiple
    processes in the process with rank 0

    Inputs:
        db: the local derep_db - will be modified in place
        my_rank: process rank
        comm_sz: the number of processes
*/
void gather_derep_db(derep_db* db, int my_rank, int comm_sz);

// void merge_derep_dbs(derep_db* db1, derep_db* db2);

#endif