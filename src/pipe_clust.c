#include <mpi.h>
#include <stdlib.h>
#include "pipe_clust.h"
#include "sequence.h"
#include "util.h"

/*
    De-replicates the fasta file fasta_fp against the de-replication
    database db

    Inputs:
        fasta_fp: fasta filepath
        db: pointer to the de-replication database
*/
void _serial_dereplication(char* fasta_fp, derep_db* db){
    // Open the FASTA file
    FILE* fd = fopen(fasta_fp, "r");
    // Check if we were able to open the file
    if(fd == NULL)
        error_handler(FATAL_ERROR, "Error opening file %s", fasta_fp);
    // Read the first sequence
    sequence* seq = read_sequence(fd);
    // Loop through all the file
    while(seq != NULL){
        // Compare if the sequence already exist on the DB
        dereplicate_db(db, seq);
        // Read next sequence
        seq = read_sequence(fd);
    }
    // Close the FASTA file
    fclose(fd);
}

/*
    De-replicates the list of files fasta_fps.

    Inputs:
        fasta_fps: list of fasta filepaths
        count: the number of fasta filepaths

    Returns a pointer to the de-replication database
*/
derep_db* serial_dereplication(char** fasta_fps, int num_files){
    int i;
    // Create the sequence DB
    derep_db* db = create_derep_db();
    // Loop through all the fasta files
    for(i = 0; i < num_files; i++){
        // Serially de-replicate current file against database
        _serial_dereplication(fasta_fps[i], db);
    }
    // Write a info message with the number of sequence read and the
    // number of unique sequences
    error_handler(INFO_MSG, "%d total sequences, %d unique sequences", db->count, db->unique);
    return db;
}
