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
    return db;
}

/*
    De-replicates the fasta file fasta_fp against the de-replication database
    db, taking into account that other processes are accessing to the same file

    Inputs:
        fasta_fp: fasta filepath
        db: pointer to the de-replication database
        first: the index of the first sequence to access in the files
        n_partners: the number of processes accessing at the file
*/
void _parallel_dereplication(char* fasta_fp, derep_db* db, int first, int n_partners){
    // Holds the current sequence
    int current = first;
    // Open the FASTA file
    FILE *fd = fopen(fasta_fp, "r");
    // Check if we were able to open the file
    if(fd == NULL)
        error_handler(FATAL_ERROR, "Error opening file %s", fasta_fp);
    // Read the first sequence
    restore_counter();
    sequence* seq = read_sequence_by_idx(fd, current);
    // Loop through all the file
    while(seq != NULL){
        // Compare if the sequence already exist on the DB
        dereplicate_db(db, seq);
        // Update current sequence
        current += n_partners;
        // Read the next sequence
        seq = read_sequence_by_idx(fd, current);
    }
    // Close the FASTA file
    fclose(fd);
}

/*
    De-replicates the list of files fasta_fps in parallel.

    Inputs:
        fasta_fps: list of fasta filepaths
        num_files: the number of fasta filepaths
        my_rank: the MPI rank of the current process
        comm_sz: the number of processes launched

    Returns a pointer to the de-replication database
*/
derep_db* parallel_dereplication(char** fasta_fps, int num_files, int my_rank, int comm_sz){
    int i;
    // Will hold the current file processed
    int current;
    // Determine how many files I have to process by myself
    int num_my_files = num_files / comm_sz;
    // Determine how many files have left unassigned
    int remaining_files = num_files % comm_sz;

    // De-replication

    // Create the sequence DB
    derep_db* db = create_derep_db();
    // Loop through all the fasta file that are assigned to me
    // and I'm going to be the only process looking at it
    current = my_rank;
    for(i = 0; i < num_my_files; i++){
        // Since I'm the only one that will visit this file
        // I can de-replicate it serially
        _serial_dereplication(fasta_fps[current], db);
        // Update file index
        current += comm_sz;
    }

    // If the number of files cannot be evenly distributed among all 
    // the processes, there are some remaining files to be de-replicated
    if(remaining_files > 0){
        // Get the number of processes that are going to be accessing each file
        int n_partners = comm_sz / remaining_files;
        // To which file I should access? Note that, at this point,
        // remaining_files < comm_sz so each process only access to a one file
        current = (num_files - remaining_files) + (my_rank % remaining_files);
        // Check if there are unassigned processes that has not been taking
        // into account in n_partners, and they are going to access to a file
        int remaining_procs = comm_sz - (n_partners * remaining_files);
        // If my file is one of the files that the 'unassigned processes' are
        // going to access, I need to update my n_partners variable
        if(current < remaining_procs)
            ++n_partners;
        // Get which is the first sequence that I need to read
        int first_sequence = my_rank / remaining_files;
        // De-replicate shared file
        _parallel_dereplication(fasta_fps[current], db, first_sequence, n_partners);
    }

    // Gather results in a single process (rank=0)
    gather_derep_db(db, my_rank, comm_sz);
    // Return the de-replicated database
    return db;
}