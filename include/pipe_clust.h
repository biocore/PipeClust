#ifndef __PIPECLUST_H__
#define __PIPECLUST_H__

#include "derep_db.h"

/*
    Serially de-replicates the list of files fasta_fps.

    Inputs:
        fasta_fps: list of fasta filepaths
        count: the number of fasta filepaths

    Returns a pointer to the de-replication database
*/
derep_db* serial_dereplication(char** fasta_fps, int count);

/*
    De-replicates the list of files fasta_fps in parallel.

    Inputs:
        fasta_fps: list of fasta filepaths
        count: the number of fasta filepaths
        my_rank: the MPI rank of the current process
        comm_sz: the number of processes launched

    Returns a pointer to the de-replication database
*/
derep_db* parallel_dereplication(char** fasta_fps, int count, int my_rank, int comm_sz);

#endif