#ifndef __PIPECLUST_H__
#define __PIPECLUST_H__

#include "derep_db.h"

/*
    De-replicates the list of files fasta_fps.

    Inputs:
        fasta_fps: list of fasta filepaths
        count: the number of fasta filepaths

    Returns a pointer to the de-replication database
*/
derep_db* serial_dereplication(char** fasta_fps, int count);

#endif