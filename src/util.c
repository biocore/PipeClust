#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "util.h"

#define MSG_LEN 4096

/*
    Prints an error through stderr
    

    Inputs:
        level: integer that indicates if the error is fatal if != 0
        format: string that contains the text to be written to the stream.
            It can optionally contain embedded format specifiers that are
            replaced by the values specified in subsequent additional
            arguments and formatted as requested.
*/
void error_handler(int level, char* format, ...){
    int my_rank;
    int comm_sz;
    char message[MSG_LEN];
    va_list args;
    // Get the process rank
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    // Turn the parameters into a character string
    va_start(args, format);
    vsnprintf(message, MSG_LEN, format, args);
    va_end(args);
    // Check which type of error we are managing
    switch(level){
        case INFO_MSG:
            // Only process 0 shows this error
            if(my_rank == 0){
                // Print the warning message
                fprintf(stderr, "%s\n", message);
                fflush(stderr);
            }
            break;
        case WARN_ERROR:
            // Print the warning message
            fprintf(stderr, "[%d/%d] WARNING: %s\n", my_rank, comm_sz, message);
            fflush(stderr);
            break;
        case FATAL_ERROR:
            // Print the fatal error message
            fprintf(stderr, "[%d/%d] FATAL ERROR: %s\n", my_rank, comm_sz, message);
            fflush(stderr);
            // Exit with an error
            MPI_Abort(MPI_COMM_WORLD, -1);
        default:
            // Error level unknown - This is a fatal error
            fprintf(stderr, "[%d/%d] FATAL ERROR: error level unknown\n", my_rank, comm_sz);
            fflush(stderr);
            // Exit with an error
            MPI_Abort(MPI_COMM_WORLD, -1);
    }
}