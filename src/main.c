#include <mpi.h>
#include <stdio.h>
#include <getopt.h>
#include <string.h>
#include "pipe_clust.h"
#include "util.h"

static char* USAGE = "USAGE: mpiexec -n <NUM PROCS> PipeClust [cmd] [cmd options] FILE1 FILE2 ...\n"
                     "Run PipeClust --help for more information";

static char* HELP = "USAGE: mpiexec -n <NUM PROCS> PipeClust [cmd] [cmd options] FILE1 FILE2 ...\n\n"
                    "  cmd\n"
                    "    --help     Print this message\n"
                    "    --derep    Execute de-replication\n"
                    "\n"
                    "  --derep options:\n"
                    "    --fasta    Path to the output FASTA file\n"
                    "    --map      Path to the output OTU-map file\n";

int main(int argc, char** argv){
    // Start MPI
    MPI_Init(&argc, &argv);
    // Get my rank among all the processes
    int my_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    // Get the number of processes
    int comm_sz;
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);

    // Set up command line options parsing
    static int derep_flag = 0;
    static int sort_flag = 1;
    static int help_flag = 0;
    char* fasta = NULL;
    char* map = NULL;
    int option_index = 0;
    int c;
    int len;
    static struct option long_options[] = {
        {"derep", no_argument, &derep_flag, 1},
        {"suppress_sort", no_argument, &sort_flag, 0},
        {"help", no_argument, &help_flag, 1},
        {"fasta", required_argument, 0, 'f'},
        {"map", required_argument, 0, 'm'},
        {0, 0, 0, 0}
    };

    // Parse the command line options
    while((c = getopt_long(argc, argv, "f:m:", long_options, &option_index)) != -1){
        switch(c){
            case 0:
                // If this is setting a flag do nothing
                if(long_options[option_index].flag != 0)
                    break;
                // There is an error
                error_handler(INFO_MSG, "Error parsing the input options\n%s", USAGE);
                // Shut down MPI
                MPI_Finalize();
                return 0;
            case 'f':
                // We got the fasta file
                len = strlen(optarg);
                fasta = (char*) malloc(sizeof(char) * (len+1));
                memcpy(fasta, optarg, len);
                fasta[len] = '\0';
                break;
            case 'm':
                // We got the otu map file
                len = strlen(optarg);
                map = (char*) malloc(sizeof(char) * (len+1));
                memcpy(map, optarg, len);
                map[len] = '\0';
                break;
            case '?':
                break;
            default:
                error_handler(INFO_MSG, "Error parsing the input options\n%s", USAGE);
                // Shut down MPI
                MPI_Finalize();
                return 0;
        }
    }

    // Check the command line options

    // Check if the user requested the help string
    if(help_flag){
        error_handler(INFO_MSG, "%s", HELP);
        // Shut down MPI
        MPI_Finalize();
        return 0;
    }

    // Get the number of input files
    int num_files = argc - optind;
    if(num_files == 0){
        // No input files provided, throw the usage error
        error_handler(INFO_MSG, "Input files not provided!\n%s", USAGE);
        // Shut down MPI
        MPI_Finalize();
        return 0;
    }

    // Check de-replication options
    if (derep_flag){
        if(!fasta || !map){
            // No output files provided, throw the usage error
            error_handler(INFO_MSG, "If doing de-replication, both the output fasta file and the output otu_map should be defined. Fasta: %s, Otu Map: %s\n%s", fasta, map, USAGE);
            // Shut down MPI
            MPI_Finalize();
            return 0;
        }
    }
    else{
        error_handler(FATAL_ERROR, "Only de-replication is currently supported");
    }
    
    // Execute the commands
    if(derep_flag){
        // Executes de-replication
        derep_db* db = NULL;
        if(comm_sz == 1)
            db = serial_dereplication(&argv[optind], num_files);
        else
            db = parallel_dereplication(&argv[optind], num_files, my_rank, comm_sz);
        // At this point, only process with rank 0 has the
        // complete de-replication database
        if(my_rank == 0){
            // Write a info message with the number of sequence read and the
            // number of unique sequences
            error_handler(INFO_MSG, "%d total sequences, %d unique sequences", db->count, db->unique);
            if(sort_flag)
                // Sort the database by abundance
                sort_db(db);
            // Write the output files
            // TODO: probably remove when implementing further clustering steps
            write_output(db, fasta, map);
            // Destroy the sequence DB
            // TODO: probably remove when implementing further clustering steps
            destroy_derep_db(db);
        }
        else{
            // Destroy the sequence DB
            destroy_derep_db(db);
        }
    }
    
    // Shut down MPI
    MPI_Finalize();
    return 0;
}