#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
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
    Creates a new seq_replicas structure with 'sequence' but with no labels on it

    Inputs:
        seq: char array of the sequence
        seq_length: the length of the sequence

    Returns a pointer to the new seq_replicas structure
*/
seq_replicas* create_empty_seq_replica(char* sequence, int seq_length){
    // Allocate memory for the new replica structure
    seq_replicas* r = (seq_replicas*) malloc(sizeof(seq_replicas));
    // Initialize counter to 0 because no labels are on it
    r->count = 0;
    // Allocate memory for the sequence string
    posix_memalign((void **) &r->sequence, 16, sizeof(char) * (seq_length+1));
    // Copy the sequence string
    memcpy(r->sequence, sequence, seq_length+1);
    // Initialize labels array
    utarray_new(r->labels, &ut_str_icd);
    return r;
}

/*
    Adds the sequence seq as a replica of r

    Inputs:
        r: pointer to the seq_replicas structure
        seq: pointer to the sequence structure
*/
void add_replica(seq_replicas* r, char* label){
    // Add the sequence's label to the labels array
    utarray_push_back(r->labels, &label);
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

/********************************************
* De-replication database private functions *
********************************************/

/* Communication functions */

/*
    Packs the de-replication db in a char array and returns it

    Inputs:
        db: de-replication database to serialize
        msg_size: output parameter - holds the length of the char array with
            the de-replicated db packed

    Returns a pointer to the char array with the packed message
*/
char* pack_derep_db(derep_db* db, int* msg_size){
    int position;
    // Re-allcoating memory is expensive, lets try to guess a good size
    // for the char array that will hold the message
    // 2 ints for seq count and unique count
    int size = 2 * sizeof(int);
    // For each sequence, allocate 500 chars + 2 ints: seq length and num labels
    size += (db->unique * (sizeof(char) * 500 + 2 * sizeof(int)));
    // For each label, allocate 100 chars + 1 int: label length
    size += (db->count * (sizeof(char) * 100 + sizeof(int)));

    // Allocate memory for the message
    char* buffer = (char*) malloc(size);

    // Start packing the message by adding the count and unique counters
    position = 0;
    MPI_Pack(&db->count, 1, MPI_INT, buffer, size, &position, MPI_COMM_WORLD);
    MPI_Pack(&db->unique, 1, MPI_INT, buffer, size, &position, MPI_COMM_WORLD);
    
    // Loop through all the unique sequences present in the db
    seq_replicas* current;
    seq_replicas* tmp;
    char** label;
    int length;
    HASH_ITER(hh, db->seqs, current, tmp){
        // Pack the length of the sequence
        length = strlen(current->sequence);
        MPI_Pack(&length, 1, MPI_INT, buffer, size, &position, MPI_COMM_WORLD);
        // Pack the sequence
        MPI_Pack(current->sequence, length, MPI_CHAR, buffer, size, &position, MPI_COMM_WORLD);
        // Pack the number of labels
        MPI_Pack(&current->count, 1, MPI_INT, buffer, size, &position, MPI_COMM_WORLD);
        // Loop through all the labels
        label = NULL;
        while((label=(char**)utarray_next(current->labels, label))){
            // Pack the label length
            length = strlen(*label);
            MPI_Pack(&length, 1, MPI_INT, buffer, size, &position, MPI_COMM_WORLD);
            // Pack the label
            MPI_Pack(*label, length, MPI_CHAR, buffer, size, &position, MPI_COMM_WORLD);
        }
    }
    // The actual message size is hold by position
    *msg_size = position;
    // Return the char array with the message
    return buffer;
}

/*
    Sends the de-replication database to process dest

    Inputs:
        db: de-replication database to send
        my_rank: this process rank
        dest: the rank of the destination process
*/
void _send_derep_db(derep_db* db, int my_rank, int dest){
    // We need to send two messages
    // The first one contains the size of second message
    // The second one contains the serialized de-replication db
    int size;
    // Pack the db into a message so we know the size of it
    char* msg = pack_derep_db(db, &size);
    // Send the message with the size info
    MPI_Send(&size, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
    // We can send now the database as the receiver has allocated enough
    // memory to receive it
    MPI_Send(msg, size, MPI_PACKED, dest, 0, MPI_COMM_WORLD);
    // We can now free up the memory allocated for the msg
    free(msg);
}

/*
    Receives a de-replication database from process source
    and merges it with the local de-replication database db

    Inputs:
        db: pointer to the local de-replication database structure
        my_rank: this process rank
        source: the rank of the source process
*/
 void _recv_derep_db(derep_db* db, int my_rank, int source){
    // We will receive two messages
    // The first one contains the size of second message
    // The second one contains the de-replication db
    int msg_size;
    // Receive the first message
    MPI_Recv(&msg_size, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Create the buffer for receiving the second message
    char* msg = (char*) malloc(sizeof(char) * msg_size);
    // Receive the second message
    MPI_Recv(msg, msg_size, MPI_PACKED, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // Merge the foreign database with the local one while unpacking it
    // This saves memory since we don't really need a derep_db structure
    // Unpack the count and unique counters
    int position = 0;
    int count;
    int unique;
    MPI_Unpack(msg, msg_size, &position, &count, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack(msg, msg_size, &position, &unique, 1, MPI_INT, MPI_COMM_WORLD);
    // We know that all the sequences present in the foreign database are going
    // to be added to the local one - so we can add already the count to it
    db->count += count;
    // Loop through all the unique sequences
    int i;
    int j;
    int length;
    int label_count;
    for(i = 0; i < unique; i++){
        // Unpack the sequence length
        MPI_Unpack(msg, msg_size, &position, &length, 1, MPI_INT, MPI_COMM_WORLD);
        // Allocate memory for the sequence
        char* sequence;
        posix_memalign((void **) &sequence, 16, sizeof(char) * (length+1));
        // Unpack the sequence
        MPI_Unpack(msg, msg_size, &position, sequence, length, MPI_CHAR, MPI_COMM_WORLD);
        sequence[length] = '\0';
        // Check if the sequence is already present on the database
        seq_replicas* r;
        HASH_FIND_STR(db->seqs, sequence, r);
        if(!r){
            // The sequence didn't exist, add as a new sequence
            r = create_empty_seq_replica(sequence, length);
            HASH_ADD_KEYPTR(hh, db->seqs, r->sequence, length, r);
            // Update unique counter
            ++db->unique;
        }
        // Free up sequence memory
        free(sequence);
        // Unpack the number of labels
        MPI_Unpack(msg, msg_size, &position, &label_count, 1, MPI_INT, MPI_COMM_WORLD);
        // Loop through all the labels
        for(j = 0; j < label_count; j++){
            // Unpack the label length
            MPI_Unpack(msg, msg_size, &position, &length, 1, MPI_INT, MPI_COMM_WORLD);
            // Allocate memory for the label
            char* label = (char*) malloc(sizeof(char) * (length+1));
            // Unpack the label
            MPI_Unpack(msg, msg_size, &position, label, length, MPI_CHAR, MPI_COMM_WORLD);
            label[length] = '\0';
            // Add the label to the replica structure
            add_replica(r, label);
            // Free up label memory
            free(label);
        }
    }
    // Free up buffer memory
    free(msg);
}

/*******************************************
* De-replication database public functions *
*******************************************/

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
        add_replica(r, seq->label);
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

/*
    Collects all the information about the derep_db spread across multiple
    processes in the process with rank 0

    Inputs:
        db: the local derep_db - will be modified in place
        my_rank: process rank
        comm_sz: the number of processes
*/
void gather_derep_db(derep_db* db, int my_rank, int comm_sz){
    // Initialize bit mask
    int bit_mask = 0x01 << (int)log2(comm_sz - 1);
    // Loop while the bit_mask is not 0
    while(bit_mask){
        // Get the rank of the partner process
        int partner = my_rank ^ bit_mask;
        if(my_rank & bit_mask){
            // I have a one on the position pointed by the bit mask
            // I am a sender
            _send_derep_db(db, my_rank, partner);
            // Once I sent my data, I'm done
            break;
        }
        else{
            // I have a zero on the position pointed by the bit mask
            // I am a receiver
            // In the first round, it's possible that not all receivers
            // have a sender, so check the partner exists
            if(partner < comm_sz){
                // Receive and merge the foreign de-replication db with
                // the local one
                _recv_derep_db(db, my_rank, partner);
            }
        }
        // Update the bitmask
        bit_mask = bit_mask >> 1;
    }
}