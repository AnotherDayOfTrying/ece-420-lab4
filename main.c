/*
    Serial Implementation of Lab 4
*/

#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "Lab4_IO.h"
#include "timer.h"

#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

int main (int argc, char* argv[]){
    // instantiate variables
    struct node *nodehead;
    int nodecount;
    double *r, *r_pre;
    int i, j;
    int iterationcount;
    double start, end;
    FILE *ip;
    /* INSTANTIATE MORE VARIABLES IF NECESSARY */
    int rank; // process id (in MPI)
    int npes; // number of processes (in MPI)
    
    // load data 
    if ((ip = fopen("data_input_meta","r")) == NULL) {
        printf("Error opening the data_input_meta file.\n");
        return 253;
    }
    fscanf(ip, "%d\n", &nodecount);
    fclose(ip);
    if (node_init(&nodehead, 0, nodecount)) return 254;
    
    // initialize variables
    r = malloc(nodecount * sizeof(double)); // current vector of page ranks
    r_pre = malloc(nodecount * sizeof(double)); // previous vector of page ranks
    double *global_r = malloc(nodecount * sizeof(double)); // vector of page ranks (used for receival)
    double *summation = malloc(nodecount * sizeof(double)); // vector of page ranks (used for eq3 summation)
    
    
    iterationcount = 0; // number of iterations
    double START_PROB = 1.0 / nodecount; // initial probability for each node
    for ( i = 0; i < nodecount; ++i)
        r[i] = START_PROB; // start with 1/N probability for each node
    /* INITIALIZE MORE VARIABLES IF NECESSARY */

    // MPI initialization
    MPI_Init(&argc, &argv); // Initialize the MPI execution. May be used to pass the command line arguments to all processes
    MPI_Comm_size(MPI_COMM_WORLD, &npes); // var holds the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // var holds current process id (different between procs)

    int local_n = nodecount/npes; // initial (floored) number of nodes per process
    int *part_sizes = malloc(npes * sizeof(int)); // array for holding the number of nodes per process
    int *displacement = malloc(npes * sizeof(int)); // array for holding the index displacement of each process
    
    // Partitioning per process (with uneven distribution, in case nodecount is not divisible by npes)
    int total = 0;
    for (i = 0; i < npes; i++) {
        if (i < nodecount % npes) {
            part_sizes[i] = local_n + 1; 
        }
        else {
            part_sizes[i] = local_n;
        }
        displacement[i] = total;
        total = total + part_sizes[i];
    }

    int p_start = displacement[rank]; // start index of this processor's partition
    int p_end = p_start + part_sizes[rank]; // end index of this processor's partition
    
    const double DAMP_START_PROB = (1.0 - DAMPING_FACTOR)/nodecount; // const used in first half of eq 3
    GET_TIME(start); // only measure the time required to run the PageRank algorithm and the necessary distributed computing overhead

    // core calculation
    do{
        ++iterationcount;
        /* IMPLEMENT ITERATIVE UPDATE */
        vec_cp(r, r_pre, nodecount);

        // perform the summation in eq 3
        for (i = p_start; i < p_end; i++) {
            struct node i_node = nodehead[i];
            r[i] = 0;
            for (j = 0; j < i_node.num_in_links; j++) {
                r[i] += summation[i_node.inlinks[j]];
            }
            r[i] += DAMP_START_PROB;
        }

        // gather the updated r values from all processes
        MPI_Allgatherv(r + p_start, p_end - p_start, MPI_DOUBLE, global_r, part_sizes, displacement, MPI_DOUBLE, MPI_COMM_WORLD);
        
        // update local_r with global_r, and calculate the summation values for the next iteration
        for (int i = 0; i < nodecount; i++) { 
            r[i] = global_r[i];
            summation[i] = DAMPING_FACTOR * r[i] / nodehead[i].num_out_links;
        }
        
    }while(rel_error(r, r_pre, nodecount) >= EPSILON); // continue until the error is less than EPSILON
    MPI_Finalize(); // end the shared execution environment
    GET_TIME(end); // only measure the time required to run the PageRank algorithm and the necessary distributed computing overhead

    Lab4_saveoutput(r, nodecount, end - start);

    // post processing
    node_destroy(nodehead, nodecount);
    free(r); free(r_pre); free(global_r); free(summation);
    return 0;
}
