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
    int rank;
    int local_n;
    int npes;
    
    // load data 
    if ((ip = fopen("data_input_meta","r")) == NULL) {
        printf("Error opening the data_input_meta file.\n");
        return 253;
    }
    fscanf(ip, "%d\n", &nodecount);
    fclose(ip);
    if (node_init(&nodehead, 0, nodecount)) return 254;
    
    // initialize variables
    r = malloc(nodecount * sizeof(double));
    r_pre = malloc(nodecount * sizeof(double));
    double *global_r = malloc(nodecount * sizeof(double));
    
    GET_TIME(start);
    iterationcount = 0;
    for ( i = 0; i < nodecount; ++i)
        r[i] = 1.0 / nodecount;
    /* INITIALIZE MORE VARIABLES IF NECESSARY */

    //MPI initialization
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    local_n = nodecount/npes;
    int *recvcounts = malloc(npes * sizeof(int));
    int *disp = malloc(npes * sizeof(int));
    int total = 0;

    // Partitioning per process
    for (i = 0; i < npes; i++) {
        if (i < nodecount % npes) {
            recvcounts[i] = local_n + 1;
        }
        else {
            recvcounts[i] = local_n;
        }
        disp[i] = total;
        total = total + recvcounts[i];
    }

    int p_start = disp[rank];
    int p_end = p_start + recvcounts[rank];

    // core calculation
    do{
        ++iterationcount;
        /* IMPLEMENT ITERATIVE UPDATE */
        vec_cp(r, r_pre, nodecount);
        for (i = p_start; i < p_end; i++) {
            struct node i_node = nodehead[i];
            double summation = 0;

            for (j = 0; j < i_node.num_in_links; j++) {
                summation += r_pre[i_node.inlinks[j]]/(nodehead[i_node.inlinks[j]].num_out_links);
            }
            r[i] = ((1 - DAMPING_FACTOR)/nodecount) + DAMPING_FACTOR * summation;
        }

        MPI_Allgatherv(r + p_start, p_end - p_start, MPI_DOUBLE, global_r, recvcounts, disp, MPI_DOUBLE, MPI_COMM_WORLD);
        
    }while(rel_error(r, r_pre, nodecount) >= EPSILON);
    MPI_Finalize();
    GET_TIME(end);

    Lab4_saveoutput(r, nodecount, end - start);

    // post processing
    node_destroy(nodehead, nodecount);
    free(r); free(r_pre); free(global_r);
    return 0;
}
