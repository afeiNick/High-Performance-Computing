//
//  This is the parallel version of the Jacobi iteration
//  
//  Yifei Sun
//  03/05/2015
//


#include <stdio.h>
/* timing in util.h requires -lrt flag to compile */
#include "util.h"
#include <unistd.h>
#include <mpi.h>

int main (int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "Function needs to specify the number of grid points and iterations as input arguments!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    

    long i, N, n, iter=0, iter_total;
    int p;
    double *u_new, *u_old, *f;
    double h;    // h is the grid size
    int rank, tag;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    MPI_Status status;

    N = atol(argv[1]);  // number of interior points
    iter_total = atol(argv[2]);  // number of iterations
    
    if ((N%p) != 0) {
        fprintf(stderr, "The number of grid points needs to be divisible by the number of cores!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    

    n = N/p + 2;
    h = 1./(N+1);
    

    
    char hostname[1024];
    gethostname(hostname, 1024);
    


    
    double message_to_left, message_to_right;
    double message_from_left, message_from_right;
    
    tag = 99;
    
    u_new = (double *) malloc(sizeof(double) * (n));
    u_old = (double *) malloc(sizeof(double) * (n));
    f = (double *) malloc(sizeof(double) * (n));
    
    // initialize
    for (i = 0; i < n; ++i)  {
        u_old[i] = 0.;
        u_new[i] = 0.;
        f[i] = 1.;
    }
    
    timestamp_type time1, time2;
    get_timestamp(&time1);
    
    do {
        
        iter++;
        
        // one Jacobi iteration
        if (rank == 0)
        {
            message_to_right = u_old[n-2];
            MPI_Send(&message_to_right, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
            MPI_Recv(&message_from_right, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
            u_old[n-1] = message_from_right;
            
            for (i = 1; i < n-1; ++i) {
                u_new[i] = 0.5*h*h*( 1 + (u_old[i-1]+u_old[i+1])/(h*h) );
            }
            
            for (i = 1; i < n-1; ++i) {
                u_old[i] = u_new[i];
            }
        }
        else if (rank == p-1)
        {
            message_to_left = u_old[1];
            MPI_Send(&message_to_left, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
            
            MPI_Recv(&message_from_left, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
            u_old[0] = message_from_left;
            
            for (i = 1; i < n-1; ++i) {
                u_new[i] = 0.5*h*h*( 1 + (u_old[i-1]+u_old[i+1])/(h*h) );
            }
            
            for (i = 1; i < n-1; ++i) {
                u_old[i] = u_new[i];
            }
        }
        else
        {
            message_to_left = u_old[1];
            message_to_right = u_old[n-2];
            MPI_Send(&message_to_left, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD);
            MPI_Send(&message_to_right, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD);
            
            MPI_Recv(&message_from_left, 1, MPI_DOUBLE, rank-1, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&message_from_right, 1, MPI_DOUBLE, rank+1, tag, MPI_COMM_WORLD, &status);
            u_old[0] = message_from_left;
            u_old[n-1] = message_from_right;
            
            for (i = 1; i < n-1; ++i) {
                u_new[i] = 0.5*h*h*( 1 + (u_old[i-1]+u_old[i+1])/(h*h) );
            }
            
            for (i = 1; i < n-1; ++i) {
                u_old[i] = u_new[i];
            }
        }
        
        
    }  while (iter < iter_total);
    
    
    double first_elem = u_new[1];
    printf("the first array element on rank %d hosted on %s is %.12f\n", rank, hostname, first_elem);

    get_timestamp(&time2);
    // compute the time elapsed
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    printf("Time elapsed is %f seconds.\n", elapsed);

    free(u_new);
    free(u_old);
    free(f);
    
    MPI_Finalize();
    return 0;
}
