//
//  MPI ring communication
//
//  Yifei Sun
//  03/05/2015
//

#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include "util.h"

int main( int argc, char *argv[])
{
    
    if (argc != 2) {
        fprintf(stderr, "Function needs to specify the number of iterations as an input argument!\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    int rank, tag, origin;
    int p;  // number of processors
    long i, N, iter=0;
    int sum = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Status status;
    
    N = atol(argv[1]);  // number of iterations
    
    char hostname[1024];
    gethostname(hostname, 1024);
    
    int message_out, message_in=0;
    tag = 99;
    
    
    
    timestamp_type time1, time2;
    get_timestamp(&time1);
    
    do {
        iter++;
        
        if(rank == 0)
        {
            message_in = message_in + rank;
            MPI_Send(&message_in, 1, MPI_INT, 1, tag, MPI_COMM_WORLD);
            //printf("rank %d sent %d to rank %d\n", rank, message_in, rank+1);
            MPI_Recv(&message_in, 1, MPI_INT, p-1, tag, MPI_COMM_WORLD, &status);
        }
        else
        {
            MPI_Recv(&message_in, 1, MPI_INT, rank-1, tag, MPI_COMM_WORLD, &status);
            message_in = message_in + rank;
            MPI_Send(&message_in, 1, MPI_INT, (rank+1)%p, tag, MPI_COMM_WORLD);
            //printf("rank %d sent %d to rank %d\n", rank, message_in, (rank+1)%p);
        }

        
    }  while (iter < N);

    get_timestamp(&time2);
    // compute the time elapsed
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    printf("Time elapsed is %f seconds.\n", elapsed);
    
    MPI_Finalize();
    return 0;
}