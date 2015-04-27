/* Parallel sample sort
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include "util.h"


static int compare(const void *a, const void *b)
{
  int *da = (int *)a;
  int *db = (int *)b;

  if (*da > *db)
    return 1;
  else if (*da < *db)
    return -1;
  else
    return 0;
}

int main( int argc, char *argv[])
{
  int rank;
  int i, j, N;
  int *vec;
  int p;    // the number of processors to use
    
  int *splitters;   // this array will contain the selected splitters from vec
  int num_splitters;  // how many splitters to select from vec, i.e. p-1
  int *combined_splitters;  // this will be on the root to store the combined splitters

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  /* Number of random numbers per processor (this should be increased
   * for actual tests or could be passed in through the command line */
  N = 100;
  num_splitters = p-1;

  vec = calloc(N, sizeof(int));
  splitters = calloc(num_splitters, sizeof(int));
    
  /* seed random number generator differently on every core */
  srand((unsigned int) (rank + 393919));

  /* fill vector with random integers */
  for (i = 0; i < N; ++i) {
    vec[i] = rand();
  }
  printf("rank: %d, first entry: %d\n", rank, vec[0]);

  timestamp_type time1, time2;
  get_timestamp(&time1);

  /* sort locally */
  qsort(vec, N, sizeof(int), compare);


  /* randomly sample s entries from vector or select local splitters,
   * i.e., every N/P-th entry of the sorted vector */
  for (i = 0; i < num_splitters; ++i) {
    splitters[i] = vec[(i+1) * (N/p)];
  }
  
  /* every processor communicates the selected entries
   * to the root processor; use for instance an MPI_Gather */
  combined_splitters = calloc(num_splitters*p, sizeof(int));
  MPI_Gather(splitters, num_splitters, MPI_INT, combined_splitters, num_splitters, MPI_INT, 0, MPI_COMM_WORLD);
  
  /* root processor does a sort, determinates splitters that
   * split the data into P buckets of approximately the same size */

  if (rank == 0) {
    qsort(combined_splitters, num_splitters*p, sizeof(int), compare);
    // we then pick (p-1) splitters from combined_splitters to put into global_splitters
    for (i = 0; i < num_splitters; ++i) {
      splitters[i] = combined_splitters[(i+1) * num_splitters];
    }
  }
    
  /* root process broadcasts splitters */
  MPI_Bcast(splitters, p-1, MPI_INT, 0, MPI_COMM_WORLD);
  

  /* every processor uses the obtained splitters to decide
   * which integers need to be sent to which other processor (local bins) */
    
  int *send_count, *receive_count;
  int *send_position, *receive_position;
  send_count = calloc(p, sizeof(int));    // corresponding to p bins
  receive_count = calloc(p, sizeof(int));
  send_position = calloc(p, sizeof(int));
  receive_position = calloc(p, sizeof(int));
    
  for (i = 0; i < N; ++i) {
    if (vec[i] <= splitters[0]) {
      send_count[0] = send_count[0] + 1;
    }
    else if (vec[i] > splitters[p-2]) {
      send_count[p-1] = send_count[p-1] + 1;
    }
    else {
      for (j = 1; j < p-1; ++j) {
        if ((vec[i] > splitters[j-1]) && (vec[i] <= splitters[j])) {
          send_count[j] = send_count[j] + 1;
        }
      }
    }
  }
  
  for (i = 0; i < 1; ++i) {
    printf("rank: %d entry %d\n", rank, send_count[i]);
  }

  /* send and receive: either you use MPI_AlltoallV, or
   * (and that might be easier), use an MPI_Alltoall to share
   * with every processor how many integers it should expect,
   * and then use MPI_Send and MPI_Recv to exchange the data */
  
  MPI_Alltoall(send_count, 1, MPI_INT, receive_count, 1, MPI_INT, MPI_COMM_WORLD);
  
  int *vec_new;   // after communication the new vector length is called vec_new
  int N_new = 0;  // length of vec_new
  
  for (i = 0; i < p; ++i) {
    N_new = N_new + receive_count[i];
  }
  
  vec_new = calloc(N_new, sizeof(int));

  for (i = 1; i < p; ++i) {
    send_position[i] = send_count[i-1] + send_position[i-1];
    receive_position[i] = receive_count[i-1] + receive_position[i-1];
  }
  
  MPI_Alltoallv(vec, send_count, send_position, MPI_INT,
                vec_new, receive_count, receive_position, MPI_INT, MPI_COMM_WORLD);
  
  /* do a local sort */
  qsort(vec_new, N_new, sizeof(int), compare);
  

  get_timestamp(&time2);
  // compute the time elapsed
  double elapsed = timestamp_diff_in_seconds(time1,time2);

  /* every processor writes its result to a file */
  FILE* fd = NULL;
  char filename[256];
  snprintf(filename, 256, "output%02d.txt", rank);
  fd = fopen(filename,"w+");
  
  if(NULL == fd)
  {
    printf("Error opening file \n");
    return 1;
  }
  
  fprintf(fd, "rank %d sorted list: \n", rank);
  for(i = 0; i < N_new; ++i)
    fprintf(fd, "  %d\n", vec_new[i]);
  fprintf(fd, "Time elapsed is %f seconds.\n", elapsed);
  
  fclose(fd);
  // writing file done
  

  free(vec);
  free(splitters);
  free(combined_splitters);
  free(send_count);
  free(receive_count);
  free(send_position);
  free(receive_position);
  free(vec_new);
  
  MPI_Finalize();
  return 0;
}
