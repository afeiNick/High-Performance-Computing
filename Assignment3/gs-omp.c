//
//  The OpenMP version of Jacobi iteration for solving Laplace's equation
//  
//
//  Created by Yifei Sun on 4/22/15
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "util.h"


int main (int argc, char **argv)
{
    if (argc != 3) {
        fprintf(stderr, "Function needs to specify the number of grid points (even) and iterations as input arguments!\n");
        abort();
    }
    
    
    long i, N, iter=0, iter_total;
    double *u, *f;
    double h;
    
    N = atol(argv[1]);  // number of interior points
    iter_total = atol(argv[2]);  // number of iterations
    h = 1./(N+1);       // grid size
    
    u = (double *) malloc(sizeof(double) * (N+2));
    f = (double *) malloc(sizeof(double) * (N+2));
    
    timestamp_type time1, time2;
    get_timestamp(&time1);
    
    // starting the parallel section
    // we must make the loop index variable private
    #pragma omp parallel private(i) shared(u,f)
    {
        // fill vectors
        // parallel for loop
        #pragma omp for schedule(static)
        for (i = 0; i <= N+1; i++) {
            u[i] = 0.;
            f[i] = 1.;
        }
        /* an implicit barrier */
        
        do {
            iter++;
        
            // here we are assuming N is even for simplicity
            // update the red (even) points
            #pragma omp for schedule(static)
            for (i = 1; i <= N/2; ++i) {
                u[2*i] = 0.5*h*h*( 1 + u[2*i-1]/(h*h) + u[2*i+1]/(h*h) );
            }
            
            // update the black (odd) points
            #pragma omp for schedule(static)
            for (i = 1; i <= N/2; ++i) {
                u[2*i-1] = 0.5*h*h*( 1 + u[2*i-2]/(h*h) + u[2*i]/(h*h) );
            }
        }
        while (iter < iter_total);
    }  // end of the parallel section

    get_timestamp(&time2);
    double elapsed = timestamp_diff_in_seconds(time1,time2);
    
    printf("Time elapsed is %f seconds.\n", elapsed);
    printf("The function value at the middle point is %e.\n", u[N/2]);
    
    free(u);
    free(f);
    return 0;
}