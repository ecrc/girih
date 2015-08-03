/*
 * Order-2, 3D 25 point stencil
 * Adapted from PLUTO and Pochoir test bench
 *
 * Tareq Malas
 */

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#ifdef LIKWID_PERFMON
#include <likwid.h>
#endif
#include "print_utils.h"
#define TESTS 2

#define MAX(a,b) ((a) > (b) ? a : b)
#define MIN(a,b) ((a) < (b) ? a : b)

#ifndef min
#define min(x,y)    ((x) < (y)? (x) : (y))
#endif

/* Subtract the `struct timeval' values X and Y,
 * storing the result in RESULT.
 *
 * Return 1 if the difference is negative, otherwise 0.
 */
int timeval_subtract(struct timeval *result, struct timeval *x, struct timeval *y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec)
  {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;

    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }

  if (x->tv_usec - y->tv_usec > 1000000)
  {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;

    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
   * tv_usec is certainly positive.
   */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}

int main(int argc, char *argv[])
{
  int t, i, j, k, test;
  int Nx, Ny, Nz, Nt;
  if (argc > 3) {
    Nx = atoi(argv[1])+8;
    Ny = atoi(argv[2])+8;
    Nz = atoi(argv[3])+8;
  }
  if (argc > 4)
    Nt = atoi(argv[4]);


  double ****A = (double ****) malloc(sizeof(double***)*2);
  double ***roc2 = (double ***) malloc(sizeof(double**));
  A[0] = (double ***) malloc(sizeof(double**)*Nz);
  A[1] = (double ***) malloc(sizeof(double**)*Nz);
  roc2 = (double ***) malloc(sizeof(double**)*Nz);
  for(i=0; i<Nz; i++){
    A[0][i] = (double**) malloc(sizeof(double*)*Ny);
    A[1][i] = (double**) malloc(sizeof(double*)*Ny);
    roc2[i] = (double**) malloc(sizeof(double*)*Ny);
    for(j=0;j<Ny;j++){
      A[0][i][j] = (double*) malloc(sizeof(double)*Nx);
      A[1][i][j] = (double*) malloc(sizeof(double)*Nx);
      roc2[i][j] = (double*) malloc(sizeof(double)*Nx);
    }
  }

  // tile size information, including extra element to decide the list length
  int *tile_size = (int*) malloc(sizeof(int));
  tile_size[0] = -1;
  // The list is modified here before source-to-source transformations 
  tile_size = (int*) realloc((void *)tile_size, sizeof(int)*5);
  tile_size[0] = 24;
  tile_size[1] = 24;
  tile_size[2] = 24;
  tile_size[3] = 2048;
  tile_size[4] = -1;


  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0, min_tdiff=1.e100;


  const int BASE = 1024;


  // initialize variables
  //
  srand(42);
    for (i = 1; i < Nz; i++) {
        for (j = 1; j < Ny; j++) {
            for (k = 1; k < Nx; k++) {
                A[0][i][j][k] = 1.0 * (rand() % BASE);
                roc2[i][j][k] = 2.0 * (rand() % BASE);
            }
        }
    }


#ifdef LIKWID_PERFMON
  LIKWID_MARKER_INIT;
  #pragma omp parallel
  {
      LIKWID_MARKER_THREADINIT;
  #pragma omp barrier
      LIKWID_MARKER_START("calc");
  }
#endif

  int num_threads = 1;
#if defined(_OPENMP)
  num_threads = omp_get_max_threads();
#endif

  const double coef0 = -0.28472;
  const double coef1 = 0.16000;
  const double coef2 = -0.02000;
  const double coef3 = 0.00254;
  const double coef4 = -0.00018;

  for(test=0; test<TESTS; test++){
      gettimeofday(&start, 0);
  // serial execution - Addition: 6 && Multiplication: 2
#pragma scop
      for (t = 0; t < Nt; t++) {
          for (i = 4; i < Nz-4; i++) {
              for (j = 4; j < Ny-4; j++) {
                  for (k = 4; k < Nx-4; k++) {

                      A[(t+1)%2][i][j][k] = 2.0*A[t%2][i][j][k] - A[(t+1)%2][i][j][k] + roc2[i][j][k]*( 
                        coef0* A[t%2][i  ][j  ][k  ] +
                        coef1*(A[t%2][i-1][j  ][k  ] + A[t%2][i+1][j  ][k  ]  +
                               A[t%2][i  ][j-1][k  ] + A[t%2][i  ][j+1][k  ]  +
                               A[t%2][i  ][j  ][k-1] + A[t%2][i  ][j  ][k+1]) +
                        coef2*(A[t%2][i-2][j  ][k  ] + A[t%2][i+2][j  ][k  ]  +
                               A[t%2][i  ][j-2][k  ] + A[t%2][i  ][j+2][k  ]  +
                               A[t%2][i  ][j  ][k-2] + A[t%2][i  ][j  ][k+2]) +
                        coef3*(A[t%2][i-3][j  ][k  ] + A[t%2][i+3][j  ][k  ]  +
                               A[t%2][i  ][j-3][k  ] + A[t%2][i  ][j+3][k  ]  +
                               A[t%2][i  ][j  ][k-3] + A[t%2][i  ][j  ][k+3]) +
                        coef4*(A[t%2][i-4][j  ][k  ] + A[t%2][i+4][j  ][k  ]  +
                               A[t%2][i  ][j-4][k  ] + A[t%2][i  ][j+4][k  ]  +
                               A[t%2][i  ][j  ][k-4] + A[t%2][i  ][j  ][k+4]) );
                  }
              }
          }
      }
#pragma endscop
      gettimeofday(&end, 0);
      ts_return = timeval_subtract(&result, &end, &start);
      tdiff = (double) (result.tv_sec + result.tv_usec * 1.0e-6);
      min_tdiff = MIN(min_tdiff, tdiff);
      printf("Rank 0 TEST# %d time: %f\n", test, tdiff);

  }

  PRINT_RESULTS(4, "constant")

#ifdef LIKWID_PERFMON
  #pragma omp parallel
  {
      LIKWID_MARKER_STOP("calc");
  }
  LIKWID_MARKER_CLOSE;
#endif

  // Free allocated arrays
 for(i=0; i<Nz; i++){
    for(j=0;j<Ny;j++){
      free(A[0][i][j]);
      free(A[1][i][j]);
      free(roc2[i][j]);
    }
    free(A[0][i]);
    free(A[1][i]);
    free(roc2[i]);
  }
  free(A[0]);
  free(A[1]);
  free(roc2);
 
  return 0;
}

