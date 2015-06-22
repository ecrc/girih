/*
 * Order-1, 3D 7 point stencil
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
    Nx = atoi(argv[1])+2;
    Ny = atoi(argv[2])+2;
    Nz = atoi(argv[3])+2;
  }
  if (argc > 4)
    Nt = atoi(argv[4]);


  double ****A = (double ****) malloc(sizeof(double***)*2);
  A[0] = (double ***) malloc(sizeof(double**)*Nz);
  A[1] = (double ***) malloc(sizeof(double**)*Nz);
  for(i=0; i<Nz; i++){
    A[0][i] = (double**) malloc(sizeof(double*)*Ny);
    A[1][i] = (double**) malloc(sizeof(double*)*Ny);
    for(j=0;j<Ny;j++){
      A[0][i][j] = (double*) malloc(sizeof(double)*Nx);
      A[1][i][j] = (double*) malloc(sizeof(double)*Nx);
    }
  }

  // tile size information, including extra element to decide the list length
  int *tile_size = (int*) malloc(sizeof(int));
  tile_size[0] = -1;
  // The list is modified here before source-to-source transformations 
  tile_size = (int*) realloc((void *)tile_size, sizeof(int)*5);
  tile_size[0] = 8;
  tile_size[1] = 8;
  tile_size[2] = 24;
  tile_size[3] = 64;
  tile_size[4] = -1;


  // for timekeeping
  int ts_return = -1;
  struct timeval start, end, result;
  double tdiff = 0.0, min_tdiff=1.e100;


  const int BASE = 1024;
  const double alpha = 0.0876;
  const double beta = 0.0765;


  // initialize variables
  //
  srand(42);
    for (i = 1; i < Nz; i++) {
        for (j = 1; j < Ny; j++) {
            for (k = 1; k < Nx; k++) {
                A[0][i][j][k] = 1.0 * (rand() % BASE);
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


  for(test=0; test<TESTS; test++){
      gettimeofday(&start, 0);
  // serial execution - Addition: 6 && Multiplication: 2
#pragma scop
      for (t = 0; t < Nt-1; t++) {
          for (i = 1; i < Nz-1; i++) {
              for (j = 1; j < Ny-1; j++) {
                  for (k = 1; k < Nx-1; k++) {

                      A[(t+1)%2][i][j][k] = alpha * (A[t%2][i][j][k])
                          + beta * (A[t%2][i - 1][j][k] + A[t%2][i][j - 1][k] + A[t%2][i][j][k - 1] +
                                  A[t%2][i + 1][j][k] + A[t%2][i][j + 1][k] + A[t%2][i][j][k + 1]);

                  }
              }
          }
      }
#pragma endscop
      gettimeofday(&end, 0);
      ts_return = timeval_subtract(&result, &end, &start);
      tdiff = (double) (result.tv_sec + result.tv_usec * 1.0e-6);
      min_tdiff = min(min_tdiff, tdiff);
      printf("Rank 0 TEST# %d time: %f\n", test, tdiff);

  }

  PRINT_RESULTS(1, "constant")

#ifdef LIKWID_PERFMON
  #pragma omp parallel
  {
      LIKWID_MARKER_STOP("calc");
  }
  LIKWID_MARKER_CLOSE;
#endif

  // Free allocated arrays (Causing performance degradation
/* for(i=0; i<Nz; i++){
    for(j=0;j<Ny;j++){
      free(A[0][i][j]);
      free(A[1][i][j]);
    }
    free(A[0][i]);
    free(A[1][i]);
  }
  free(A[0]);
  free(A[1]);
*/
  return 0;
}

