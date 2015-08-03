#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

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
  tile_size[0] = 16;
  tile_size[1] = 16;
  tile_size[2] = 8;
  tile_size[3] = 256;
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
/* Copyright (C) 1991-2014 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */
/* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((Nt >= 1) && (Nx >= 9) && (Ny >= 9) && (Nz >= 9)) {
  for (t1=-1;t1<=floord(Nt-1,2);t1++) {
    lbp=max(ceild(t1,2),ceild(4*t1-Nt+2,4));
    ubp=min(floord(4*Nt+Nz-9,16),floord(8*t1+Nz+2,16));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(max(max(0,ceild(16*t2-Nz+5,8)),t1),2*t1-2*t2+1);t3<=min(min(min(floord(4*Nt+Ny-9,8),floord(8*t1+Ny+7,8)),floord(16*t2+Ny+3,8)),floord(16*t1-16*t2+Nz+Ny+5,8));t3++) {
        for (t4=max(max(max(0,ceild(t1-31,32)),ceild(16*t2-Nz-243,256)),ceild(8*t3-Ny-243,256));t4<=min(min(min(min(floord(4*Nt+Nx-9,256),floord(8*t1+Nx+7,256)),floord(16*t2+Nx+3,256)),floord(8*t3+Nx-5,256)),floord(16*t1-16*t2+Nz+Nx+5,256));t4++) {
          for (t5=max(max(max(max(max(0,ceild(16*t2-Nz+5,4)),ceild(8*t3-Ny+5,4)),ceild(256*t4-Nx+5,4)),2*t1),4*t1-4*t2+1);t5<=min(min(min(min(min(floord(16*t1-16*t2+Nz+10,4),2*t3),Nt-1),2*t1+3),4*t2+2),64*t4+62);t5++) {
            for (t6=max(max(16*t2,4*t5+4),-16*t1+16*t2+8*t5-15);t6<=min(min(16*t2+15,-16*t1+16*t2+8*t5),4*t5+Nz-5);t6++) {
              for (t7=max(8*t3,4*t5+4);t7<=min(8*t3+7,4*t5+Ny-5);t7++) {
                lbv=max(256*t4,4*t5+4);
                ubv=min(256*t4+255,4*t5+Nx-5);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[( t5 + 1) % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8)] = (((2.0 * A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8)]) - A[( t5 + 1) % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8)]) + (roc2[ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8)] * (((((coef0 * A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8)]) + (coef1 * (((((A[ t5 % 2][ (-4*t5+t6) - 1][ (-4*t5+t7)][ (-4*t5+t8)] + A[ t5 % 2][ (-4*t5+t6) + 1][ (-4*t5+t7)][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7) - 1][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7) + 1][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8) - 1]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8) + 1]))) + (coef2 * (((((A[ t5 % 2][ (-4*t5+t6) - 2][ (-4*t5+t7)][ (-4*t5+t8)] + A[ t5 % 2][ (-4*t5+t6) + 2][ (-4*t5+t7)][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7) - 2][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7) + 2][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8) - 2]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8) + 2]))) + (coef3 * (((((A[ t5 % 2][ (-4*t5+t6) - 3][ (-4*t5+t7)][ (-4*t5+t8)] + A[ t5 % 2][ (-4*t5+t6) + 3][ (-4*t5+t7)][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7) - 3][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7) + 3][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8) - 3]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8) + 3]))) + (coef4 * (((((A[ t5 % 2][ (-4*t5+t6) - 4][ (-4*t5+t7)][ (-4*t5+t8)] + A[ t5 % 2][ (-4*t5+t6) + 4][ (-4*t5+t7)][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7) - 4][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7) + 4][ (-4*t5+t8)]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8) - 4]) + A[ t5 % 2][ (-4*t5+t6)][ (-4*t5+t7)][ (-4*t5+t8) + 4])))));;
                }
              }
            }
          }
        }
      }
    }
  }
}
/* End of CLooG code */
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

