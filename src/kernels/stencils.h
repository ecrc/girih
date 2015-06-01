
#ifndef STENCILS_H_
#define STENCILS_H_


#include "data_structures.h"
#include <stdlib.h>
#include <stdio.h>
//#include <xmmintrin.h>

#define U(i,j,k)         (u[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])
#define V(i,j,k)         (v[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])
#define ROC2(i,j,k)   (roc2[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])
#define COEF(m,i,j,k) (coef[((k)*(nny)+(j))*(nnx)+(i)+((ln_domain)*(m))])

#define CAT(X,Y) X##_##Y
#define TEMPLATE(X,Y) CAT(X,Y)

#define BLOCK_COND(x) while((x)){asm("");}

//function for the unsupported features
void not_supported_mwd KERNEL_MWD_SIG{
  printf("ERROR: unsupported configuration for the selected stencil\n");
  exit(1); 
}
void not_supported KERNEL_SIG{
  printf("ERROR: unsupported configuration for the selected stencil\n");
  exit(1); 
}

// Source: http://stackoverflow.com/questions/13772567/get-cpu-cycle-count
static inline uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

#endif /* STENCILS_H_ */
