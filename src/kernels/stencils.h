
#ifndef STENCILS_H_
#define STENCILS_H_


#include "data_structures.h"
#include <stdlib.h>
#include <stdio.h>

#define U(i,j,k)         (u[((k)*(nny)+(j))*(nnx)+(i)])
#define V(i,j,k)         (v[((k)*(nny)+(j))*(nnx)+(i)])
#define ROC2(i,j,k)   (roc2[((k)*(nny)+(j))*(nnx)+(i)])
#define COEF(m,i,j,k) (coef[((k)*(nny)+(j))*(nnx)+(i)+(((unsigned long) (ln_domain))*(m))])

#define CAT(X,Y) X##_##Y
#define TEMPLATE(X,Y) CAT(X,Y)


//function for the unsupported features
void not_supported_mwd KERNEL_MWD_SIG{
  printf("ERROR: unsupported configuration for the selected stencil\n");
  exit(1); 
}
void not_supported KERNEL_SIG{
  printf("ERROR: unsupported configuration for the selected stencil\n");
  exit(1); 
}
#endif /* STENCILS_H_ */
