
#ifndef STENCILS_H_
#define STENCILS_H_


#include "data_structures.h"
#include <stdlib.h>

#define U(i,j,k)         (u[((k)*(nny)+(j))*(nnx)+(i)])
#define V(i,j,k)         (v[((k)*(nny)+(j))*(nnx)+(i)])
#define ROC2(i,j,k)   (roc2[((k)*(nny)+(j))*(nnx)+(i)])
#define COEF(m,i,j,k) (coef[((k)*(nny)+(j))*(nnx)+(i)+(((unsigned long) (ln_domain))*(m))])

#define CAT(X,Y) X##_##Y
#define TEMPLATE(X,Y) CAT(X,Y)

#endif /* STENCILS_H_ */
