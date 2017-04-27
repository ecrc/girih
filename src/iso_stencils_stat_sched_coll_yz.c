#include "data_structures.h"

#define U(i,j,k)         (u[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])
#define V(i,j,k)         (v[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])
#define ROC2(i,j,k)   (roc2[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])
#define COEF(m,i,j,k) (coef[((k)*(nny)+(j))*(nnx)+(i)+((ln_domain)*(m))])

// ISO stencil 8th-order-in-space-2nd-order-in-time with constant coefficient
void iso_ref( const int shape[3], const int xb, const int yb, const int zb, const int xe, const int ye, const int ze,
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2, stencil_CTX stencil_ctx) {
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i,j,k, jb, je;
  int nny =shape[1];
  int nnx =shape[0];
  
  int  n, nz, nj;
  nz = ze -zb;

  real_t lap;


  for(jb=yb; jb<ye; jb+=stencil_ctx.bs_y) // blocking in Y
  {
    je = (jb+stencil_ctx.bs_y)<(ye)?(jb+stencil_ctx.bs_y):(ye);
#pragma omp parallel num_threads(stencil_ctx.thread_group_size)
    {
      
      nj = je-jb;
#pragma omp for private(n,k,j,i,lap) schedule(static)
      for(n=0; n<nz*nj; n++){
          k = n/nj;
          j = n - nj*k;
          k+= zb;
          j+= jb;
#pragma simd
          for(i=xb; i<xe; i++) {

            lap=coef[0]*V(i,j,k)
            +coef[1]*(V(i+1,j  ,k  )+V(i-1,j  ,k  ))
            +coef[1]*(V(i  ,j+1,k  )+V(i  ,j-1,k  ))
            +coef[1]*(V(i  ,j  ,k+1)+V(i  ,j  ,k-1))
            +coef[2]*(V(i+2,j  ,k  )+V(i-2,j  ,k  ))
            +coef[2]*(V(i  ,j+2,k  )+V(i  ,j-2,k  ))
            +coef[2]*(V(i  ,j  ,k+2)+V(i  ,j  ,k-2))
            +coef[3]*(V(i+3,j  ,k  )+V(i-3,j  ,k  ))
            +coef[3]*(V(i  ,j+3,k  )+V(i  ,j-3,k  ))
            +coef[3]*(V(i  ,j  ,k+3)+V(i  ,j  ,k-3))
            +coef[4]*(V(i+4,j  ,k  )+V(i-4,j  ,k  ))
            +coef[4]*(V(i  ,j+4,k  )+V(i  ,j-4,k  ))
            +coef[4]*(V(i  ,j  ,k+4)+V(i  ,j  ,k-4));

#if DP
            U(i,j,k) = 2.*V(i,j,k) - U(i,j,k) + ROC2(i,j,k)*lap;
#else
            U(i,j,k) = 2.f*V(i,j,k) - U(i,j,k) + ROC2(i,j,k)*lap;
#endif
          }
        }
      
    }
  }
}


// ISO stencil 2nd-order-in-space-1st-order-in-time with constant coefficient
void iso_ref_2space_1time( const int shape[3], const int xb, const int yb, const int zb, const int xe, const int ye, const int ze,
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2, stencil_CTX stencil_ctx) {
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i,j,k, jb, je;
  int nny =shape[1];
  int nnx =shape[0];

  int  n, nz, nj;
  nz = ze -zb;

  for(jb=yb; jb<ye; jb+=stencil_ctx.bs_y) // blocking in Y
  {
    je = (jb+stencil_ctx.bs_y)<(ye)?(jb+stencil_ctx.bs_y):(ye);
#pragma omp parallel num_threads(stencil_ctx.thread_group_size)
    {
      
      nj = je-jb;
#pragma omp for private(n,k,j,i) schedule(static)
      for(n=0; n<nz*nj; n++){
          k = n/nj;
          j = n - nj*k;
          k+= zb;
          j+= jb;
#pragma simd
          for(i=xb; i<xe; i++) {

            U(i,j,k) = coef[0]*V(i,j,k)
                      +coef[1]*(V(i+1,j  ,k  )+V(i-1,j  ,k  ))
                      +coef[1]*(V(i  ,j+1,k  )+V(i  ,j-1,k  ))
                      +coef[1]*(V(i  ,j  ,k+1)+V(i  ,j  ,k-1));
          }
        }
      
    }
  }

}


// ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
void iso_ref_2space_1time_var( const int shape[3], const int xb, const int yb, const int zb, const int xe, const int ye, const int ze,
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2, stencil_CTX stencil_ctx) {
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i,j,k, jb, je;
  int nny =shape[1];
  int nnx =shape[0];
  uint64_t ln_domain = ((uint64_t) 1)* shape[0]*shape[1]*shape[2];

  int  n, nz, nj;
  nz = ze -zb;

  for(jb=yb; jb<ye; jb+=stencil_ctx.bs_y) // blocking in Y
  {
    je = (jb+stencil_ctx.bs_y)<(ye)?(jb+stencil_ctx.bs_y):(ye);
#pragma omp parallel num_threads(stencil_ctx.thread_group_size)
    {
      
      nj = je-jb;
#pragma omp for private(n,k,j,i) schedule(static)
      for(n=0; n<nz*nj; n++){
          k = n/nj;
          j = n - nj*k;
          k+= zb;
          j+= jb;
#pragma simd
          for(i=xb; i<xe; i++) {

            U(i,j,k) = COEF(0,i,j,k)*V(i,j,k)
                      +COEF(1,i,j,k)*(V(i+1,j  ,k  )+V(i-1,j  ,k  ))
                      +COEF(1,i,j,k)*(V(i  ,j+1,k  )+V(i  ,j-1,k  ))
                      +COEF(1,i,j,k)*(V(i  ,j  ,k+1)+V(i  ,j  ,k-1));
          }
      }
      
    }
  }

}


// ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
void iso_ref_2space_1time_var_axsym( const int shape[3], const int xb, const int yb, const int zb, const int xe, const int ye, const int ze,
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2, stencil_CTX stencil_ctx) {
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i,j,k, jb, je;
  int nny =shape[1];
  int nnx =shape[0];
  uint64_t ln_domain = ((uint64_t) 1)* shape[0]*shape[1]*shape[2];

  int  n, nz, nj;
  nz = ze -zb;

  for(jb=yb; jb<ye; jb+=stencil_ctx.bs_y) // blocking in Y
  {
    je = (jb+stencil_ctx.bs_y)<(ye)?(jb+stencil_ctx.bs_y):(ye);
#pragma omp parallel num_threads(stencil_ctx.thread_group_size)
    {
      
      nj = je-jb;
#pragma omp for private(n,k,j,i) schedule(static)
      for(n=0; n<nz*nj; n++){
          k = n/nj;
          j = n - nj*k;
          k+= zb;
          j+= jb;
#pragma simd
          for(i=xb; i<xe; i++) {

            U(i,j,k) = COEF(0,i,j,k)*V(i,j,k)
                      +COEF(1,i,j,k)*(V(i+1,j  ,k  )+V(i-1,j  ,k  ))
                      +COEF(2,i,j,k)*(V(i  ,j+1,k  )+V(i  ,j-1,k  ))
                      +COEF(3,i,j,k)*(V(i  ,j  ,k+1)+V(i  ,j  ,k-1));
          }
      }
      
    }
  }
}


// ISO stencil 8th-order-in-space-1st-order-in-time with variable axis symmetric coefficient
void iso_ref_8space_1time_var_axsym( const int shape[3], const int xb, const int yb, const int zb, const int xe, const int ye, const int ze,
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2, stencil_CTX stencil_ctx) {
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i,j,k, jb, je;
  int nny =shape[1];
  int nnx =shape[0];
  uint64_t ln_domain = ((uint64_t) 1)* shape[0]*shape[1]*shape[2];

  int  n, nz, nj;
  nz = ze -zb;

  for(jb=yb; jb<ye; jb+=stencil_ctx.bs_y) // blocking in Y
  {
    je = (jb+stencil_ctx.bs_y)<(ye)?(jb+stencil_ctx.bs_y):(ye);
#pragma omp parallel num_threads(stencil_ctx.thread_group_size)
    {
      
      nj = je-jb;
#pragma omp for private(n,k,j,i) schedule(static)
      for(n=0; n<nz*nj; n++){
          k = n/nj;
          j = n - nj*k;
          k+= zb;
          j+= jb;
#pragma simd
          for(i=xb; i<xe; i++) {
            U(i,j,k) = COEF(0 ,i,j,k)*V(i,j,k)
                      +COEF(1 ,i,j,k)*(V(i+1,j  ,k  )+V(i-1,j  ,k  ))
                      +COEF(2 ,i,j,k)*(V(i  ,j+1,k  )+V(i  ,j-1,k  ))
                      +COEF(3 ,i,j,k)*(V(i  ,j  ,k+1)+V(i  ,j  ,k-1))
                      +COEF(4 ,i,j,k)*(V(i+2,j  ,k  )+V(i-2,j  ,k  ))
                      +COEF(5 ,i,j,k)*(V(i  ,j+2,k  )+V(i  ,j-2,k  ))
                      +COEF(6 ,i,j,k)*(V(i  ,j  ,k+2)+V(i  ,j  ,k-2))
                      +COEF(7 ,i,j,k)*(V(i+3,j  ,k  )+V(i-3,j  ,k  ))
                      +COEF(8 ,i,j,k)*(V(i  ,j+3,k  )+V(i  ,j-3,k  ))
                      +COEF(9 ,i,j,k)*(V(i  ,j  ,k+3)+V(i  ,j  ,k-3))
                      +COEF(10,i,j,k)*(V(i+4,j  ,k  )+V(i-4,j  ,k  ))
                      +COEF(11,i,j,k)*(V(i  ,j+4,k  )+V(i  ,j-4,k  ))
                      +COEF(12,i,j,k)*(V(i  ,j  ,k+4)+V(i  ,j  ,k-4));
          }
      }
      
    }
  }
}


// ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
void iso_ref_2space_1time_var_nosym( const int shape[3], const int xb, const int yb, const int zb, const int xe, const int ye, const int ze,
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2, stencil_CTX stencil_ctx) {
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i,j,k, jb, je;
  int nny =shape[1];
  int nnx =shape[0];
  uint64_t ln_domain = ((uint64_t) 1)* shape[0]*shape[1]*shape[2];

  int  n, nz, nj;
  nz = ze -zb;

  for(jb=yb; jb<ye; jb+=stencil_ctx.bs_y) // blocking in Y
  {
    je = (jb+stencil_ctx.bs_y)<(ye)?(jb+stencil_ctx.bs_y):(ye);
#pragma omp parallel num_threads(stencil_ctx.thread_group_size)
    {
      
      nj = je-jb;
#pragma omp for private(n,k,j,i) schedule(static)
      for(n=0; n<nz*nj; n++){
          k = n/nj;
          j = n - nj*k;
          k+= zb;
          j+= jb;
#pragma simd
          for(i=xb; i<xe; i++) {

            U(i,j,k) = COEF(0,i,j,k)*V(i,j,k)
                      +COEF(1,i,j,k)*V(i-1,j  ,k  )
                      +COEF(2,i,j,k)*V(i+1,j  ,k  )
                      +COEF(3,i,j,k)*V(i  ,j-1,k  )
                      +COEF(4,i,j,k)*V(i  ,j+1,k  )
                      +COEF(5,i,j,k)*V(i  ,j  ,k-1)
                      +COEF(6,i,j,k)*V(i  ,j  ,k+1);
          }
      }
      
    }
  }
}
