#include "data_structures.h"

#define U(i)         (u_s[(i)])
#define V(i,j,k)     (v_s[((k)*(clu_ctx.nny)+(j))*(clu_ctx.nnx)+(i)])
#define ROC2(i)      (roc2_s[(i)])
#define COEF(m,i,j,k) (coef[((k)*(clu_ctx.nny)+(j))*(clu_ctx.nnx)+(i)+(((unsigned long) (clu_ctx.ln_domain))*(m))])


// ISO stencils
void iso_ref_split( const int shape[3], const int xb, const int yb, const int zb, const int xe, const int ye, const int ze,
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2, stencil_CTX stencil_ctx) {
  int j,k, jb, je;

  CLU_CTX clu_ctx;
  clu_ctx.nnx = shape[0];
  clu_ctx.nny = shape[1];
  clu_ctx.nnz = shape[2];
  clu_ctx.xb = xb;
  clu_ctx.xe = xe;
  clu_ctx.ln_domain = shape[0]*shape[1]*shape[2];

  for(jb=yb; jb<ye; jb+=stencil_ctx.bs_y) // blocking in Y
  {
    je = (jb+stencil_ctx.bs_y)<(ye)?(jb+stencil_ctx.bs_y):(ye);
#pragma omp parallel num_threads(stencil_ctx.thread_group_size)
    {
#pragma omp for private(k,j) schedule(static,1)
      for(k=zb; k<ze; k++) {
        for(j=jb; j<je; j++) {
          stencil_ctx.ref_stride(clu_ctx, j, k, coef, u, v, roc2);
        }
      }
    }
  }
}

// 1-stide of ISO stencil 8th-order-in-space-2nd-order-in-time with constant coefficient
void iso_ref_8space_2time_stride STRIDE1_SIG{

  int i;
  FLOAT_PRECISION lap;
  FLOAT_PRECISION * restrict u_s    =    &u[(j + k*clu_ctx.nny)*clu_ctx.nnx];
  const FLOAT_PRECISION * restrict v_s    =    &v[(j + k*clu_ctx.nny)*clu_ctx.nnx];
  const FLOAT_PRECISION * restrict roc2_s = &roc2[(j + k*clu_ctx.nny)*clu_ctx.nnx];

#pragma simd
  for(i=clu_ctx.xb; i<clu_ctx.xe; i++) {

    lap=coef[0]*V(i,0,0)                    // MOVE THE MUL OUT
    +coef[1]*(V(i+1,0  ,0  )+V(i-1,0  ,0  ))
    +coef[1]*(V(i  , +1,0  )+V(i  , -1,0  ))
    +coef[1]*(V(i  ,0  , +1)+V(i  ,0  , -1))
    +coef[2]*(V(i+2,0  ,0  )+V(i-2,0  ,0  ))
    +coef[2]*(V(i  , +2,0  )+V(i  , -2,0  ))
    +coef[2]*(V(i  ,0  , +2)+V(i  ,0  , -2))
    +coef[3]*(V(i+3,0  ,0  )+V(i-3,0  ,0  ))
    +coef[3]*(V(i  , +3,0  )+V(i  , -3,0  ))
    +coef[3]*(V(i  ,0  , +3)+V(i  ,0  , -3))
    +coef[4]*(V(i+4,0  ,0  )+V(i-4,0  ,0  ))
    +coef[4]*(V(i  , +4,0  )+V(i  , -4,0  ))
    +coef[4]*(V(i  ,0  , +4)+V(i  ,0  , -4));

#if DP
    U(i) = 2.*V(i,0,0) - U(i) + ROC2(i)*lap;
#else
    U(i) = 2.f*V(i,0,0) - U(i) + ROC2(i)*lap;
#endif
  }
}


// 1-stride of ISO stencil 2nd-order-in-space-1st-order-in-time with constant coefficient
void iso_ref_2space_1time_stride STRIDE1_SIG{

  int i;
  FLOAT_PRECISION * restrict u_s       = &u[(j + k*clu_ctx.nny)*clu_ctx.nnx];
  const FLOAT_PRECISION * restrict v_s = &v[(j + k*clu_ctx.nny)*clu_ctx.nnx];

#pragma simd
  for(i=clu_ctx.xb; i<clu_ctx.xe; i++) {
    U(i) = coef[0]*V(i,0,0)
          +coef[1]*(V(i+1,0  ,0  )+V(i-1,0  ,0  ))
          +coef[1]*(V(i  , +1,0  )+V(i  , -1,0  ))
          +coef[1]*(V(i  ,0  , +1)+V(i  ,0  , -1));
  }
}


// 1-stride of ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
void iso_ref_2space_1time_var_stride STRIDE1_SIG{

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = j*clu_ctx.nnx + k*nxny;
  FLOAT_PRECISION * restrict u_s       = &u[stride_start];
  const FLOAT_PRECISION * restrict v_s = &v[stride_start];
  const FLOAT_PRECISION * restrict coef0_s = &coef[stride_start];
  const FLOAT_PRECISION * restrict coef1_s = &coef[stride_start + clu_ctx.ln_domain];
#pragma simd
  for(i=clu_ctx.xb; i<clu_ctx.xe; i++) {
    U(i) = coef0_s[i]*V(i,0,0)
          +coef1_s[i]*(V(i+1,0  ,0  )+V(i-1,0  ,0  ))
          +coef1_s[i]*(V(i  , +1,0  )+V(i  , -1,0  ))
          +coef1_s[i]*(V(i  ,0  , +1)+V(i  ,0  , -1));
  }
}

// 1-stride of ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
void iso_ref_2space_1time_var_axsym_stride STRIDE1_SIG{

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = j*clu_ctx.nnx + k*nxny;
  FLOAT_PRECISION * restrict u_s       = &u[stride_start];
  const FLOAT_PRECISION * restrict v_s = &v[stride_start];
  const FLOAT_PRECISION * restrict coef0_s = &coef[stride_start];
  const FLOAT_PRECISION * restrict coef1_s = &coef[stride_start + clu_ctx.ln_domain];
  const FLOAT_PRECISION * restrict coef2_s = &coef[stride_start + clu_ctx.ln_domain*2];
  const FLOAT_PRECISION * restrict coef3_s = &coef[stride_start + clu_ctx.ln_domain*3];
#pragma simd
  for(i=clu_ctx.xb; i<clu_ctx.xe; i++) {
    U(i) = coef0_s[i]*V(i,0,0)
          +coef1_s[i]*(V(i+1,0  ,0  )+V(i-1,0  ,0  ))
          +coef2_s[i]*(V(i  , +1,0  )+V(i  , -1,0  ))
          +coef3_s[i]*(V(i  ,0  , +1)+V(i  ,0  , -1));
  }
}
