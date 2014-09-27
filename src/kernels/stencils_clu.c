#include "data_structures.h"

#define U(i)         (u_s[(i)])
#define V(i,j,k)     (v_s[((k)*(clu_ctx.nny)+(j))*(clu_ctx.nnx)+(i)])
#define ROC2(i)      (roc2_s[(i)])
#define COEF(m,i,j,k) (coef[((k)*(clu_ctx.nny)+(j))*(clu_ctx.nnx)+(i)+(((unsigned long) (clu_ctx.ln_domain))*(m))])

// 1-stide of ISO stencil 8th-order-in-space-2nd-order-in-time with constant coefficient
void iso_ref_clu_8space_2time CLU_SIG{

  int i;
  FLOAT_PRECISION two=2.0;
  FLOAT_PRECISION lap;
  FLOAT_PRECISION * restrict u_s    =    &u[(clu_ctx.j + clu_ctx.k*clu_ctx.nny)*clu_ctx.nnx];
  const FLOAT_PRECISION * restrict v_s    =    &v[(clu_ctx.j + clu_ctx.k*clu_ctx.nny)*clu_ctx.nnx];
  const FLOAT_PRECISION * restrict roc2_s = &roc2[(clu_ctx.j + clu_ctx.k*clu_ctx.nny)*clu_ctx.nnx];

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

    U(i) = two*V(i,0,0) - U(i) + ROC2(i)*lap;
  }
}


// 1-stride of ISO stencil 2nd-order-in-space-1st-order-in-time with constant coefficient
void iso_ref_clu_2space_1time CLU_SIG{

  int i;
  FLOAT_PRECISION * restrict u_s       = &u[(clu_ctx.j + clu_ctx.k*clu_ctx.nny)*clu_ctx.nnx];
  const FLOAT_PRECISION * restrict v_s = &v[(clu_ctx.j + clu_ctx.k*clu_ctx.nny)*clu_ctx.nnx];

#pragma simd
  for(i=clu_ctx.xb; i<clu_ctx.xe; i++) {
    U(i) = coef[0]*V(i,0,0)
          +coef[1]*(V(i+1,0  ,0  )+V(i-1,0  ,0  ))
          +coef[1]*(V(i  , +1,0  )+V(i  , -1,0  ))
          +coef[1]*(V(i  ,0  , +1)+V(i  ,0  , -1));
  }
}


// 1-stride of ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
void iso_ref_clu_2space_1time_var CLU_SIG{

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = clu_ctx.j*clu_ctx.nnx + clu_ctx.k*nxny;
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
void iso_ref_clu_2space_1time_var_axsym CLU_SIG{

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = clu_ctx.j*clu_ctx.nnx + clu_ctx.k*nxny;
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


void iso_ref_clu_8space_1time_var_axsym CLU_SIG{

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = clu_ctx.j*clu_ctx.nnx + clu_ctx.k*nxny;
  FLOAT_PRECISION * restrict u_s       = &u[stride_start];
  const FLOAT_PRECISION * restrict v_s = &v[stride_start];
  const FLOAT_PRECISION * restrict coef0_s = &coef[stride_start];
  const FLOAT_PRECISION * restrict coef1_s = &coef[stride_start + clu_ctx.ln_domain];
  const FLOAT_PRECISION * restrict coef2_s = &coef[stride_start + clu_ctx.ln_domain*2];
  const FLOAT_PRECISION * restrict coef3_s = &coef[stride_start + clu_ctx.ln_domain*3];
  const FLOAT_PRECISION * restrict coef4_s = &coef[stride_start + clu_ctx.ln_domain*4];
  const FLOAT_PRECISION * restrict coef5_s = &coef[stride_start + clu_ctx.ln_domain*5];
  const FLOAT_PRECISION * restrict coef6_s = &coef[stride_start + clu_ctx.ln_domain*6];
  const FLOAT_PRECISION * restrict coef7_s = &coef[stride_start + clu_ctx.ln_domain*7];
  const FLOAT_PRECISION * restrict coef8_s = &coef[stride_start + clu_ctx.ln_domain*8];
  const FLOAT_PRECISION * restrict coef9_s = &coef[stride_start + clu_ctx.ln_domain*9];
  const FLOAT_PRECISION * restrict coef10_s= &coef[stride_start + clu_ctx.ln_domain*10];
  const FLOAT_PRECISION * restrict coef11_s= &coef[stride_start + clu_ctx.ln_domain*11];
  const FLOAT_PRECISION * restrict coef12_s= &coef[stride_start + clu_ctx.ln_domain*12];

#pragma simd
  for(i=clu_ctx.xb; i<clu_ctx.xe; i++) {
    U(i) = coef0_s[i]*V(i,0,0)
          +coef1_s[i]*(V(i+1,0  ,0  )+V(i-1,0  ,0  ))
          +coef2_s[i]*(V(i  , +1,0  )+V(i  , -1,0  ))
          +coef3_s[i]*(V(i  ,0  , +1)+V(i  ,0  , -1))
          +coef4_s[i]*(V(i+2,0  ,0  )+V(i-2,0  ,0  ))
          +coef5_s[i]*(V(i  , +2,0  )+V(i  , -2,0  ))
          +coef6_s[i]*(V(i  ,0  , +2)+V(i  ,0  , -2))
          +coef7_s[i]*(V(i+3,0  ,0  )+V(i-3,0  ,0  ))
          +coef8_s[i]*(V(i  , +3,0  )+V(i  , -3,0  ))
          +coef9_s[i]*(V(i  ,0  , +3)+V(i  ,0  , -3))
         +coef10_s[i]*(V(i+4,0  ,0  )+V(i-4,0  ,0  ))
         +coef11_s[i]*(V(i  , +4,0  )+V(i  , -4,0  ))
         +coef12_s[i]*(V(i  ,0  , +4)+V(i  ,0  , -4));


  }
}

// 1-stride of ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
void iso_ref_clu_2space_1time_var_nosym CLU_SIG{

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = clu_ctx.j*clu_ctx.nnx + clu_ctx.k*nxny;
  FLOAT_PRECISION * restrict u_s       = &u[stride_start];
  const FLOAT_PRECISION * restrict v_s = &v[stride_start];
  const FLOAT_PRECISION * restrict coef0_s = &coef[stride_start];
  const FLOAT_PRECISION * restrict coef1_s = &coef[stride_start + clu_ctx.ln_domain];
  const FLOAT_PRECISION * restrict coef2_s = &coef[stride_start + clu_ctx.ln_domain*2];
  const FLOAT_PRECISION * restrict coef3_s = &coef[stride_start + clu_ctx.ln_domain*3];
  const FLOAT_PRECISION * restrict coef4_s = &coef[stride_start + clu_ctx.ln_domain*4];
  const FLOAT_PRECISION * restrict coef5_s = &coef[stride_start + clu_ctx.ln_domain*5];
  const FLOAT_PRECISION * restrict coef6_s = &coef[stride_start + clu_ctx.ln_domain*6];
#pragma simd
  for(i=clu_ctx.xb; i<clu_ctx.xe; i++) {
    U(i) = coef0_s[i]*V(i  ,0 ,0 )
          +coef1_s[i]*V(i-1,0 ,0 )
          +coef2_s[i]*V(i+1,0 ,0 )
          +coef3_s[i]*V(i  ,-1,0 )
          +coef4_s[i]*V(i  ,+1,0 )
          +coef5_s[i]*V(i  ,0 ,-1)
          +coef6_s[i]*V(i  ,0 ,+1);
  }
}

clu_func_t clu_func_list[] = {
    iso_ref_clu_8space_2time,
    iso_ref_clu_2space_1time,
    iso_ref_clu_2space_1time_var,
    iso_ref_clu_2space_1time_var_axsym,
    iso_ref_clu_8space_1time_var_axsym,
    iso_ref_clu_2space_1time_var_nosym,
};
