#include "data_structures.h"

#define U(i)         (u_s[(i)])
#define V(i,j,k)     (v_s[((k)*(clu_ctx.nny)+(j))*(clu_ctx.nnx)+(i)])
#define ROC2(i)      (roc2_s[(i)])
#define COEF(m,i,j,k) (coef[((k)*(clu_ctx.nny)+(j))*(clu_ctx.nnx)+(i)+(((uint64_t) (clu_ctx.ln_domain))*(m))])

// 1-stide of ISO stencil 8th-order-in-space-2nd-order-in-time with constant coefficient
void iso_ref_clu_8space_2time CLU_SIG{

  int i;
  real_t two=2.0;
  real_t lap;
  real_t * restrict u_s    =    &u[(j + k*clu_ctx.nny)*clu_ctx.nnx];
  const real_t * restrict v_s    =    &v[(j + k*clu_ctx.nny)*clu_ctx.nnx];
  const real_t * restrict roc2_s = &roc2[(j + k*clu_ctx.nny)*clu_ctx.nnx];

  // TODO Add the source contribution
  real_t customroc = 16; //4000^2 0.001^2   @KADIR
  printf("diamond. ");
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

#pragma simd
  for(i=xb; i<xe; i++) {

    lap=3.*coef[0]*V(i,0,0)                    // MOVE THE MUL OUT
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

    //U(i) = two*V(i,0,0) - U(i) + ROC2(i)*lap; //@KADIR
    U(i) = two*V(i,0,0) - U(i) + customroc*lap/400.; //@KADIR
  }
}


// 1-stride of ISO stencil 2nd-order-in-space-1st-order-in-time with constant coefficient
void iso_ref_clu_2space_1time CLU_SIG{
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i;
  real_t * restrict u_s       = &u[(j + k*clu_ctx.nny)*clu_ctx.nnx];
  const real_t * restrict v_s = &v[(j + k*clu_ctx.nny)*clu_ctx.nnx];

#pragma simd
  for(i=xb; i<xe; i++) {
    U(i) = coef[0]*V(i,0,0)
          +coef[1]*(V(i+1,0  ,0  )+V(i-1,0  ,0  ))
          +coef[1]*(V(i  , +1,0  )+V(i  , -1,0  ))
          +coef[1]*(V(i  ,0  , +1)+V(i  ,0  , -1));
  }
}


// 1-stride of ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
void iso_ref_clu_2space_1time_var CLU_SIG{
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = j*clu_ctx.nnx + k*nxny;
  real_t * restrict u_s       = &u[stride_start];
  const real_t * restrict v_s = &v[stride_start];
  const real_t * restrict coef0_s = &coef[stride_start];
  const real_t * restrict coef1_s = &coef[stride_start + clu_ctx.ln_domain];
#pragma simd
  for(i=xb; i<xe; i++) {
    U(i) = coef0_s[i]*V(i,0,0)
          +coef1_s[i]*(V(i+1,0  ,0  )+V(i-1,0  ,0  ))
          +coef1_s[i]*(V(i  , +1,0  )+V(i  , -1,0  ))
          +coef1_s[i]*(V(i  ,0  , +1)+V(i  ,0  , -1));
  }
}

// 1-stride of ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
void iso_ref_clu_2space_1time_var_axsym CLU_SIG{
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = j*clu_ctx.nnx + k*nxny;
  real_t * restrict u_s       = &u[stride_start];
  const real_t * restrict v_s = &v[stride_start];
  const real_t * restrict coef0_s = &coef[stride_start];
  const real_t * restrict coef1_s = &coef[stride_start + clu_ctx.ln_domain];
  const real_t * restrict coef2_s = &coef[stride_start + clu_ctx.ln_domain*2];
  const real_t * restrict coef3_s = &coef[stride_start + clu_ctx.ln_domain*3];
#pragma simd
  for(i=xb; i<xe; i++) {
    U(i) = coef0_s[i]*V(i,0,0)
          +coef1_s[i]*(V(i+1,0  ,0  )+V(i-1,0  ,0  ))
          +coef2_s[i]*(V(i  , +1,0  )+V(i  , -1,0  ))
          +coef3_s[i]*(V(i  ,0  , +1)+V(i  ,0  , -1));
  }
}


void iso_ref_clu_8space_1time_var_axsym CLU_SIG{
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = j*clu_ctx.nnx + k*nxny;
  real_t * restrict u_s       = &u[stride_start];
  const real_t * restrict v_s = &v[stride_start];
  const real_t * restrict coef0_s = &coef[stride_start];
  const real_t * restrict coef1_s = &coef[stride_start + clu_ctx.ln_domain];
  const real_t * restrict coef2_s = &coef[stride_start + clu_ctx.ln_domain*2];
  const real_t * restrict coef3_s = &coef[stride_start + clu_ctx.ln_domain*3];
  const real_t * restrict coef4_s = &coef[stride_start + clu_ctx.ln_domain*4];
  const real_t * restrict coef5_s = &coef[stride_start + clu_ctx.ln_domain*5];
  const real_t * restrict coef6_s = &coef[stride_start + clu_ctx.ln_domain*6];
  const real_t * restrict coef7_s = &coef[stride_start + clu_ctx.ln_domain*7];
  const real_t * restrict coef8_s = &coef[stride_start + clu_ctx.ln_domain*8];
  const real_t * restrict coef9_s = &coef[stride_start + clu_ctx.ln_domain*9];
  const real_t * restrict coef10_s= &coef[stride_start + clu_ctx.ln_domain*10];
  const real_t * restrict coef11_s= &coef[stride_start + clu_ctx.ln_domain*11];
  const real_t * restrict coef12_s= &coef[stride_start + clu_ctx.ln_domain*12];

#pragma simd
  for(i=xb; i<xe; i++) {
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
  printf("\n%s %s %d  %d\n", __FILE__, __func__, __LINE__, __COUNTER__);

  int i;
  int nxny = clu_ctx.nnx*clu_ctx.nny;
  int stride_start = j*clu_ctx.nnx + k*nxny;
  real_t * restrict u_s       = &u[stride_start];
  const real_t * restrict v_s = &v[stride_start];
  const real_t * restrict coef0_s = &coef[stride_start];
  const real_t * restrict coef1_s = &coef[stride_start + clu_ctx.ln_domain];
  const real_t * restrict coef2_s = &coef[stride_start + clu_ctx.ln_domain*2];
  const real_t * restrict coef3_s = &coef[stride_start + clu_ctx.ln_domain*3];
  const real_t * restrict coef4_s = &coef[stride_start + clu_ctx.ln_domain*4];
  const real_t * restrict coef5_s = &coef[stride_start + clu_ctx.ln_domain*5];
  const real_t * restrict coef6_s = &coef[stride_start + clu_ctx.ln_domain*6];
#pragma simd
  for(i=xb; i<xe; i++) {
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
