#include "data_structures.h"

// ISO stencil
void iso_ref_all_wf_split( const int shape[3], const int xb, const int yb_r, const int zb, const int xe, const int ye_r, const int ze,
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2, int t_dim, int b_inc, int e_inc, int NHALO, stencil_CTX stencil_ctx, int mtid) {
  int j, k, t, yb, ye, zi, yb_r2;

  int yb_m, ye_m;
  int time_blk = t_dim*2+1; //temporal block size
  yb_m = yb_r - 2*t_dim*b_inc;
  ye_m = ye_r + 2*t_dim*e_inc;

  CLU_CTX clu_ctx;
  clu_ctx.nnx = shape[0];
  clu_ctx.nny = shape[1];
  clu_ctx.nnz = shape[2];
  clu_ctx.xb = xb;
  clu_ctx.xe = xe;
  clu_ctx.ln_domain = shape[0]*shape[1]*shape[2];

  for(zi=zb; zi<ze; zi++) {
    #pragma omp parallel for private(t,j,yb,ye,k) schedule(static,1) //num_threads(stencil_ctx.thread_group_size)
    for(t=0; t< time_blk; t++){

      if(t <= t_dim){ // inverted trapezoid (or lower half of the diamond)
        yb = yb_r - t*b_inc;
        ye = ye_r + t*e_inc;
      }else{ // trapezoid  (or upper half of the diamond)
        yb = yb_m + t*b_inc;
        ye = ye_m - t*e_inc;
      }
      k = zi - t*(NHALO+1);

      if((t)%2 == 1){

        for(j=yb; j<ye; j++) {
          stencil_ctx.ref_stride(clu_ctx, j, k, coef, u, v, roc2);
        }

      }else{

        for(j=yb; j<ye; j++) {
          stencil_ctx.ref_stride(clu_ctx, j, k, coef, v, u, roc2);
        }
      }
    }
  }
}
