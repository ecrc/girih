#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include "verification.h"

void verify(Parameters *p){

  verification_printing(*p);
  mpi_halo_init(p);

  // allocate and initialize the required arrays
  arrays_allocate(p);
  init_coeff(p);
  domain_data_fill(p);

  // compute data in parallel using the time stepper to be tested
  TSList[p->target_ts].func(p);

  // aggregate all subdomains into rank zero to compare with the serial results
  FLOAT_PRECISION * restrict aggr_domain = NULL;
  if(p->stencil.type == REGULAR) aggr_domain = aggregate_subdomains(*p);
  // compute data serially using reference time stepper
  //    compare and report results
  if(p->mpi_rank==0) {
    verify_serial_generic(aggr_domain, *p);
    if(p->stencil.type == REGULAR) free(aggr_domain);
  }

  arrays_free(p);
  mpi_halo_finalize(p);
}
void verify_serial_generic(FLOAT_PRECISION * target_domain, Parameters p) {
  int male;
  // function pointer to select the reference stencil operator
  void (*std_kernel)( const int [],
      const FLOAT_PRECISION * restrict, FLOAT_PRECISION * restrict,
      const FLOAT_PRECISION * restrict, const FLOAT_PRECISION * restrict);

  FLOAT_PRECISION * restrict u, * restrict v, * restrict roc2, * restrict coef;
  unsigned long it, i, j , k, f, ax;
  int nnx = p.stencil_shape[0]+2*p.stencil.r;
  int nny = p.stencil_shape[1]+2*p.stencil.r;
  int nnz = p.stencil_shape[2]+2*p.stencil.r;
  unsigned long n_domain = nnx*nny*nnz;
  unsigned long domain_size;

  switch(p.stencil.type){
    case REGULAR:
      male = posix_memalign((void **)&(u),    p.alignment, sizeof(FLOAT_PRECISION)*n_domain); check_merr(male);
      male = posix_memalign((void **)&(v),    p.alignment, sizeof(FLOAT_PRECISION)*n_domain); check_merr(male);
    if(p.stencil.time_order == 2)
      male = posix_memalign((void **)&(roc2), p.alignment, sizeof(FLOAT_PRECISION)*n_domain); check_merr(male);
     break;

    case SOLAR:
      domain_size = n_domain*12lu*2lu;
      male = posix_memalign((void **)&(u), p.alignment, sizeof(FLOAT_PRECISION)*domain_size); check_merr(male);
      v = 0; 
      break;

   default:
    printf("ERROR: unknown type of stencil\n");
    exit(1);
    break;
  }


  // allocate the the coefficients according to the stencil type
  // initialize the coefficients
  // select the desired stencil operator function
  unsigned long coef_size, idx;
  switch(p.stencil.coeff){
  case CONSTANT_COEFFICIENT:
    coef_size = 11;
    male = posix_memalign((void **)&(coef), p.alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male);
    for(i=0;i<p.stencil.r+1;i++) coef[i] = p.g_coef[i];

    // select the stencil type
    switch(p.stencil.r){
    case 1:
      std_kernel = &std_kernel_2space_1time;
      break;
    case 4:
      std_kernel = &std_kernel_8space_2time;
      break;
    default:
      printf("ERROR: unknown type of stencil\n");
      exit(1);
      break;
    }
    break;


  case VARIABLE_COEFFICIENT:
    coef_size = n_domain*(1 + p.stencil.r);
    male = posix_memalign((void **)&(coef), p.alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male);
    for(k=0; k <= p.stencil.r; k++){
      for(i=0; i<n_domain; i++){
        coef[i + k*n_domain] = p.g_coef[k];
      }
    }
    // select the stencil type
    switch(p.stencil.r){
    case 1:
      std_kernel = &std_kernel_2space_1time_var;
      break;
    default:
      printf("ERROR: unknown type of stencil\n");
      exit(1);
      break;
    }
    break;


  case VARIABLE_COEFFICIENT_AXSYM:
    coef_size = n_domain*(1 + 3*p.stencil.r);
    male = posix_memalign((void **)&(coef), p.alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male);

    // central point coeff
    for(i=0; i<n_domain; i++){
      coef[i] = p.g_coef[0];
    }
    for(k=0; k < p.stencil.r; k++){
      for(ax=0; ax<3; ax++){
        for(i=0; i<n_domain; i++){
          coef[i + n_domain + 3*k*n_domain + ax*n_domain] = p.g_coef[k+1];
        }
      }
    }

    // select the stencil type
    switch(p.stencil.r){
    case 1:
      std_kernel = &std_kernel_2space_1time_var_axsym;
      break;
    case 4:
      std_kernel = &std_kernel_8space_1time_var_axsym;
      break;
    default:
      printf("ERROR: unknown type of stencil\n");
      exit(1);
      break;
    }
    break;


  case VARIABLE_COEFFICIENT_NOSYM:
    coef_size = n_domain*(1 + 6*p.stencil.r);
    male = posix_memalign((void **)&(coef), p.alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male);

    // central point coeff
    for(i=0; i<n_domain; i++){
      coef[i] = p.g_coef[0];
    }
    for(k=0; k < p.stencil.r; k++){
      for(ax=0; ax<3; ax++){
        for(i=0; i<n_domain; i++){
          coef[i + n_domain + 6*k*n_domain +  2*ax   *n_domain] = p.g_coef[k+1];
          coef[i + n_domain + 6*k*n_domain + (2*ax+1)*n_domain] = p.g_coef[k+1];
        }
      }
    }

    // select the stencil type
    switch(p.stencil.r){
    case 1:
      std_kernel = &std_kernel_2space_1time_var_nosym;
      break;
    default:
      printf("ERROR: unknown type of stencil\n");
      exit(1);
      break;
    }
    break;


  case SOLAR_COEFFICIENT:
    std_kernel = &solar_kernel;
    coef_size = n_domain*28lu*2lu;
    male = posix_memalign((void **)&(coef), p.alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male); 
    for(f=0; f<28;f++){
      for(idx=0;idx<n_domain*2lu; idx++){
        coef[idx+f*n_domain*2lu] = p.g_coef[(idx+f*n_domain*2lu)%10];
      }
    }
    break;


  default:
    printf("ERROR: unknown type of stencil\n");
    exit(1);
    break;
  }


  switch(p.stencil.type){
    case REGULAR:
      //-- Initialize u,v,roc2
      // fill the array points according to their location in space and pad the boundary with zeroes
      for(i=0; i<n_domain;i++){
        u[i] = 0.0;
        v[i] = 0.0;
        if(p.stencil.time_order == 2)
          roc2[i]= 0.0;
      }
      FLOAT_PRECISION r;
      for(k=p.stencil.r; k<p.stencil_shape[2]+p.stencil.r; k++){
        for(j=p.stencil.r; j<p.stencil_shape[1]+p.stencil.r; j++){
          for(i=p.stencil.r; i<p.stencil_shape[0]+p.stencil.r; i++){
            r = 1.0/3 * (1.0*i/p.stencil_shape[0] + 1.0*j/p.stencil_shape[1]  + 1.0*k/p.stencil_shape[2]);
            U(i, j, k)    = r*1.845703;
            V(i, j, k)    = r*1.845703;

            if(p.stencil.time_order == 2)
              ROC2(i, j, k) = r*1.845703;
          }
        }
      }
    //   set source points at the boundary of the leading dimension
      for(k=0; k<nnz; k++){
        for(j=0; j<nny; j++){
          U(0, j, k) += BOUNDARY_SRC_VAL;
          V(0, j, k) += BOUNDARY_SRC_VAL;
          U(nnx-1, j, k) += BOUNDARY_SRC_VAL;
          V(nnx-1, j, k) += BOUNDARY_SRC_VAL;
        }
      }
      break;
  
    case SOLAR:
      for(i=0; i<n_domain*24lu;i++){
        u[i] = 0.0;
      }
      for(f=0; f<12; f++){
        for(k=0; k<nnz; k++){
          for(j=0; j<nny; j++){
            for(i=0; i<nnx; i++){
              idx = 2*((k*nny+j)*nnx + i +n_domain*f);
              r = 1.0/(3.0) * (1.0*i/p.stencil_shape[0] + 1.0*j/p.stencil_shape[1]  + 1.0*k/p.stencil_shape[2]);
              u[idx] = r*1.845703;
              u[idx+1] = r*1.845703/3.0;
            }
          }
        }
      }
      break;
    
   default:
    printf("ERROR: unknown type of stencil\n");
    exit(1);
    break;
  } 
  //-- Reference kernel Main Loop -
  int domain_shape[3];
  domain_shape[0] = nnx;
  domain_shape[1] = nny;
  domain_shape[2] = nnz;
  for(it=0; it<p.nt; it+=2){
    std_kernel(domain_shape, coef, u, v, roc2);
    std_kernel(domain_shape, coef, v, u, roc2);
  }
//  print_3Darray("u"   , u, nnx, nny, nnz, 0);
//  u[(p.stencil.r+1)*(nnx * nny + nny + 1)] += 100.1;

  // compare results
  compare_results(u, target_domain, p.alignment, p.stencil_shape[0], p.stencil_shape[1], p.stencil_shape[2], p.stencil.r, p);


  switch(p.stencil.type){
    case REGULAR:
      free(u);
      free(v);
      free(coef);
      if(p.stencil.time_order == 2)
        free(roc2);
      break;

    case SOLAR:
      free(u);
      free(coef);
      break;

   default:
    printf("ERROR: unknown type of stencil\n");
    exit(1);
    break;
  }

}

// This is the standard ISO 25-points stencil kernel with constant coefficients,
void std_kernel_8space_2time( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2) {

  int i,j,k;
  int nnz =shape[2];
  int nny =shape[1];
  int nnx =shape[0];

  FLOAT_PRECISION lap;

  for(k=4; k<nnz-4; k++) {
    for(j=4; j<nny-4; j++) {
      for(i=4; i<nnx-4; i++) {

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
        U(i,j,k) = 2. *V(i,j,k) - U(i,j,k) + ROC2(i,j,k)*lap;
#else
        U(i,j,k) = 2.f*V(i,j,k) - U(i,j,k) + ROC2(i,j,k)*lap;
#endif
      }
    }
  }
}

// This is the standard ISO 7-points stencil kernel with constant coefficient
void std_kernel_2space_1time( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2) {

  int i,j,k;
  int nnz =shape[2];
  int nny =shape[1];
  int nnx =shape[0];

  for(k=1; k<nnz-1; k++) {
    for(j=1; j<nny-1; j++) {
      for(i=1; i<nnx-1; i++) {
        U(i,j,k) = coef[0]*V(i,j,k)
                  +coef[1]*(V(i+1,j  ,k  )+V(i-1,j  ,k  ))
                  +coef[1]*(V(i  ,j+1,k  )+V(i  ,j-1,k  ))
                  +coef[1]*(V(i  ,j  ,k+1)+V(i  ,j  ,k-1));
      }
    }
  }

}
// This is the standard ISO 7-points stencil kernel with variable coefficient
void std_kernel_2space_1time_var( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2) {

  int i,j,k;
  int nnz =shape[2];
  int nny =shape[1];
  int nnx =shape[0];
  unsigned long ln_domain = shape[0]*shape[1]*shape[2];

  for(k=1; k<nnz-1; k++) {
    for(j=1; j<nny-1; j++) {
      for(i=1; i<nnx-1; i++) {
        U(i,j,k) = COEF(0,i,j,k)*V(i,j,k)
                  +COEF(1,i,j,k)*(V(i+1,j  ,k  )+V(i-1,j  ,k  ))
                  +COEF(1,i,j,k)*(V(i  ,j+1,k  )+V(i  ,j-1,k  ))
                  +COEF(1,i,j,k)*(V(i  ,j  ,k+1)+V(i  ,j  ,k-1));
      }
    }
  }

}
// This is the standard ISO 7-points stencil kernel with variable coefficient
void std_kernel_2space_1time_var_axsym( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2) {

  int i,j,k;
  int nnz =shape[2];
  int nny =shape[1];
  int nnx =shape[0];
  unsigned long ln_domain = shape[0]*shape[1]*shape[2];

  for(k=1; k<nnz-1; k++) {
    for(j=1; j<nny-1; j++) {
      for(i=1; i<nnx-1; i++) {
        U(i,j,k) = COEF(0,i,j,k)*V(i,j,k)
                  +COEF(1,i,j,k)*(V(i+1,j  ,k  )+V(i-1,j  ,k  ))
                  +COEF(2,i,j,k)*(V(i  ,j+1,k  )+V(i  ,j-1,k  ))
                  +COEF(3,i,j,k)*(V(i  ,j  ,k+1)+V(i  ,j  ,k-1));
      }
    }
  }

}
// This is the standard ISO 25-points stencil kernel with variable axis symmetric coefficients,
void std_kernel_8space_1time_var_axsym( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2) {

  int i,j,k;
  int nnz =shape[2];
  int nny =shape[1];
  int nnx =shape[0];
  unsigned long ln_domain = shape[0]*shape[1]*shape[2];

  for(k=4; k<nnz-4; k++) {
    for(j=4; j<nny-4; j++) {
      for(i=4; i<nnx-4; i++) {

        U(i,j,k) =COEF(0 ,i,j,k)*V(i,j,k)
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
// This is the standard ISO 7-points stencil kernel with variable coefficient
void std_kernel_2space_1time_var_nosym( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2) {

  int i,j,k;
  int nnz =shape[2];
  int nny =shape[1];
  int nnx =shape[0];
  unsigned long ln_domain = shape[0]*shape[1]*shape[2];

  for(k=1; k<nnz-1; k++) {
    for(j=1; j<nny-1; j++) {
      for(i=1; i<nnx-1; i++) {
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

void solar_h_field_ref( const int shape[3], const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u){
  int i,j,k;
  int nnz =shape[2];
  int nny =shape[1];
  int nnx =shape[0];
  unsigned long ln_domain = 2ul*shape[0]*shape[1]*shape[2];
  unsigned long ixmin, ixmax;

  FLOAT_PRECISION * restrict Hyxd = &(u[0ul*ln_domain]);
  FLOAT_PRECISION * restrict Hzxd = &(u[1ul*ln_domain]);
  FLOAT_PRECISION * restrict Hxyd = &(u[2ul*ln_domain]);
  FLOAT_PRECISION * restrict Hzyd = &(u[3ul*ln_domain]);
  FLOAT_PRECISION * restrict Hxzd = &(u[4ul*ln_domain]);
  FLOAT_PRECISION * restrict Hyzd = &(u[5ul*ln_domain]);

  FLOAT_PRECISION * restrict Exzd = &(u[6ul*ln_domain]);
  FLOAT_PRECISION * restrict Eyzd = &(u[7ul*ln_domain]);
  FLOAT_PRECISION * restrict Eyxd = &(u[8ul*ln_domain]);
  FLOAT_PRECISION * restrict Ezxd = &(u[9ul*ln_domain]);
  FLOAT_PRECISION * restrict Exyd = &(u[10ul*ln_domain]);
  FLOAT_PRECISION * restrict Ezyd = &(u[11ul*ln_domain]);


  const FLOAT_PRECISION * restrict cHyxd = &(coef[0ul*ln_domain]);
  const FLOAT_PRECISION * restrict cHzxd = &(coef[1ul*ln_domain]);
  const FLOAT_PRECISION * restrict cHxyd = &(coef[2ul*ln_domain]);
  const FLOAT_PRECISION * restrict cHzyd = &(coef[3ul*ln_domain]);
  const FLOAT_PRECISION * restrict cHxzd = &(coef[4ul*ln_domain]);
  const FLOAT_PRECISION * restrict cHyzd = &(coef[5ul*ln_domain]);

  const FLOAT_PRECISION * restrict tHyxd = &(coef[6ul*ln_domain]);
  const FLOAT_PRECISION * restrict tHzxd = &(coef[7ul*ln_domain]);
  const FLOAT_PRECISION * restrict tHxyd = &(coef[8ul*ln_domain]);
  const FLOAT_PRECISION * restrict tHzyd = &(coef[9ul*ln_domain]);
  const FLOAT_PRECISION * restrict tHxzd = &(coef[10ul*ln_domain]);
  const FLOAT_PRECISION * restrict tHyzd = &(coef[11ul*ln_domain]);

  const FLOAT_PRECISION * restrict Hxbndd = &(coef[12ul*ln_domain]);
  const FLOAT_PRECISION * restrict Hybndd = &(coef[13ul*ln_domain]);

  unsigned long isub;
  FLOAT_PRECISION stagDiffR, stagDiffI, asgn;

  // Update H-field
  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) {
        // Hy_x: csc, Ex_{y,z}: css
        isub      = i + 2 * ( -nnx * nny );
        stagDiffR = Exyd[i] - Exyd[isub] + Exzd[i] - Exzd[isub];
        stagDiffI = Exyd[i + 1] - Exyd[isub + 1] + Exzd[i + 1] - Exzd[isub + 1];
        asgn      = Hyxd[i] * tHyxd[i] - Hyxd[i + 1] * tHyxd[i + 1] + Hybndd[i] - cHyxd[i] * stagDiffR + cHyxd[i + 1] * stagDiffI;
        Hyxd[i + 1]  = Hyxd[i] * tHyxd[i + 1] + Hyxd[i + 1] * tHyxd[i] + Hybndd[i + 1] - cHyxd[i] * stagDiffI - cHyxd[i + 1] * stagDiffR;
        Hyxd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) {
        // Hz_x: ccs, Ex_{y,z}: css
        isub      = i + 2 * ( -nnx );
        stagDiffR = Exyd[isub] - Exyd[i] + Exzd[isub] - Exzd[i];
        stagDiffI = Exyd[isub + 1] - Exyd[i + 1] + Exzd[isub + 1] - Exzd[i + 1];
        asgn      = Hzxd[i] * tHzxd[i] - Hzxd[i + 1] * tHzxd[i + 1] - cHzxd[i] * stagDiffR + cHzxd[i + 1] * stagDiffI;
        Hzxd[i + 1]  = Hzxd[i] * tHzxd[i + 1] + Hzxd[i + 1] * tHzxd[i] - cHzxd[i] * stagDiffI - cHzxd[i + 1] * stagDiffR;
        Hzxd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) {
        // Hx_y: scc, Ey_{x,z}: scs
        isub      = i + 2 * ( -nnx * nny );
        stagDiffR = Eyxd[isub] - Eyxd[i] + Eyzd[isub] - Eyzd[i];
        stagDiffI = Eyxd[isub + 1] - Eyxd[i + 1] + Eyzd[isub + 1] - Eyzd[i + 1];
        asgn      = Hxyd[i] * tHxyd[i] - Hxyd[i + 1] * tHxyd[i + 1] + Hxbndd[i] - cHxyd[i] * stagDiffR + cHxyd[i + 1] * stagDiffI;
        Hxyd[i + 1]  = Hxyd[i] * tHxyd[i + 1] + Hxyd[i + 1] * tHxyd[i] + Hxbndd[i + 1] - cHxyd[i] * stagDiffI - cHxyd[i + 1] * stagDiffR;
        Hxyd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) {
        // Hz_y: ccs, Ey_{x,z}: scs
        isub      = i + 2 * ( -1 );
        stagDiffR = Eyxd[i] - Eyxd[isub] + Eyzd[i] - Eyzd[isub];
        stagDiffI = Eyxd[i + 1] - Eyxd[isub + 1] + Eyzd[i + 1] - Eyzd[isub + 1];
        asgn      = Hzyd[i] * tHzyd[i] - Hzyd[i + 1] * tHzyd[i + 1] - cHzyd[i] * stagDiffR + cHzyd[i + 1] * stagDiffI;
        Hzyd[i + 1]  = Hzyd[i] * tHzyd[i + 1] + Hzyd[i + 1] * tHzyd[i] - cHzyd[i] * stagDiffI - cHzyd[i + 1] * stagDiffR;
        Hzyd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) {
        // Hx_z: scc, Ez_{x,y}: ssc
        isub      = i + 2 * ( -nnx );
        stagDiffR = Ezxd[i] - Ezxd[isub] + Ezyd[i] - Ezyd[isub];
        stagDiffI = Ezxd[i + 1] - Ezxd[isub + 1] + Ezyd[i + 1] - Ezyd[isub + 1];
        asgn      = Hxzd[i] * tHxzd[i] - Hxzd[i + 1] * tHxzd[i + 1] - cHxzd[i] * stagDiffR + cHxzd[i + 1] * stagDiffI;
        Hxzd[i + 1]  = Hxzd[i] * tHxzd[i + 1] + Hxzd[i + 1] * tHxzd[i] - cHxzd[i] * stagDiffI - cHxzd[i + 1] * stagDiffR;
        Hxzd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) {
        // Hy_z: csc, Ez_{x,y}: ssc
        isub      = i + 2 * ( -1 );
        stagDiffR = Ezxd[isub] + Ezyd[isub] - Ezxd[i] - Ezyd[i];
        stagDiffI = Ezxd[isub + 1] + Ezyd[isub + 1] - Ezxd[i + 1] - Ezyd[i + 1];
        asgn      = Hyzd[i] * tHyzd[i] - Hyzd[i + 1] * tHyzd[i + 1] - cHyzd[i] * stagDiffR + cHyzd[i + 1] * stagDiffI;
        Hyzd[i + 1]  = Hyzd[i] * tHyzd[i + 1] + Hyzd[i + 1] * tHyzd[i] - cHyzd[i] * stagDiffI - cHyzd[i + 1] * stagDiffR;
        Hyzd[i]      = asgn;
      }
    }
  }
}
void solar_e_field_ref( const int shape[3], const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u){
  int i,j,k;
  int nnz =shape[2];
  int nny =shape[1];
  int nnx =shape[0];
  unsigned long ln_domain = 2ul*shape[0]*shape[1]*shape[2];
  unsigned long ixmin, ixmax;

  FLOAT_PRECISION * restrict Hyxd = &(u[0ul*ln_domain]);
  FLOAT_PRECISION * restrict Hzxd = &(u[1ul*ln_domain]);
  FLOAT_PRECISION * restrict Hxyd = &(u[2ul*ln_domain]);
  FLOAT_PRECISION * restrict Hzyd = &(u[3ul*ln_domain]);
  FLOAT_PRECISION * restrict Hxzd = &(u[4ul*ln_domain]);
  FLOAT_PRECISION * restrict Hyzd = &(u[5ul*ln_domain]);

  FLOAT_PRECISION * restrict Exzd = &(u[6ul*ln_domain]);
  FLOAT_PRECISION * restrict Eyzd = &(u[7ul*ln_domain]);
  FLOAT_PRECISION * restrict Eyxd = &(u[8ul*ln_domain]);
  FLOAT_PRECISION * restrict Ezxd = &(u[9ul*ln_domain]);
  FLOAT_PRECISION * restrict Exyd = &(u[10ul*ln_domain]);
  FLOAT_PRECISION * restrict Ezyd = &(u[11ul*ln_domain]);

  const FLOAT_PRECISION * restrict cExzd = &(coef[14ul*ln_domain]);
  const FLOAT_PRECISION * restrict cEyzd = &(coef[15ul*ln_domain]);
  const FLOAT_PRECISION * restrict cEyxd = &(coef[16ul*ln_domain]);
  const FLOAT_PRECISION * restrict cEzxd = &(coef[17ul*ln_domain]);
  const FLOAT_PRECISION * restrict cExyd = &(coef[18ul*ln_domain]);
  const FLOAT_PRECISION * restrict cEzyd = &(coef[19ul*ln_domain]);

  const FLOAT_PRECISION * restrict tExzd = &(coef[20ul*ln_domain]);
  const FLOAT_PRECISION * restrict tEyzd = &(coef[21ul*ln_domain]);
  const FLOAT_PRECISION * restrict tEyxd = &(coef[22ul*ln_domain]);
  const FLOAT_PRECISION * restrict tEzxd = &(coef[23ul*ln_domain]);
  const FLOAT_PRECISION * restrict tExyd = &(coef[24ul*ln_domain]);
  const FLOAT_PRECISION * restrict tEzyd = &(coef[25ul*ln_domain]);

  const FLOAT_PRECISION * restrict Exbndd = &(coef[26ul*ln_domain]);
  const FLOAT_PRECISION * restrict Eybndd = &(coef[27ul*ln_domain]);

  unsigned long isub;
  FLOAT_PRECISION stagDiffR, stagDiffI, asgn;

  // Update E-field
  // ---------------------------------------------------------------------------------------
  // -----  Ex_z = Cex_tz * Ex_z + Cex_z * (N(Hz_x + Hz_y) - S(Hz_x + Hz_y) ) --------------
  // ---------------------------------------------------------------------------------------
  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);

      for(i=ixmin; i<ixmax; i+=2) {
        // Ex_z: css, Hz_{x,y}: ccs
        isub      = i + 2 * ( +nnx );
        stagDiffR = Hzxd[isub] - Hzxd[i] + Hzyd[isub] - Hzyd[i];
        stagDiffI = Hzxd[isub + 1] - Hzxd[i + 1] + Hzyd[isub + 1] - Hzyd[i + 1];
        asgn      = Exzd[i] * tExzd[i] - Exzd[i + 1] * tExzd[i + 1] + cExzd[i] * stagDiffR - cExzd[i + 1] * stagDiffI;
        Exzd[i + 1]  = Exzd[i] * tExzd[i + 1] + Exzd[i + 1] * tExzd[i] + cExzd[i] * stagDiffI + cExzd[i + 1] * stagDiffR;
        Exzd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) { 
        // Ey_{x,z}: scs, Hz_{x,y}: ccs
        isub      = i + 2 * ( +1 );
        stagDiffR = Hzxd[i] + Hzyd[i] - Hzxd[isub] - Hzyd[isub];
        stagDiffI = Hzxd[i + 1] + Hzyd[i + 1] - Hzxd[isub + 1] - Hzyd[isub + 1];
        asgn      = Eyzd[i] * tEyzd[i] - Eyzd[i + 1] * tEyzd[i + 1] + cEyzd[i] * stagDiffR - cEyzd[i + 1] * stagDiffI;
        Eyzd[i + 1]  = Eyzd[i] * tEyzd[i + 1] + Eyzd[i + 1] * tEyzd[i] + cEyzd[i] * stagDiffI + cEyzd[i + 1] * stagDiffR;
        Eyzd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) { 
        // Ey_{x,z}: scs, Hx_{y,z}: scc
        isub      = i + 2 * ( +nnx * nny );
        stagDiffR = Hxyd[isub] - Hxyd[i] + Hxzd[isub] - Hxzd[i];
        stagDiffI = Hxyd[isub + 1] - Hxyd[i + 1] + Hxzd[isub + 1] - Hxzd[i + 1];
        asgn      = Eyxd[i] * tEyxd[i] - Eyxd[i + 1] * tEyxd[i + 1] + Eybndd[i] + cEyxd[i] * stagDiffR - cEyxd[i + 1] * stagDiffI;
        Eyxd[i + 1]  = Eyxd[i] * tEyxd[i + 1] + Eyxd[i + 1] * tEyxd[i] + Eybndd[i + 1] + cEyxd[i] * stagDiffI + cEyxd[i + 1] * stagDiffR;
        Eyxd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) { 
        // Ez_x: ssc, Hx_{y,z}: scc
        isub      = i + 2 * ( +nnx );
        stagDiffR = Hxyd[i] + Hxzd[i] - Hxyd[isub] - Hxzd[isub];
        stagDiffI = Hxyd[i + 1] + Hxzd[i + 1] - Hxyd[isub + 1] - Hxzd[isub + 1];
        asgn      = Ezxd[i] * tEzxd[i] - Ezxd[i + 1] * tEzxd[i + 1] + cEzxd[i] * stagDiffR - cEzxd[i + 1] * stagDiffI;
        Ezxd[i + 1]  = Ezxd[i] * tEzxd[i + 1] + Ezxd[i + 1] * tEzxd[i] + cEzxd[i] * stagDiffI + cEzxd[i + 1] * stagDiffR;
        Ezxd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) { 
        // Ex_y: css, Hy_{x,z}: csc
        isub      = i + 2 * ( +nnx * nny );
        stagDiffR = Hyxd[i] - Hyxd[isub] + Hyzd[i] - Hyzd[isub];
        stagDiffI = Hyxd[i + 1] - Hyxd[isub + 1] + Hyzd[i + 1] - Hyzd[isub + 1];
        asgn      = Exyd[i] * tExyd[i] - Exyd[i + 1] * tExyd[i + 1] + Exbndd[i] + cExyd[i] * stagDiffR - cExyd[i + 1] * stagDiffI;
        Exyd[i + 1]  = Exyd[i] * tExyd[i + 1] + Exyd[i + 1] * tExyd[i] + Exbndd[i + 1] + cExyd[i] * stagDiffI + cExyd[i + 1] * stagDiffR;
        Exyd[i]      = asgn;
      }
    }
  }

  for(k=2; k<nnz-2; k++) {
    for(j=2; j<nny-2; j++) {

      ixmin  = 2 * ( ( k * nny + j ) * nnx + 2);
      ixmax  = 2 * ( ( k * nny + j ) * nnx + nnx-2);
      for(i=ixmin; i<ixmax; i+=2) { 
        // Ez_y: ssc, Hy_{x,z}: csc
        isub      = i + 2 * ( +1 );
        stagDiffR = Hyxd[isub] - Hyxd[i] + Hyzd[isub] - Hyzd[i];
        stagDiffI = Hyxd[isub + 1] - Hyxd[i + 1] + Hyzd[isub + 1] - Hyzd[i + 1];
        asgn      = Ezyd[i] * tEzyd[i] - Ezyd[i + 1] * tEzyd[i + 1] + cEzyd[i] * stagDiffR - cEzyd[i + 1] * stagDiffI;
        Ezyd[i + 1]  = Ezyd[i] * tEzyd[i + 1] + Ezyd[i + 1] * tEzyd[i] + cEzyd[i] * stagDiffI + cEzyd[i + 1] * stagDiffR;
        Ezyd[i]      = asgn;
      }      
    }
  } 

}
void solar_kernel( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2) {

  if(u==0) u= (FLOAT_PRECISION * restrict) v; // in odd iterations the solution domain will be passed in v array (temporary patch)
  solar_h_field_ref(shape, coef, u);
  solar_e_field_ref(shape, coef, u);
}




void compare_results_std(FLOAT_PRECISION *restrict u, FLOAT_PRECISION *restrict target_domain, int alignment, int nx, int ny, int nz, int NHALO){
  int nnx=nx+2*NHALO, nny=ny+2*NHALO, nnz=nz+2*NHALO;
  unsigned long n_domain = nnx*nny*nnz;
  int i, j, k, male;

  double ref_l1 = 0.0;
  FLOAT_PRECISION diff_l1=0.0, abs_diff, max_error=0.0;
  FLOAT_PRECISION * restrict snapshot_error;
  male = posix_memalign((void **)&(snapshot_error), alignment, sizeof(FLOAT_PRECISION)*n_domain); check_merr(male);
  for (k=0; k<nz; k++) {
    for (j=0; j<ny; j++) {
      for(i=0; i<nx; i++) {
        abs_diff = fabs(U(i+NHALO, j+NHALO, k+NHALO) - TARGET_DOMAIN(i, j, k));
        SNAPSHOT_ERROR(i, j, k) = U(i+NHALO, j+NHALO, k+NHALO) - TARGET_DOMAIN(i, j, k);
        if (abs_diff > max_error) max_error = abs_diff;
        diff_l1 += abs_diff;
      }
    }
  }
  if ( (diff_l1 > 0.0) || (diff_l1*0 != 0) || (diff_l1 != diff_l1) ){
    printf("Max snapshot abs. err.:%e  L1 norm:%e\n", max_error, diff_l1);
    fprintf(stderr,"BROKEN KERNEL\n");

    print_3Darray("error_snapshot", snapshot_error, nnx, nny, nnz, NHALO);
    print_3Darray("reference_snapshot", u , nnx, nny, nnz, NHALO);
    print_3Darray("target_snapshot"   , target_domain, nx, ny, nz, 0);
    exit(1);
  } else {
    printf("eMax:%.3e|eL1:%.3e", max_error, diff_l1);
    printf("-PASSED\n");

    for (k=0; k<n_domain; k++) ref_l1 += fabs(u[k]);
    if(ref_l1 < 1e-6) printf("##WARNING: The L1 norm of the solution domain is too small (%e). This verification is not sufficient to discover errors\n", ref_l1);

  }

  free(snapshot_error);
}
void compare_results_solar(FLOAT_PRECISION *restrict u, Parameters p){
  FLOAT_PRECISION *restrict target_domain = p.U1;

  int nnx = p.ldomain_shape[0];
  int nny = p.ldomain_shape[1];
  int nnz = p.ldomain_shape[2];
  unsigned long idx, n_domain = nnx*nny*nnz;
  int f, i, j, k, male;
  double ref_l1 = 0.0;
  FLOAT_PRECISION diff_l1=0.0, abs_diff, max_error=0.0;
  FLOAT_PRECISION * restrict snapshot_error;
  male = posix_memalign((void **)&(snapshot_error), p.alignment, sizeof(FLOAT_PRECISION)*n_domain*24lu); check_merr(male);
  for(f=0; f<12;f++){
    for (k=0; k<nnz; k++) {
      for (j=0; j<nny; j++) {
        for(i=0; i<nnx; i++) {
          idx = 2*( n_domain*f + ( k   *nny +  j   ) * nnx +   i);
          abs_diff = fabs( u[idx] - target_domain[idx]);
          snapshot_error[idx] = abs_diff;
          if (abs_diff > max_error) max_error = abs_diff;
          diff_l1 += abs_diff;
        }
      }
    }
  }
  if ( (diff_l1 > 0.0) || (diff_l1*0 != 0) || (diff_l1 != diff_l1) ){
    printf("Max snapshot abs. err.:%e  L1 norm:%e\n", max_error, diff_l1);
    fprintf(stderr,"BROKEN KERNEL\n");

    print_3Darray_solar("error_snapshot", snapshot_error, nnx, nny, nnz, 0);
    print_3Darray_solar("reference_snapshot", u , nnx, nny, nnz, 0);
    print_3Darray_solar("target_snapshot"   , target_domain, nnx, nny, nnz, 0);
    exit(1);
  } else {
    printf("eMax:%.3e|eL1:%.3e", max_error, diff_l1);
    printf("-PASSED\n");

    for (k=0; k<n_domain; k++) ref_l1 += fabs(u[k]);
    if(ref_l1 < 1e-6) printf("##WARNING: The L1 norm of the solution domain is too small (%e). This verification is not sufficient to discover errors\n", ref_l1);

  }

  free(snapshot_error); 
}
void compare_results(FLOAT_PRECISION *restrict u, FLOAT_PRECISION *restrict target_domain, int alignment, int nx, int ny, int nz, int NHALO, Parameters p){
  if(p.stencil.type == REGULAR){
     compare_results_std(u, target_domain, alignment, nx, ny, nz, NHALO);
     return;
  }
  compare_results_solar(u, p);
}


void verification_printing(Parameters vp){
  char *coeff_type, *precision, *concat;
  if(vp.mpi_rank==0) {
    switch (vp.stencil.coeff){
    case CONSTANT_COEFFICIENT:
      coeff_type = "const    ";
      break;
    case VARIABLE_COEFFICIENT:
      coeff_type = "var      ";
      break;
    case VARIABLE_COEFFICIENT_AXSYM:
      coeff_type = "var_axsym";
      break;
    case VARIABLE_COEFFICIENT_NOSYM:
      coeff_type = "var_nosym";
      break;
    case SOLAR_COEFFICIENT:
      coeff_type = "Solar kernel";
      break;
    }
    precision = ((sizeof(FLOAT_PRECISION)==4)?"SP":"DP");
    concat = ((vp.halo_concat==0)?"no-concat":"   concat");
    printf("#ts:%s stencil:%s|R:%d|T:%d|%s nt:%03d thrd:%d prec.:%s %s"
        " dom:(%d,%03d,%03d) top:(%d,%d,%d) ",
        TSList[vp.target_ts].name, vp.stencil.name,
        vp.stencil.r,
        vp.stencil.time_order, coeff_type,
        vp.nt, vp.num_threads, precision, concat,
        vp.lstencil_shape[0], vp.lstencil_shape[1], vp.lstencil_shape[2],
        vp.t.shape[0], vp.t.shape[1], vp.t.shape[2]);
    if(vp.target_ts == 2)
      printf("TB:%d wf:%d ",vp.t_dim, vp.wavefront);
    if(vp.target_ts == 2){
      printf("|thrd_group|:%d ",vp.stencil_ctx.thread_group_size);
      printf("num-wf:%d ", vp.stencil_ctx.num_wf);
    }
    if(vp.verbose ==1) print_param(vp);
  }
}


FLOAT_PRECISION *restrict aggregate_subdomains(Parameters vp){
  // aggregate the domains if there are more than 1 MPI ranks
  int i, j, k, ind, sind, male;
  FLOAT_PRECISION * restrict aggr_domain;

  if(vp.mpi_rank==0) {
    male = posix_memalign((void **)&(aggr_domain), vp.alignment,
        sizeof(FLOAT_PRECISION)*vp.n_stencils); check_merr(male);
  }
  if(vp.mpi_size > 1) {
    // aggregate the domains
    aggregate_MPI_subdomains(vp, aggr_domain);
  } else { // copy the domain without the halo points
    for(i=0; i<vp.stencil_shape[2]; i++) {
      for(j=0; j<vp.stencil_shape[1]; j++) {
        for(k=0; k<vp.stencil_shape[0]; k++) {
          ind = ((i+vp.stencil.r)*vp.ldomain_shape[1] +
              (j+vp.stencil.r)) * vp.ldomain_shape[0] + (k+vp.stencil.r);
          sind = (i*vp.stencil_shape[1] + j) * vp.stencil_shape[0] + k;
          aggr_domain[sind] = vp.U1[ind];
        }
      }
    }
  }
  return aggr_domain;
}

void aggregate_MPI_subdomains(Parameters vp, FLOAT_PRECISION * restrict aggr_domain){
  int i;
  // create MPI type to send the data without the halo points
  MPI_Datatype stencil_type, aggr_stencil_type;
  MPI_Status status;
  MPI_Request wait_req;
  MPI_Status wait_stat;
  int stencil_starts[3];
  int subdomain_info[6];
  stencil_starts[0] = stencil_starts[1] = stencil_starts[2] = vp.stencil.r;
  ierr= MPI_Type_create_subarray(3, vp.ldomain_shape, vp.lstencil_shape,
     stencil_starts, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &stencil_type); CHKERR(ierr);
  ierr = MPI_Type_commit(&stencil_type); CHKERR(ierr);


  for(i=0; i<vp.mpi_size; i++) {
    // default initialization for the other ranks
    subdomain_info[0] = vp.lstencil_shape[0];
    subdomain_info[1] = vp.lstencil_shape[1];
    subdomain_info[2] = vp.lstencil_shape[2];
    subdomain_info[3] = vp.gb[0];
    subdomain_info[4] = vp.gb[1];
    subdomain_info[5] = vp.gb[2];

    // send the coordinates information to rank 0
    if(i > 0){
      if (vp.mpi_rank == i) {
        ierr = MPI_Send(subdomain_info, 6, MPI_INT, 0, 0,
            vp.t.cart_comm); CHKERR(ierr);
      }
      if (vp.mpi_rank==0) {
        ierr = MPI_Recv(subdomain_info, 6, MPI_INT, i, MPI_ANY_TAG,
            vp.t.cart_comm, &status); CHKERR(ierr);
      }
    }

    // create data type to put the results in place at rank 0
    ierr = MPI_Type_create_subarray(3, vp.stencil_shape, subdomain_info,
        &(subdomain_info[3]), MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION,
        &aggr_stencil_type); CHKERR(ierr);
    ierr = MPI_Type_commit(&aggr_stencil_type); CHKERR(ierr);

    if (vp.mpi_rank==i) {
      ierr = MPI_Isend(vp.U1, 1, stencil_type, 0, 0, vp.t.cart_comm,
          &wait_req); CHKERR(ierr);
    }
    if (vp.mpi_rank==0) {
      ierr = MPI_Recv(aggr_domain, 1, aggr_stencil_type, i, MPI_ANY_TAG,
          vp.t.cart_comm, &status); CHKERR(ierr);
    }
    if (vp.mpi_rank==i) {
      ierr = MPI_Waitall(1, &wait_req, &wait_stat); CHKERR(ierr);
    }

    MPI_Barrier(vp.t.cart_comm);
    ierr = MPI_Type_free(&aggr_stencil_type); CHKERR(ierr);
  }
  ierr = MPI_Type_free(&stencil_type); CHKERR(ierr);
}
