#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include "verification.h"

void verify(Parameters *p){

//printf("[%d]  Here\n", p->mpi_rank);
//printf("[%d]1  Here\n", p->mpi_rank);

//printf("[%d]2  Here\n", p->mpi_rank);
  verification_printing(*p);
  mpi_halo_init(p);

//printf("[%d]3  Here\n", p->mpi_rank);
  // allocate and initialize the required arrays
  arrays_allocate(p);
  init_coeff(p);
//printf("[%d]4  Here\n", p->mpi_rank);
  domain_data_fill(p);

//printf("[%d]5.5  Here\n", p->mpi_rank);
  // compute data in parallel using the time stepper to be tested
//printf("[%d] target    nt:%d\n", p->mpi_rank, p->nt);
  TSList[p->target_ts].func(p);
//printf("[%d]6  Here\n", p->mpi_rank);

  // aggregate all subdomains into rank zero to compare with the serial results
  FLOAT_PRECISION * restrict aggr_domain = aggregate_subdomains(*p);
//printf("[%d]7  Here\n", p->mpi_rank);

  // compute data serially using reference time stepper
  //    compare and report results
  if(p->mpi_rank==0) {
    verify_serial_generic(aggr_domain, *p);
    free(aggr_domain);
  }
//printf("[%d]8  Here\n", p->mpi_rank);

  mpi_halo_finalize(p);
  arrays_free(p);
//printf("[%d]9  Here\n", p->mpi_rank);
}
void verify_serial_generic(FLOAT_PRECISION * target_domain, Parameters p) {
//  int nt = p.nt;
//  int nt2;
  int male;
  // function pointer to select the reference stencil operator
  void (*std_kernel)( const int [],
      const FLOAT_PRECISION * restrict, FLOAT_PRECISION * restrict,
      const FLOAT_PRECISION * restrict, const FLOAT_PRECISION * restrict);

  FLOAT_PRECISION * restrict u, * restrict v, * restrict roc2, * restrict coef;
  int it, i, j , k, ax;
  int nnx = p.stencil_shape[0]+2*NHALO;
  int nny = p.stencil_shape[1]+2*NHALO;
  int nnz = p.stencil_shape[2]+2*NHALO;
  unsigned long n_domain = nnx*nny*nnz;

  male = posix_memalign((void **)&(u),    p.alignment, sizeof(FLOAT_PRECISION)*n_domain); check_merr(male);
  male = posix_memalign((void **)&(v),    p.alignment, sizeof(FLOAT_PRECISION)*n_domain); check_merr(male);

  if(KernelList[p.target_kernel].time_order == 2)
    male = posix_memalign((void **)&(roc2), p.alignment, sizeof(FLOAT_PRECISION)*n_domain); check_merr(male);


  // allocate the the coefficients according to the stencil type
  // initialize the coefficients
  // select the desired stencil operator function
  unsigned long coef_size;
  switch(KernelList[p.target_kernel].stencil_coeff){
  case CONSTANT_COEFFICIENT:
    coef_size = 11;
    male = posix_memalign((void **)&(coef), p.alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male);
    for(i=0;i<NHALO+1;i++) coef[i] = p.g_coef[i];

    // select the stencil type
    switch(KernelList[p.target_kernel].stencil_radius){
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
    coef_size = n_domain*(1 + KernelList[p.target_kernel].stencil_radius);
    male = posix_memalign((void **)&(coef), p.alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male);
    for(k=0; k <= KernelList[p.target_kernel].stencil_radius; k++){
      for(i=0; i<n_domain; i++){
        coef[i + k*n_domain] = p.g_coef[k];
      }
    }
    // select the stencil type
    switch(KernelList[p.target_kernel].stencil_radius){
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
    coef_size = n_domain*(1 + 3*KernelList[p.target_kernel].stencil_radius);
    male = posix_memalign((void **)&(coef), p.alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male);

    // central point coeff
    for(i=0; i<n_domain; i++){
      coef[i] = p.g_coef[0];
    }
    for(k=0; k < KernelList[p.target_kernel].stencil_radius; k++){
      for(ax=0; ax<3; ax++){
        for(i=0; i<n_domain; i++){
          coef[i + n_domain + 3*k*n_domain + ax*n_domain] = p.g_coef[k+1];
        }
      }
    }

    // select the stencil type
    switch(KernelList[p.target_kernel].stencil_radius){
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
      coef_size = n_domain*(1 + 6*KernelList[p.target_kernel].stencil_radius);
      male = posix_memalign((void **)&(coef), p.alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male);

      // central point coeff
      for(i=0; i<n_domain; i++){
        coef[i] = p.g_coef[0];
      }
      for(k=0; k < KernelList[p.target_kernel].stencil_radius; k++){
        for(ax=0; ax<3; ax++){
          for(i=0; i<n_domain; i++){
            coef[i + n_domain + 6*k*n_domain +  2*ax   *n_domain] = p.g_coef[k+1];
            coef[i + n_domain + 6*k*n_domain + (2*ax+1)*n_domain] = p.g_coef[k+1];
          }
        }
      }

      // select the stencil type
      switch(KernelList[p.target_kernel].stencil_radius){
      case 1:
        std_kernel = &std_kernel_2space_1time_var_nosym;
        break;
      default:
        printf("ERROR: unknown type of stencil\n");
        exit(1);
        break;
      }
      break;


  default:
    printf("ERROR: unknown type of stencil\n");
    exit(1);
    break;
  }


  //-- Initialize u,v,roc2
  // fill the array points according to their location in space and pad the boundary with zeroes
  for(i=0; i<n_domain;i++){
    u[i] = 0.0;
    v[i] = 0.0;
    if(KernelList[p.target_kernel].time_order == 2)
      roc2[i]= 0.0;
  }
  FLOAT_PRECISION r;
  for(k=NHALO; k<p.stencil_shape[2]+NHALO; k++){
    for(j=NHALO; j<p.stencil_shape[1]+NHALO; j++){
      for(i=NHALO; i<p.stencil_shape[0]+NHALO; i++){
        r = 1.0/3 * (1.0*i/p.stencil_shape[0] + 1.0*j/p.stencil_shape[1]  + 1.0*k/p.stencil_shape[2]);
        U(i, j, k)    = r*1.845703;
        V(i, j, k)    = r*1.845703;

        if(KernelList[p.target_kernel].time_order == 2)
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
//  u[(NHALO+1)*(nnx * nny + nny + 1)] += 100.1;

  // compare results
  compare_results(u, target_domain, p.alignment, p.stencil_shape[0], p.stencil_shape[1], p.stencil_shape[2]);

  free(u);
  free(v);

  if(KernelList[p.target_kernel].time_order == 2)
    free(roc2);
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

void compare_results(FLOAT_PRECISION *restrict u, FLOAT_PRECISION *restrict target_domain, int alignment, int nx, int ny, int nz){
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


void verification_printing(Parameters vp){
  char *coeff_type, *precision, *concat;
  if(vp.mpi_rank==0) {
    switch (KernelList[vp.target_kernel].stencil_coeff){
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
    }
    precision = ((sizeof(FLOAT_PRECISION)==4)?"SP":"DP");
    concat = ((vp.halo_concat==0)?"no-concat":"   concat");
    printf("#ts:%s stencil:%s|R:%d|T:%d|%s nt:%03d thrd:%d prec.:%s %s"
        " dom:(%d,%03d,%03d) top:(%d,%d,%d) ",
        TSList[vp.target_ts].name, KernelList[vp.target_kernel].name,
        KernelList[vp.target_kernel].stencil_radius,
        KernelList[vp.target_kernel].time_order, coeff_type,
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
          ind = ((i+NHALO)*vp.ldomain_shape[1] +
              (j+NHALO)) * vp.ldomain_shape[0] + (k+NHALO);
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
  stencil_starts[0] = stencil_starts[1] = stencil_starts[2] = NHALO;
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
