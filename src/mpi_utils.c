#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include "data_structures.h"

extern void copyv(int *a, int * b, int n);

#define SRC_BUF(i,j,k)         (src_buf[((k)*(src_size[1])+(j))*(src_size[0])+(i)])
#define DST_BUF(i,j,k)         (dst_buf[((k)*(dst_size[1])+(j))*(dst_size[0])+(i)])

void sub_array_copy(const FLOAT_PRECISION * restrict src_buf, FLOAT_PRECISION * restrict dst_buf, int *src_size, int *dst_size, int *cpy_size, int *src_offs, int *dst_offs){
    int *ds = dst_size;
    int i,j,k;
#pragma omp parallel
  {
#pragma omp for private(k,j,i) schedule(static)
    for(i=0; i<cpy_size[2]; i++){
      for(j=0; j<cpy_size[1]; j++){
        for(k=0; k<cpy_size[0]; k++){
          DST_BUF(k+dst_offs[0],j+dst_offs[1],i+dst_offs[2]) = SRC_BUF(k+src_offs[0], j+src_offs[1], i+src_offs[2]);
        }
      }
    }
  }
}

void sub_array_copy_tg(const FLOAT_PRECISION * restrict src_buf, FLOAT_PRECISION * restrict dst_buf, int *src_size, int *dst_size, int *cpy_size, int *src_offs, int *dst_offs, int tgs){ 
    int *ds = dst_size;
    int i,j,k;
#pragma omp parallel num_threads(tgs)
  {
#pragma omp for private(k,j,i) schedule(static) 
    for(i=0; i<cpy_size[2]; i++){ 
      for(j=0; j<cpy_size[1]; j++){ 
        for(k=0; k<cpy_size[0]; k++){
          DST_BUF(k+dst_offs[0],j+dst_offs[1],i+dst_offs[2]) = SRC_BUF(k+src_offs[0], j+src_offs[1], i+src_offs[2]);
        }
      }
    }
  }
}

void mpi_topology_init(Parameters *p) {
  int old_rank;
  old_rank = p->mpi_rank;

//  if(p->is_diamond_ts == 1){ // diamond methods
//    p->t.is_periodic[2] = 1; // diamonds across the z-direction perform periodic communication
//  }
  if(p->target_ts == 2){ // intra-diamond methods
    p->t.is_periodic[1] = 1; // diamonds across the y-direction perform periodic communication
  }

  ierr = MPI_Cart_create(MPI_COMM_WORLD, 3, p->t.shape, p->t.is_periodic, 1, &(p->t.cart_comm)); CHKERR(ierr);
  // MPI_Errhandler_set(p->t.cart_comm, MPI_ERRORS_RETURN);
  MPI_Comm_rank (p->t.cart_comm, &(p->mpi_rank));
  ierr = MPI_Cart_coords(p->t.cart_comm, p->mpi_rank, 3, p->t.rank_coords); CHKERR(ierr);
  ierr = MPI_Cart_shift(p->t.cart_comm, 0, 1, &(p->t.left), &(p->t.right)); CHKERR(ierr); // in X direction
  ierr = MPI_Cart_shift(p->t.cart_comm, 1, 1, &(p->t.down), &(p->t.up)); CHKERR(ierr); // in Y direction
  ierr = MPI_Cart_shift(p->t.cart_comm, 2, 1, &(p->t.back), &(p->t.front)); CHKERR(ierr); // in Z direction

  if(p->debug ==1){
    MPI_Barrier(MPI_COMM_WORLD);
    if(p->mpi_rank == 0) {
      printf("\n******************************************************\n"); fflush(stdout);
      printf("DEBUG topology initialization information BEGIN\n"); fflush(stdout);
      printf("******************************************************\n"); fflush(stdout); sleep(1);
    }
    int j;
    for(j=0; j<p->mpi_size; j++){
      if(j == old_rank){
        printf("Rank %03d -> %03d\n", old_rank, p->mpi_rank); fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    sleep(1);
    MPI_Barrier(MPI_COMM_WORLD);

    for(j=0; j<p->mpi_size; j++){
      if(j == p->mpi_rank){
        printf("Rank %03d topology (npx, npy, npz):(%02d,%02d,%02d) | left:%03d right:%03d | down:%03d up:%03d | back:%03d front:%03d\n",
            p->mpi_rank, p->t.rank_coords[0], p->t.rank_coords[1], p->t.rank_coords[2], p->t.left, p->t.right, p->t.down, p->t.up, p->t.back, p->t.front); fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    sleep(1);
    if(p->mpi_rank == 0) {
      printf("******************************************************\n"); fflush(stdout);
      printf("DEBUG topology initialization information END\n"); fflush(stdout);
      printf("******************************************************\n"); fflush(stdout);
    }
  }
}


void standard_mpi_halo_init(Parameters *p){
  if (p->t.shape[0] > 1) {
    p->h[0].shape[0]= p->stencil.r;
    p->h[0].shape[1]= p->ldomain_shape[1]-2*p->stencil.r;
    p->h[0].shape[2]= p->ldomain_shape[2]-2*p->stencil.r;

    // Y and Z direction have fixed Halo beginning across X
    p->h[0].recv_b[0]= 0;     p->h[0].recv_e[0]= p->lstencil_shape[0]+p->stencil.r;
    p->h[0].recv_b[1]= p->stencil.r; p->h[0].recv_e[1]= p->stencil.r;
    p->h[0].recv_b[2]= p->stencil.r; p->h[0].recv_e[2]= p->stencil.r;

    p->h[0].send_b[0]= p->stencil.r; p->h[0].send_e[0]= p->lstencil_shape[0];
    p->h[0].send_b[1]= p->stencil.r; p->h[0].send_e[1]= p->stencil.r;
    p->h[0].send_b[2]= p->stencil.r; p->h[0].send_e[2]= p->stencil.r;

    ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->h[0].shape, p->h[0].recv_b, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->h[0].recv_hb)); CHKERR(ierr);
    ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->h[0].shape, p->h[0].recv_e  , MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->h[0].recv_he)); CHKERR(ierr);
    ierr = MPI_Type_commit(&(p->h[0].recv_hb)); CHKERR(ierr);
    ierr = MPI_Type_commit(&(p->h[0].recv_he)); CHKERR(ierr);

    ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->h[0].shape, p->h[0].send_b, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->h[0].send_hb)); CHKERR(ierr);
    ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->h[0].shape, p->h[0].send_e  , MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->h[0].send_he)); CHKERR(ierr);
    ierr = MPI_Type_commit(&(p->h[0].send_hb)); CHKERR(ierr);
    ierr = MPI_Type_commit(&(p->h[0].send_he)); CHKERR(ierr);

    p->h[0].size = p->h[0].shape[0] * p->h[0].shape[1] * p->h[0].shape[2];
  }

  if (p->t.shape[1] > 1) {
    p->h[1].shape[0]= p->lstencil_shape[0];
    p->h[1].shape[1]= p->stencil.r;
    p->h[1].shape[2]= p->ldomain_shape[2]-2*p->stencil.r;

    // X and Z direction have fixed Halo beginning across Y
    p->h[1].recv_b[0]= p->stencil.r; p->h[1].recv_e[0]= p->stencil.r;
    p->h[1].recv_b[1]= 0;     p->h[1].recv_e[1]= p->ldomain_shape[1]-p->stencil.r;
    p->h[1].recv_b[2]= p->stencil.r; p->h[1].recv_e[2]= p->stencil.r;

    p->h[1].send_b[0]= p->stencil.r; p->h[1].send_e[0]= p->stencil.r;
    p->h[1].send_b[1]= p->stencil.r; p->h[1].send_e[1]= p->ldomain_shape[1]-2*p->stencil.r;
    p->h[1].send_b[2]= p->stencil.r; p->h[1].send_e[2]= p->stencil.r;

    ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->h[1].shape, p->h[1].recv_b, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->h[1].recv_hb)); CHKERR(ierr);
    ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->h[1].shape, p->h[1].recv_e  , MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->h[1].recv_he)); CHKERR(ierr);
    ierr = MPI_Type_commit(&(p->h[1].recv_hb)); CHKERR(ierr);
    ierr = MPI_Type_commit(&(p->h[1].recv_he)); CHKERR(ierr);

    ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->h[1].shape, p->h[1].send_b, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->h[1].send_hb)); CHKERR(ierr);
    ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->h[1].shape, p->h[1].send_e  , MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->h[1].send_he)); CHKERR(ierr);
    ierr = MPI_Type_commit(&(p->h[1].send_hb)); CHKERR(ierr);
    ierr = MPI_Type_commit(&(p->h[1].send_he)); CHKERR(ierr);

    p->h[1].size = p->h[1].shape[0] * p->h[1].shape[1] * p->h[1].shape[2];
  }

  int z_halo_size, xy_plain;
  if (p->t.shape[2] > 1) {
    if (p->h[2].is_contiguous==1){
      p->h[2].shape[0]= p->ldomain_shape[0];
      p->h[2].shape[1]= p->ldomain_shape[1];
    } else{
      p->h[2].shape[0]= p->lstencil_shape[0];
      p->h[2].shape[1]= p->lstencil_shape[1];
    }

    p->h[2].shape[2]= p->stencil.r;

    p->h[2].recv_b[0]= p->stencil.r;
    p->h[2].recv_b[1]= p->stencil.r;

    xy_plain = p->ldomain_shape[0] * p->ldomain_shape[1];
    p->h[2].recv_b[2] = 0;
    p->h[2].send_b[2] = xy_plain * p->stencil.r;
    p->h[2].send_e[2] = xy_plain * (p->ldomain_shape[2]-2*p->stencil.r);
    p->h[2].recv_e[2] = xy_plain * (p->ldomain_shape[2]-p->stencil.r);

    if(p->h[2].is_contiguous == 0){
      ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->h[2].shape, p->h[2].recv_b, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->h[2].halo)); CHKERR(ierr);
    } else{ // contiguous type is required
      z_halo_size = p->h[2].shape[0] * p->h[2].shape[1] * p->h[2].shape[2];
      ierr = MPI_Type_contiguous(z_halo_size, MPI_FLOAT_PRECISION, &(p->h[2].halo)); CHKERR(ierr);
    }
    ierr = MPI_Type_commit(&(p->h[2].halo)); CHKERR(ierr);

    p->h[2].size = p->h[2].shape[0] * p->h[2].shape[1] * p->h[2].shape[2];
  }
}
void intra_diamond_mpi_halo_init(Parameters *p) {
  if(p->t.shape[1] == 1) return;

  // Halo data types for the buffer u
  p->hu[1].shape[0]= p->lstencil_shape[0];
  p->hu[1].shape[1]= p->stencil.r * (p->t_dim+1);
  p->hu[1].shape[2]= p->lstencil_shape[2];
  p->hu[1].size = p->hu[1].shape[0] * p->hu[1].shape[1] * p->hu[1].shape[2];

  // X and Z directions have fixed Halo beginning across Y
  p->hu[1].recv_b[0]= p->stencil.r;
  p->hu[1].recv_b[1]= p->stencil.r; // NOTE: shifted by halo size
  p->hu[1].recv_b[2]= p->stencil.r;

  p->hu[1].recv_e[0]= p->stencil.r;
  p->hu[1].recv_e[1]= p->ldomain_shape[1] - p->stencil.r*(p->t_dim +2);
  p->hu[1].recv_e[2]= p->stencil.r;

  p->hu[1].send_b[0]= p->stencil.r;
  p->hu[1].send_b[1]= p->stencil.r;
  p->hu[1].send_b[2]= p->stencil.r;

  p->hu[1].send_e[0]= p->stencil.r;
  p->hu[1].send_e[1]= p->ldomain_shape[1] - p->stencil.r*(p->t_dim +2);
  p->hu[1].send_e[2]= p->stencil.r;

  if(p->debug ==1){
    MPI_Barrier(MPI_COMM_WORLD);
    if(p->mpi_rank == 0) {
      printf("\n******************************************************\n"); fflush(stdout);
      printf("DEBUG u halo information BEGIN\n"); fflush(stdout);
      printf("******************************************************\n"); fflush(stdout); sleep(1);
    }
    int j;
    for(j=0; j<p->mpi_size; j++){
      if(j == p->mpi_rank){

        printf("[%02d]:top(%02d,%02d,%02d)\n", p->mpi_rank, p->t.rank_coords[0], p->t.rank_coords[1], p->t.rank_coords[2]); fflush(stdout);
        printf("  halo shape:(%03d,%03d,%03d)\n", p->hu[1].shape[0], p->hu[1].shape[1], p->hu[1].shape[2]); fflush(stdout);
        printf("  Recv begin:(%03d,%03d,%03d)\n", p->hu[1].recv_b[0], p->hu[1].recv_b[1], p->hu[1].recv_b[2]); fflush(stdout);
        printf("  Recv end  :(%03d,%03d,%03d)\n", p->hu[1].recv_e[0], p->hu[1].recv_e[1], p->hu[1].recv_e[2]); fflush(stdout);
        printf("  Send begin:(%03d,%03d,%03d)\n", p->hu[1].send_b[0], p->hu[1].send_b[1], p->hu[1].send_b[2]); fflush(stdout);
        printf("  Send end  :(%03d,%03d,%03d)\n", p->hu[1].send_e[0], p->hu[1].send_e[1], p->hu[1].send_e[2]); fflush(stdout);
        printf("\n"); fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if(p->mpi_rank == 0) {
      printf("******************************************************\n"); fflush(stdout);
      printf("DEBUG u halo information END\n"); fflush(stdout);
      printf("******************************************************\n\n"); fflush(stdout);
    }
  }

  ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->hu[1].shape, p->hu[1].recv_b, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->hu[1].recv_hb)); CHKERR(ierr);
  ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->hu[1].shape, p->hu[1].recv_e  , MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->hu[1].recv_he)); CHKERR(ierr);
  ierr = MPI_Type_commit(&(p->hu[1].recv_hb)); CHKERR(ierr);
  ierr = MPI_Type_commit(&(p->hu[1].recv_he)); CHKERR(ierr);

  ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->hu[1].shape, p->hu[1].send_b, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->hu[1].send_hb)); CHKERR(ierr);
  ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->hu[1].shape, p->hu[1].send_e  , MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->hu[1].send_he)); CHKERR(ierr);
  ierr = MPI_Type_commit(&(p->hu[1].send_hb)); CHKERR(ierr);
  ierr = MPI_Type_commit(&(p->hu[1].send_he)); CHKERR(ierr);


  // Halo data types for the buffer v
  copyv(p->hu[1].shape , p->hv[1].shape, 3);

  copyv(p->hu[1].recv_b, p->hv[1].recv_b, 3);
  p->hv[1].recv_b[1]= 0;

  copyv(p->hu[1].recv_e, p->hv[1].recv_e, 3);
  p->hv[1].recv_e[1]= p->ldomain_shape[1] - p->stencil.r*(p->t_dim + 1);

  copyv(p->hu[1].send_b, p->hv[1].send_b, 3);
  p->hv[1].send_b[1]= 2*p->stencil.r;

  copyv(p->hu[1].send_e, p->hv[1].send_e, 3);
  p->hv[1].send_e[1]= p->ldomain_shape[1] - p->stencil.r*(p->t_dim +3);

  p->hv[1].size = p->hu[1].size;

  if(p->debug ==1){
    MPI_Barrier(MPI_COMM_WORLD);
    if(p->mpi_rank == 0) {
      printf("\n******************************************************\n"); fflush(stdout);
      printf("DEBUG v halo information BEGIN\n"); fflush(stdout);
      printf("******************************************************\n"); fflush(stdout); sleep(1);
    }
    int j;
    for(j=0; j<p->mpi_size; j++){
      if(j == p->mpi_rank){

        printf("[%02d]:top(%02d,%02d,%02d)\n", p->mpi_rank, p->t.rank_coords[0], p->t.rank_coords[1], p->t.rank_coords[2]); fflush(stdout);
        printf("  halo shape:(%03d,%03d,%03d)\n", p->hv[1].shape[0], p->hv[1].shape[1], p->hv[1].shape[2]); fflush(stdout);
        printf("  Recv begin:(%03d,%03d,%03d)\n", p->hv[1].recv_b[0], p->hv[1].recv_b[1], p->hv[1].recv_b[2]); fflush(stdout);
        printf("  Recv end  :(%03d,%03d,%03d)\n", p->hv[1].recv_e[0], p->hv[1].recv_e[1], p->hv[1].recv_e[2]); fflush(stdout);
        printf("  Send begin:(%03d,%03d,%03d)\n", p->hv[1].send_b[0], p->hv[1].send_b[1], p->hv[1].send_b[2]); fflush(stdout);
        printf("  Send end  :(%03d,%03d,%03d)\n", p->hv[1].send_e[0], p->hv[1].send_e[1], p->hv[1].send_e[2]); fflush(stdout);
        printf("\n"); fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if(p->mpi_rank == 0) {
      printf("******************************************************\n"); fflush(stdout);
      printf("DEBUG v halo information END\n"); fflush(stdout);
      printf("******************************************************\n\n"); fflush(stdout);
    }
  }

  ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->hv[1].shape, p->hv[1].recv_b, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->hv[1].recv_hb)); CHKERR(ierr);
  ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->hv[1].shape, p->hv[1].recv_e  , MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->hv[1].recv_he)); CHKERR(ierr);
  ierr = MPI_Type_commit(&(p->hv[1].recv_hb)); CHKERR(ierr);
  ierr = MPI_Type_commit(&(p->hv[1].recv_he)); CHKERR(ierr);

  ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->hv[1].shape, p->hv[1].send_b, MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->hv[1].send_hb)); CHKERR(ierr);
  ierr = MPI_Type_create_subarray(3, p->ldomain_shape, p->hv[1].shape, p->hv[1].send_e  , MPI_ORDER_FORTRAN, MPI_FLOAT_PRECISION, &(p->hv[1].send_he)); CHKERR(ierr);
  ierr = MPI_Type_commit(&(p->hv[1].send_hb)); CHKERR(ierr);
  ierr = MPI_Type_commit(&(p->hv[1].send_he)); CHKERR(ierr);

}
void mpi_halo_init(Parameters *p) {
  switch(p->target_ts){
  case 0: //standard methods
  case 1:
    standard_mpi_halo_init(p);
    break;
  case 2: // intra-diamond methods
    intra_diamond_mpi_halo_init(p);
    break;
  }
}


void standard_mpi_halo_finalize(Parameters *p){
  int i;
  for(i=0; i<2; i++){
    if (p->t.shape[i] > 1) {
      ierr = MPI_Type_free(&(p->h[i].recv_hb)); CHKERR(ierr);
      ierr = MPI_Type_free(&(p->h[i].recv_he)); CHKERR(ierr);
      ierr = MPI_Type_free(&(p->h[i].send_hb)); CHKERR(ierr);
      ierr = MPI_Type_free(&(p->h[i].send_he)); CHKERR(ierr);
    }
  }
  if (p->t.shape[2] > 1)
    ierr = MPI_Type_free(&(p->h[i].halo)); CHKERR(ierr);

}
void intra_diamond_mpi_halo_finalize(Parameters *p){
  if(p->t.shape[1] > 1){
    ierr = MPI_Type_free(&(p->hu[1].recv_hb)); CHKERR(ierr);
    ierr = MPI_Type_free(&(p->hu[1].recv_he)); CHKERR(ierr);
    ierr = MPI_Type_free(&(p->hu[1].send_hb)); CHKERR(ierr);
    ierr = MPI_Type_free(&(p->hu[1].send_he)); CHKERR(ierr);

    ierr = MPI_Type_free(&(p->hv[1].recv_hb)); CHKERR(ierr);
    ierr = MPI_Type_free(&(p->hv[1].recv_he)); CHKERR(ierr);
    ierr = MPI_Type_free(&(p->hv[1].send_hb)); CHKERR(ierr);
    ierr = MPI_Type_free(&(p->hv[1].send_he)); CHKERR(ierr);
  }
}
void mpi_halo_finalize(Parameters *p){
  switch(p->target_ts){
  case 0: //standard methods
  case 1:
    standard_mpi_halo_finalize(p);
    break;
  case 2: // intra-diamond methods
    intra_diamond_mpi_halo_finalize(p);
    break;
  }
}

