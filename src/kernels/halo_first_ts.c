#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include "data_structures.h"

#define U(i,j,k)          (    u[((k)*(p->ldomain_shape[1])+(j))*(p->ldomain_shape[0])+(i)])

extern void sub_array_copy(const real_t * restrict src_buf, real_t * restrict dst_buf, int *src_size, int *dst_size, int *cpy_size, int *src_offs, int *dst_offs);

static inline void exchange_halo_concat_start(Parameters *p, real_t * restrict u, MPI_Request wait_req[3][4], real_t * restrict x_send_buf,
    real_t * restrict x_recv_buf, real_t * restrict y_send_buf, real_t * restrict y_recv_buf) {
  int i;

  int z_offs[] = {0,0,0};

  // Asynch receive in X
  if (p->t.shape[0] > 1) {
    ierr = MPI_Irecv(x_recv_buf                 , p->h[0].size, MPI_real_t, p->t.left , MPI_ANY_TAG, p->t.cart_comm, &(wait_req[0][0])); CHKERR(ierr);
    ierr = MPI_Irecv(&(x_recv_buf[p->h[0].size]), p->h[0].size, MPI_real_t, p->t.right, MPI_ANY_TAG, p->t.cart_comm, &(wait_req[0][1])); CHKERR(ierr);
  }
  // Asynch receive in Y
  if (p->t.shape[1] > 1) {
    ierr = MPI_Irecv(y_recv_buf                 , p->h[1].size, MPI_real_t, p->t.down , MPI_ANY_TAG, p->t.cart_comm, &(wait_req[1][0])); CHKERR(ierr);
    ierr = MPI_Irecv(&(y_recv_buf[p->h[1].size]), p->h[1].size, MPI_real_t, p->t.up, MPI_ANY_TAG, p->t.cart_comm, &(wait_req[1][1])); CHKERR(ierr);
  }
  // Asynch receive in Z
  if (p->t.shape[2] > 1) {
    ierr = MPI_Irecv(&(u[p->h[2].recv_b[2]]), 1, p->h[2].halo, p->t.back , MPI_ANY_TAG, p->t.cart_comm, &(wait_req[2][0])); CHKERR(ierr);
    ierr = MPI_Irecv(&(u[p->h[2].recv_e[2]]), 1, p->h[2].halo, p->t.front, MPI_ANY_TAG, p->t.cart_comm, &(wait_req[2][1])); CHKERR(ierr);
  }


  // Asynch send in X
  if (p->t.shape[0] > 1) {
    if( p->t.left != MPI_PROC_NULL)
      sub_array_copy(u, &(x_send_buf[p->h[0].size]), p->ldomain_shape, p->h[0].shape, p->h[0].shape, p->h[0].send_b, z_offs);
    ierr = MPI_Isend(&(x_send_buf[p->h[0].size]), p->h[0].size, MPI_real_t, p->t.left , 0, p->t.cart_comm, &(wait_req[0][2])); CHKERR(ierr);

    if( p->t.right != MPI_PROC_NULL)
      sub_array_copy(u, x_send_buf, p->ldomain_shape, p->h[0].shape, p->h[0].shape, p->h[0].send_e, z_offs);
    ierr = MPI_Isend(x_send_buf, p->h[0].size, MPI_real_t, p->t.right, 0, p->t.cart_comm, &(wait_req[0][3])); CHKERR(ierr);
  }
  // Asynch send in Y
  if (p->t.shape[1] > 1) {
    if( p->t.down != MPI_PROC_NULL)
      sub_array_copy(u, &(y_send_buf[p->h[1].size]), p->ldomain_shape, p->h[1].shape, p->h[1].shape, p->h[1].send_b, z_offs);
    ierr = MPI_Isend(&(y_send_buf[p->h[1].size]), p->h[1].size, MPI_real_t, p->t.down , 0, p->t.cart_comm, &(wait_req[1][2])); CHKERR(ierr);

    if( p->t.up != MPI_PROC_NULL)
      sub_array_copy(u, y_send_buf, p->ldomain_shape, p->h[1].shape, p->h[1].shape, p->h[1].send_e, z_offs);
    ierr = MPI_Isend(y_send_buf, p->h[1].size, MPI_real_t, p->t.up, 0, p->t.cart_comm, &(wait_req[1][3])); CHKERR(ierr);
  }
  // Asynch send in Z
  if (p->t.shape[2] > 1) {
    ierr = MPI_Isend(&(u[p->h[2].send_b[2]]), 1, p->h[2].halo, p->t.back , 0, p->t.cart_comm, &(wait_req[2][2])); CHKERR(ierr);
    ierr = MPI_Isend(&(u[p->h[2].send_e[2]]), 1, p->h[2].halo, p->t.front, 0, p->t.cart_comm, &(wait_req[2][3])); CHKERR(ierr);
  }
}
static inline void exchange_halo_srtided(Parameters *p, real_t * restrict u, MPI_Request wait_req[3][4]) {
  if(p->t.shape[0]>1){
    ierr = MPI_Irecv(u, 1, p->h[0].recv_hb, p->t.left , MPI_ANY_TAG, p->t.cart_comm, &(wait_req[0][0])); CHKERR(ierr);
    ierr = MPI_Irecv(u, 1, p->h[0].recv_he, p->t.right, MPI_ANY_TAG, p->t.cart_comm, &(wait_req[0][1])); CHKERR(ierr);

    ierr = MPI_Isend(u, 1, p->h[0].send_hb, p->t.left , 0, p->t.cart_comm, &(wait_req[0][2])); CHKERR(ierr);
    ierr = MPI_Isend(u, 1, p->h[0].send_he, p->t.right, 0, p->t.cart_comm, &(wait_req[0][3])); CHKERR(ierr);
  }
  if(p->t.shape[1]>1){
    ierr = MPI_Irecv(u, 1, p->h[1].recv_hb, p->t.down, MPI_ANY_TAG, p->t.cart_comm, &(wait_req[1][0])); CHKERR(ierr);
    ierr = MPI_Irecv(u, 1, p->h[1].recv_he, p->t.up  , MPI_ANY_TAG, p->t.cart_comm, &(wait_req[1][1])); CHKERR(ierr);

    ierr = MPI_Isend(u, 1, p->h[1].send_hb, p->t.down, 0, p->t.cart_comm, &(wait_req[1][2])); CHKERR(ierr);
    ierr = MPI_Isend(u, 1, p->h[1].send_he, p->t.up  , 0, p->t.cart_comm, &(wait_req[1][3])); CHKERR(ierr);
  }
  if(p->t.shape[2]>1){
    ierr = MPI_Irecv(&(u[p->h[2].recv_b[2]]), 1, p->h[2].halo, p->t.back , MPI_ANY_TAG, p->t.cart_comm, &(wait_req[2][0])); CHKERR(ierr);
    ierr = MPI_Irecv(&(u[p->h[2].recv_e[2]]), 1, p->h[2].halo, p->t.front, MPI_ANY_TAG, p->t.cart_comm, &(wait_req[2][1])); CHKERR(ierr);

    ierr = MPI_Isend(&(u[p->h[2].send_b[2]]), 1, p->h[2].halo, p->t.back , 0, p->t.cart_comm, &(wait_req[2][2])); CHKERR(ierr);
    ierr = MPI_Isend(&(u[p->h[2].send_e[2]]), 1, p->h[2].halo, p->t.front, 0, p->t.cart_comm, &(wait_req[2][3])); CHKERR(ierr);
  }
}
static inline void exchange_halo_start(Parameters *p, real_t * restrict u, MPI_Request wait_req[3][4], real_t * restrict x_send_buf,
    real_t * restrict x_recv_buf, real_t * restrict y_send_buf, real_t * restrict y_recv_buf) {
  if(p->halo_concat == 1){
    exchange_halo_concat_start(p, u, wait_req, x_send_buf, x_recv_buf, y_send_buf, y_recv_buf);
  } else {
    exchange_halo_srtided(p, u, wait_req);
  }
}

static inline void exchange_halo_concat_finish(Parameters *p, real_t * restrict u, MPI_Request wait_req[3][4],
    real_t * restrict x_recv_buf, real_t * restrict y_recv_buf) {
  MPI_Status wait_stat[4];
  int z_offs[] = {0,0,0};

  // Wait communication across X to complete
  if (p->t.shape[0] > 1) {
    // wait for data to arrive and copy to the strided buffer
    ierr = MPI_Waitall(2, wait_req[0], wait_stat); CHKERR(ierr);
    if( p->t.left != MPI_PROC_NULL)
      sub_array_copy(x_recv_buf,                  u, p->h[0].shape, p->ldomain_shape, p->h[0].shape, z_offs, p->h[0].recv_b);
    if( p->t.right != MPI_PROC_NULL)
      sub_array_copy(&(x_recv_buf[p->h[0].size]), u, p->h[0].shape, p->ldomain_shape, p->h[0].shape, z_offs, p->h[0].recv_e);

    // wait for data to be sent
    ierr = MPI_Waitall(2, &(wait_req[0][2]), wait_stat); CHKERR(ierr);
  }

  // Wait communication across Y to complete
  if (p->t.shape[1] > 1){
    // wait for data to arrive and copy to the strided buffer
    ierr = MPI_Waitall(2, wait_req[1], wait_stat); CHKERR(ierr);
    if( p->t.down != MPI_PROC_NULL)
      sub_array_copy(y_recv_buf,                  u, p->h[1].shape, p->ldomain_shape, p->h[1].shape, z_offs, p->h[1].recv_b);
    if( p->t.up != MPI_PROC_NULL)
      sub_array_copy(&(y_recv_buf[p->h[1].size]), u, p->h[1].shape, p->ldomain_shape, p->h[1].shape, z_offs, p->h[1].recv_e);

    // wait for data to be sent
    ierr = MPI_Waitall(2, &(wait_req[1][2]), wait_stat); CHKERR(ierr);
  }

  // Wait communication across Z to complete
  if (p->t.shape[2] > 1) ierr = MPI_Waitall(4, wait_req[2], wait_stat); CHKERR(ierr);
}
static inline void exchange_halo_finish(Parameters *p, real_t * restrict u, MPI_Request wait_req[3][4], real_t * restrict x_recv_buf, real_t * restrict y_recv_buf) {
  int x;
  MPI_Status wait_stat[4];
  if(p->halo_concat == 1){
    exchange_halo_concat_finish(p, u, wait_req, x_recv_buf, y_recv_buf);
  } else {
    for(x=0; x<3; x++)
      if (p->t.shape[x] > 1) ierr = MPI_Waitall(4, wait_req[x], wait_stat); CHKERR(ierr);
  }
}

static inline void halo_first_step(Parameters *p, int it, real_t * restrict u, real_t * restrict v,
    int start_b[3][3], int start_e[3][3], int finish_b[3][3], int finish_e[3][3], int middle_b[3], int middle_e[3],
    int side_source, int middle_source, real_t * restrict x_send_buf, real_t * restrict x_recv_buf,
    real_t * restrict y_send_buf, real_t * restrict y_recv_buf) {
  double t1, t2, t3, t4;
  int x;
  MPI_Request wait_req[3][4];

  t1 = MPI_Wtime();

  for(x=0; x<3; x++){
    if(p->t.shape[x]>1){
      // compute begin
      p->stencil.spt_blk_func(p->ldomain_shape, start_b[x][0], start_b[x][1], start_b[x][2], start_e[x][0], start_e[x][1], start_e[x][2], p->coef, u, v, p->U3, ALL_FIELDS, p->stencil_ctx);
      // compute finish
      p->stencil.spt_blk_func(p->ldomain_shape, finish_b[x][0], finish_b[x][1], finish_b[x][2], finish_e[x][0], finish_e[x][1], finish_e[x][2], p->coef, u, v, p->U3, ALL_FIELDS, p->stencil_ctx);
    }
  }
  if( (side_source==1) && (p->source_point_enabled==1)) U(p->lsource_pt[0],p->lsource_pt[1],p->lsource_pt[2]) += p->source[it];
  t2 = MPI_Wtime();

  // communicate halo
  exchange_halo_start(p, u, wait_req, x_send_buf, x_recv_buf, y_send_buf, y_recv_buf);

  t3 = MPI_Wtime();
  // compute middle
  p->stencil.spt_blk_func(p->ldomain_shape, middle_b[0], middle_b[1], middle_b[2], middle_e[0], middle_e[1], middle_e[2], p->coef, u, v, p->U3, ALL_FIELDS, p->stencil_ctx);
  if( (middle_source==1) && (p->source_point_enabled==1)) U(p->lsource_pt[0],p->lsource_pt[1],p->lsource_pt[2]) += p->source[it];

  t4 = MPI_Wtime();

  // wait for the rest of communication
  exchange_halo_finish(p, u, wait_req, x_recv_buf, y_recv_buf);

  p->prof.wait += (MPI_Wtime() - t4);
  p->prof.compute += (t2-t1 + t4-t3);
  p->prof.communicate += (t3-t2);

}

void halo_first_ts(Parameters *p) {
  int it;

  // create buffers to aggregate halo data for communication
  real_t *x_send_buf, *x_recv_buf, *y_send_buf, *y_recv_buf;
  int comm_buf_size;
  if (p->halo_concat == 1){
    if (p->t.shape[0] > 1){
      // assuming same halo size for both U and V buffers
      comm_buf_size = 2 * p->h[0].size;
      posix_memalign((void **)&(x_recv_buf), p->alignment, sizeof(real_t)*comm_buf_size);
      posix_memalign((void **)&(x_send_buf), p->alignment, sizeof(real_t)*comm_buf_size);
    }
    if (p->t.shape[1] > 1){
      // assuming same halo size for both U and V buffers
      comm_buf_size = 2 * p->h[1].size;
      posix_memalign((void **)&(y_recv_buf), p->alignment, sizeof(real_t)*comm_buf_size);
      posix_memalign((void **)&(y_send_buf), p->alignment, sizeof(real_t)*comm_buf_size);
    }
  }

  int start_b[3][3], start_e[3][3], finish_b[3][3], finish_e[3][3], middle_b[3], middle_e[3];
  int send_b[3][3],send_e[3][3], shape[3][3];
  int x,y;

  // sides shapes and boundaries
  send_b[0][0]= p->stencil.r; send_e[0][0]= p->lstencil_shape[0];
  send_b[0][1]= p->stencil.r; send_e[0][1]= p->stencil.r;
  send_b[0][2]= p->stencil.r; send_e[0][2]= p->stencil.r;
  send_b[1][0]= p->stencil.r; send_e[1][0]= p->stencil.r;
  send_b[1][1]= p->stencil.r; send_e[1][1]= p->ldomain_shape[1]-2*p->stencil.r;
  send_b[1][2]= p->stencil.r; send_e[1][2]= p->stencil.r;
  send_b[2][0]= p->stencil.r; send_e[2][0]= p->stencil.r;
  send_b[2][1]= p->stencil.r; send_e[2][1]= p->stencil.r;
  send_b[2][2]= p->stencil.r; send_e[2][2]= p->ldomain_shape[2]-2*p->stencil.r;
  shape[0][0]= p->stencil.r;
  shape[0][1]= p->ldomain_shape[1]-2*p->stencil.r;
  shape[0][2]= p->ldomain_shape[2]-2*p->stencil.r;
  shape[1][0]= p->lstencil_shape[0];
  shape[1][1]= p->stencil.r;
  shape[1][2]= p->ldomain_shape[2]-2*p->stencil.r;
  shape[2][0]= p->lstencil_shape[0];
  shape[2][1]= p->ldomain_shape[1]-2*p->stencil.r;
  shape[2][2]= p->stencil.r;

  // compute the coordinates of the sides computations
  for(x=0; x<3; x++){
    for(y=0; y<3; y++){
      start_b[x][y] = send_b[x][y];
      start_e[x][y] = send_b[x][y] + shape[x][y];
      finish_b[x][y] = send_e[x][y];
      finish_e[x][y] = send_e[x][y] + shape[x][y];
    }
  }
  // trim the overlapping sides computations across dimensions
  for(x=0; x<3; x++){
    for(y=0; y<3; y++){
      if(p->t.shape[y]>1){
        if(y>x){ // domain dimension > topology dimension
          start_b[x][y] += p->stencil.r;
          start_e[x][y] -= p->stencil.r;

          finish_b[x][y] += p->stencil.r;
          finish_e[x][y] -= p->stencil.r;
        }
      }
    }
  }

  // compute the coordinates of the middle computations
  for(x=0; x<3; x++){
    if(p->t.shape[x]>1){
      middle_b[x] = start_e[x][x];
      middle_e[x] = finish_b[x][x];
    } else {
      middle_b[x] = p->stencil.r;
      middle_e[x] = p->lstencil_shape[x] + p->stencil.r;
    }
  }

  // locate the source point
  int side_source=0, middle_source=0;
  if(p->has_source==1) {
    if( (p->lsource_pt[0] <  2*p->stencil.r) || (p->lsource_pt[0] >= p->lstencil_shape[0]) ||
        (p->lsource_pt[1] <  2*p->stencil.r) || (p->lsource_pt[1] >= p->ldomain_shape[1]- 2*p->stencil.r) ||
        (p->lsource_pt[2] <  2*p->stencil.r) || (p->lsource_pt[2] >= p->ldomain_shape[2]- 2*p->stencil.r) )
      side_source=1;
    else
      middle_source=1;
  }

//  //Debug
//  for(x=0; x<3; x++){
//    if(p->mpi_rank ==0){
//      if(p->t.shape[x]>1){
//        printf("start[%d]:  (%03d, %03d, %03d)-(%03d, %03d, %03d)\n", x,
//            start_b[x][0],start_b[x][1],start_b[x][2],
//            start_e[x][0],start_e[x][1],start_e[x][2]);
//        printf("finish[%d]: (%03d, %03d, %03d)-(%03d, %03d, %03d)\n\n", x,
//            finish_b[x][0],finish_b[x][1],finish_b[x][2],
//            finish_e[x][0],finish_e[x][1],finish_e[x][2]);
//      }
//    }
//  }
//  if(p->mpi_rank ==0){
//    printf("Middle:   (%03d, %03d, %03d)-(%03d, %03d, %03d)\n",
//        middle_b[0],middle_b[1],middle_b[2],
//        middle_e[0],middle_e[1],middle_e[2]);
//    printf("side_source:%d | middle_source:%d\n", side_source, middle_source);
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//  MPI_Finalize();
//  exit(1);

  for(it=0; it<p->nt; it+=2){
    halo_first_step(p, it, p->U1, p->U2,
        start_b, start_e, finish_b, finish_e, middle_b, middle_e, side_source, middle_source,
        x_send_buf, x_recv_buf, y_send_buf, y_recv_buf);
    halo_first_step(p, it+1, p->U2, p->U1,
        start_b, start_e, finish_b, finish_e, middle_b, middle_e, side_source, middle_source,
        x_send_buf, x_recv_buf, y_send_buf, y_recv_buf);
  }

  // clean up the buffers that aggregate halo data for communication
  if (p->halo_concat == 1){
    if (p->t.shape[0] > 1){
      free(x_recv_buf);
      free(x_send_buf);
    }
    if (p->t.shape[1] > 1){
      free(y_recv_buf);
      free(y_send_buf);
    }
  }
}
