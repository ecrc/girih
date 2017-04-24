#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"
#include "data_structures.h"

#define U(i,j,k)          (    u[((k)*(p->ldomain_shape[1])+(j))*(p->ldomain_shape[0])+(i)])
#define U1(i,j,k)         (p->U1[((k)*(p->ldomain_shape[1])+(j))*(p->ldomain_shape[0])+(i)])
#define U2(i,j,k)         (p->U2[((k)*(p->ldomain_shape[1])+(j))*(p->ldomain_shape[0])+(i)])

extern void sub_array_copy(const real_t * restrict src_buf, real_t * restrict dst_buf, int *src_size, int *dst_size, int *cpy_size, int *src_offs, int *dst_offs);

static inline void exchange_halo_srtided_asynch(Parameters *p, real_t * restrict u) {
  MPI_Request wait_req[3][4];
  MPI_Status wait_stat[3][4];
  int i;

  // Asynch receive in X
  if (p->t.shape[0] > 1) {
    ierr = MPI_Irecv(u, 1, p->h[0].recv_hb, p->t.left , MPI_ANY_TAG, p->t.cart_comm, &(wait_req[0][0])); CHKERR(ierr);
    ierr = MPI_Irecv(u, 1, p->h[0].recv_he  , p->t.right, MPI_ANY_TAG, p->t.cart_comm, &(wait_req[0][1])); CHKERR(ierr);
  }
  // Asynch receive in Y
  if (p->t.shape[1] > 1) {
    ierr = MPI_Irecv(u, 1, p->h[1].recv_hb, p->t.down , MPI_ANY_TAG, p->t.cart_comm, &(wait_req[1][0])); CHKERR(ierr);
    ierr = MPI_Irecv(u, 1, p->h[1].recv_he  , p->t.up, MPI_ANY_TAG, p->t.cart_comm, &(wait_req[1][1])); CHKERR(ierr);
  }
  // Asynch receive in Z
  if (p->t.shape[2] > 1) {
    ierr = MPI_Irecv(&(u[p->h[2].recv_b[2]]), 1, p->h[2].halo, p->t.back , MPI_ANY_TAG, p->t.cart_comm, &(wait_req[2][0])); CHKERR(ierr);
    ierr = MPI_Irecv(&(u[p->h[2].recv_e[2]]), 1, p->h[2].halo, p->t.front, MPI_ANY_TAG, p->t.cart_comm, &(wait_req[2][1])); CHKERR(ierr);
  }


  // Asynch send in X
  if (p->t.shape[0] > 1) {
    ierr = MPI_Isend(u, 1, p->h[0].send_hb, p->t.left , 0, p->t.cart_comm, &(wait_req[0][2])); CHKERR(ierr);
    ierr = MPI_Isend(u, 1, p->h[0].send_he  , p->t.right, 0, p->t.cart_comm, &(wait_req[0][3])); CHKERR(ierr);
  }
  // Asynch send in Y
  if (p->t.shape[1] > 1) {
    ierr = MPI_Isend(u, 1, p->h[1].send_hb, p->t.down , 0, p->t.cart_comm, &(wait_req[1][2])); CHKERR(ierr);
    ierr = MPI_Isend(u, 1, p->h[1].send_he  , p->t.up, 0, p->t.cart_comm, &(wait_req[1][3])); CHKERR(ierr);
  }
  // Asynch send in Z
  if (p->t.shape[2] > 1) {
    ierr = MPI_Isend(&(u[p->h[2].send_b[2]]), 1, p->h[2].halo, p->t.back , 0, p->t.cart_comm, &(wait_req[2][2])); CHKERR(ierr);
    ierr = MPI_Isend(&(u[p->h[2].send_e[2]]), 1, p->h[2].halo, p->t.front, 0, p->t.cart_comm, &(wait_req[2][3])); CHKERR(ierr);
  }


  // Wait all communication to complete
  for(i=0; i<3; i++)
    if (p->t.shape[i] > 1) ierr = MPI_Waitall(4, wait_req[i], wait_stat[i]); CHKERR(ierr);
}
static inline void exchange_halo_concat_asynch(Parameters *p, real_t * restrict u, real_t * restrict x_send_buf, real_t * restrict x_recv_buf, real_t * restrict y_send_buf, real_t * restrict y_recv_buf) {
  MPI_Request wait_req[3][4];
  MPI_Status wait_stat[4];
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
static inline void exchange_halo_asynch(Parameters *p, real_t * restrict u, real_t * restrict x_send_buf, real_t * restrict x_recv_buf, real_t * restrict y_send_buf, real_t * restrict y_recv_buf) {
  if(p->halo_concat == 1){
    exchange_halo_concat_asynch(p, u, x_send_buf, x_recv_buf, y_send_buf, y_recv_buf);
  } else {
    exchange_halo_srtided_asynch(p, u);
  }
}


void naive_nonblocking_ts(Parameters *p) {
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

  double t1,t2,t3,t4,t5;
  for(it=0; it<p->nt; it+=2){
    t1 = MPI_Wtime();
    p->stencil.spt_blk_func(p->ldomain_shape, p->stencil.r, p->stencil.r, p->stencil.r, p->lstencil_shape[0]+p->stencil.r, p->ldomain_shape[1]-p->stencil.r, p->ldomain_shape[2]-p->stencil.r, p->coef, p->U1, p->U2, p->U3, ALL_FIELDS, p->stencil_ctx);
    if( (p->has_source==1) && (p->source_point_enabled==1)) U1(p->lsource_pt[0],p->lsource_pt[1],p->lsource_pt[2]) += p->src_exc_coef[it]; //@KADIR Source exicitations are moved before stencil computations
    printf("%s %d: Source updated: p->has_source:%d p->source_point_enabled:%d\n", __FILE__, __LINE__, p->has_source, p->source_point_enabled); //@KADIR
    t2 = MPI_Wtime();
    exchange_halo_asynch(p, p->U1, x_send_buf, x_recv_buf, y_send_buf, y_recv_buf);

    t3 = MPI_Wtime();

    p->stencil.spt_blk_func(p->ldomain_shape, p->stencil.r, p->stencil.r, p->stencil.r, p->lstencil_shape[0]+p->stencil.r, p->ldomain_shape[1]-p->stencil.r, p->ldomain_shape[2]-p->stencil.r, p->coef, p->U2, p->U1, p->U3, ALL_FIELDS, p->stencil_ctx);
    if( (p->has_source==1) && (p->source_point_enabled==1)) U2(p->lsource_pt[0],p->lsource_pt[1],p->lsource_pt[2]) += p->src_exc_coef[it+1]; //@KADIR Source exicitations are moved before stencil computations
    t4 = MPI_Wtime();
    exchange_halo_asynch(p, p->U2, x_send_buf, x_recv_buf, y_send_buf, y_recv_buf);
    t5 = MPI_Wtime();

    p->prof.communicate += (t3 - t2) + (t5 - t4);
    p->prof.compute += (t2-t1) + (t4-t3);
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
