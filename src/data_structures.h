#ifndef DATA_STRUCTURES_H_
#define DATA_STRUCTURES_H_

#ifdef LIKWID_PERFMON
#include "likwid.h"
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_START(reg)
#define LIKWID_MARKER_STOP(reg)
#define LIKWID_MARKER_CLOSE
#endif


#ifndef _MPI_INCLUDE
#include "mpi.h"
#endif

// For thread binding
#define _GNU_SOURCE
#define __USE_GNU
#include <sched.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <errno.h>

#ifndef restrict
#define restrict __restrict__
#endif

// constants
#define MAX_CACHE_SIZE (70*1024)  // cache size in kB
#define FLUSH_SIZE (50*1024*1024) // in bytes
#define BOUNDARY_SRC_VAL (100.1)

#ifndef ENABLE_CPU_BIND
#define ENABLE_CPU_BIND (1)
#endif

#ifndef DP
#define DP (0)
#endif

// Experimental features
#ifndef USE_SPLIT_STRIDE
#define USE_SPLIT_STRIDE (0)
#endif

#if DP
typedef double real_t;
#define MPI_real_t MPI_DOUBLE
#else
typedef float real_t;
#define MPI_real_t MPI_FLOAT
#endif

#define CHKERR(error_stat) {elen=0; estr[0]='d';}//do{ if(error_stat != MPI_SUCCESS){MPI_Error_string(ierr,estr,&elen); fprintf(stderr,estr); MPI_Finalize();} } while(0)

#include <stdint.h>

int ierr;
char estr[MPI_MAX_ERROR_STRING]; int elen; // for MPI error hadling macro

enum Stencil_Shapes{
  STAR,
  TTI,
  BOX
};

enum Stencil_Coefficients{
  CONSTANT_COEFFICIENT,
  VARIABLE_COEFFICIENT,
  VARIABLE_COEFFICIENT_AXSYM,
  VARIABLE_COEFFICIENT_NOSYM,
  SOLAR_COEFFICIENT
};

enum Stencil_Type{
  REGULAR,
  SOLAR
};

// Fields in solar kernel grid cell
enum Solar_Fields{
  ALL_FIELDS,
  H_FIELD,
  E_FIELD,
};

// Profiling
typedef struct{
  double compute, communicate, send_recv, wait, total, others, ts_main, ts_others;
}Profile;

// Halo information
typedef struct{
  int shape[3], recv_b[3], recv_e[3], send_b[3], send_e[3], is_contiguous, size;
  MPI_Datatype recv_hb, recv_he, send_hb, send_he;

  // new datatypes for the new structure
  MPI_Datatype halo;
}Halo;

// MPI topology information
typedef struct{
  int right, // x+
      left,  // x-
      up,    // y+
      down,  // y-
      front, // z+
      back;  // z-
  int shape[3], is_periodic[3], rank_coords[3];
  MPI_Comm cart_comm;
//  MPI_Request wait_req[8];
//  MPI_Status wait_stat16[16],wait_stat8[8],wait_stat4[4];
}mpi_topology;

typedef struct{
  int nnx, nny, nnz;
  uint64_t ln_domain;
}CLU_CTX;

#define CLU_SIG (const CLU_CTX clu_ctx, const int xb, const int xe, const int j, const int k, \
const real_t * restrict coef, real_t * restrict u, \
const real_t * restrict v, const real_t * restrict roc2)

typedef void (*clu_func_t)CLU_SIG;

// context information passed to the stencil kernel
typedef struct{
  int bs_x; //depricated
  int bs_y; // for spatial blocking in Y at the standard methods
  int thread_group_size;
  int th_x, th_y, th_z, th_c; // number of threads per dimension in x, y, and z, and per component

  // cpu binding masks
  cpu_set_t **bind_masks;
  int setsize;

  // for separate stride-1 functions
  clu_func_t clu_func;

  int num_wf; // number of wavefront updats per iteration
  // wavefront profiling
  double *t_wf_comm, *t_wait, *t_wf_prologue, *t_wf_main, *t_wf_epilogue, *wf_num_resolved_diamonds, *t_group_wait;

  // likwid markers usage enable/disable 
  int enable_likwid_m;
}stencil_CTX;


// Kernels and time steppers data structures
#define KERNEL_SIG     ( const int shape[3], const int xb, const int yb,  const int zb, const int xe, const int ye, const int ze,\
    const real_t * restrict coef, real_t * restrict u, const real_t * restrict v, const real_t * restrict roc2, int field, stencil_CTX stencil_ctx)
#define KERNEL_MWD_SIG ( const int shape[3], const int xb, const int yb_r, const int zb, const int xe, const int ye_r, const int ze, \
    const real_t * restrict coef, real_t * restrict u, \
    real_t * restrict v, const real_t * restrict roc2, int t_dim, int b_inc, int e_inc, int NHALO, int tb, int te, stencil_CTX stencil_ctx, int mtid)

typedef void (*spt_blk_func_t)KERNEL_SIG;
typedef void (*mwd_func_t)KERNEL_MWD_SIG;

struct Stencil {
  const char *name;
  int r;
  int time_order;
  int nd;
  enum Stencil_Shapes shape;
  enum Stencil_Coefficients coeff;
  enum Stencil_Type type;

  spt_blk_func_t spt_blk_func;
  spt_blk_func_t stat_sched_func;
  mwd_func_t mwd_func;
};

struct StencilInfo {
  const char *name;
  int r;
  int time_order;
  int nd;
  enum Stencil_Shapes shape;
  enum Stencil_Coefficients coeff;
  enum Stencil_Type type;
};

// context information
typedef struct{
  int alignment, verbose, stencil_shape[3];
  uint64_t n_stencils, ln_domain, ln_stencils;
  int target_ts, target_kernel;
  int mpi_rank, mpi_size;
  int n_tests, nt;
  int verify;
  int source_pt[3];
  int debug;
  int num_threads;
  int use_omp_stat_sched;
  int lstencil_shape[3], ldomain_shape[3], gb[3], ge[3], lsource_pt[3], has_source; //MPI ranks' global indices, and local source locations

//  int stencil_radius, is_constant_coefficient;
//  enum Stencil_Types stencil_type;

  real_t * restrict U1, * restrict U2, * restrict U3, * restrict source;
  real_t * restrict coef;

  // Holds the value of cache blocking across Y axis
  stencil_CTX stencil_ctx;

  // to enable/disable source point update
  int source_point_enabled;

  int cache_size; // Last level cache usable size in KB for blocking

  // to enable concatenating halo information before communication
  int halo_concat;

  // Specific data for the diamond method
  int t_dim, larger_t_dim, is_last, mwd_type;
  Halo hu[3], hv[3];

  int wavefront;
  uint64_t idiamond_pro_epi_logue_updates;
  uint64_t wf_blk_size, wf_larger_blk_size;

  Halo h[3]; // Halo information for X,Y, and Z directions
  mpi_topology t;
  Profile prof;

  struct Stencil stencil;

  // list of coefficients to be used in stencil operators
  real_t g_coef[11];

  int array_padding;

}Parameters;

struct time_stepper {
  const char *name;
  void (*func)(Parameters *p);
};

#endif /* DATA_STRUCTURES_H_ */
