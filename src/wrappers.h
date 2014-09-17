#ifndef WRAPPERS_H_
#define WRAPPERS_H_

#include "data_structures.h"

// 1 core multi-wavefront kernels
extern void swd_iso_ref KERNEL_MWD_SIG;
extern void swd_iso_ref_2space_1time KERNEL_MWD_SIG;
extern void swd_iso_ref_2space_1time_var KERNEL_MWD_SIG;
extern void swd_iso_ref_2space_1time_var_axsym KERNEL_MWD_SIG;
extern void swd_iso_ref_8space_1time_var_axsym KERNEL_MWD_SIG;
extern void swd_iso_ref_2space_1time_var_nosym KERNEL_MWD_SIG;

// standard kernels

// naive cases with static "chunk" openmp schedule
extern void stat_sched_iso_ref KERNEL_SIG;
extern void stat_sched_iso_ref_2space_1time KERNEL_SIG;
extern void stat_sched_iso_ref_2space_1time_var KERNEL_SIG;
extern void stat_sched_iso_ref_2space_1time_var_axsym KERNEL_SIG;
extern void stat_sched_iso_ref_8space_1time_var_axsym KERNEL_SIG;
extern void stat_sched_iso_ref_2space_1time_var_nosym KERNEL_SIG;

#if USE_SPLIT_STRIDE
// naive cases with split function call to the 1-stride
extern void iso_ref_split KERNEL_SIG;

#define iso_refi iso_ref_split
#define iso_ref_2space_1timei iso_ref_split
#define iso_ref_2space_1time_vari iso_ref_split
#define iso_ref_2space_1time_var_axsymi iso_ref_split

#else
// naive cases
extern void iso_ref KERNEL_SIG;
extern void iso_ref_2space_1time KERNEL_SIG;
extern void iso_ref_2space_1time_var KERNEL_SIG;
extern void iso_ref_2space_1time_var_axsym KERNEL_SIG;
extern void iso_ref_8space_1time_var_axsym KERNEL_SIG;
extern void iso_ref_2space_1time_var_nosym KERNEL_SIG;

#define iso_refi iso_ref
#define iso_ref_2space_1timei iso_ref_2space_1time
#define iso_ref_2space_1time_vari iso_ref_2space_1time_var
#define iso_ref_2space_1time_var_axsymi iso_ref_2space_1time_var_axsym
#define iso_ref_8space_1time_var_axsymi iso_ref_8space_1time_var_axsym
#define iso_ref_2space_1time_var_nosymi iso_ref_2space_1time_var_nosym
#endif


extern void naive_nonblocking_ts(Parameters *p);
extern void halo_first_ts(Parameters *p);
extern void dynamic_intra_diamond_ts(Parameters *p);

struct time_stepper TSList[] = {
/* 0 */    {"Spatial Blocking", naive_nonblocking_ts},
    // Naive Halo exchange implementation with non-locking MPI communication
/* 1 */    {"Halo-first", halo_first_ts},
    // Compute the sides then communicate the halo while computing the middle
/* 2 */    {"Diamond", dynamic_intra_diamond_ts},
    // Dynamic scheduling Intra-Diamond method implementation with no communication computation overlapping
/* end */    {0, 0},
};


struct Kernel KernelList[] = {
    // parameters: name string, function pointer, stencil's semi-bandwidth,
    //    stencil's order in time, stencil's shape, coefficients variability

    /* 0 */    {"star", 4, 2, STAR, CONSTANT_COEFFICIENT,
        iso_refi,
        stat_sched_iso_ref,
        swd_iso_ref},
    // basic stencil operator kernel with openMP static,1 schedule, intel's SIMD directive, and cache blocking in the Y direction

    /* 1 */    {"star", 1, 1, STAR, CONSTANT_COEFFICIENT,
        iso_ref_2space_1timei,
        stat_sched_iso_ref_2space_1time,
        swd_iso_ref_2space_1time},
    // basic stencil operator kernel with openMP static,1 schedule, intel's SIMD directive, and cache blocking in the Y direction

    /* 2 */    {"star", 1, 1, STAR, VARIABLE_COEFFICIENT,
        iso_ref_2space_1time_vari,
        stat_sched_iso_ref_2space_1time_var,
        swd_iso_ref_2space_1time_var},
    // basic stencil operator kernel with openMP static,1 schedule, intel's SIMD directive, and cache blocking in the Y direction

    /* 3 */    {"star", 1, 1, STAR, VARIABLE_COEFFICIENT_AXSYM,
        iso_ref_2space_1time_var_axsymi,
        stat_sched_iso_ref_2space_1time_var_axsym,
        swd_iso_ref_2space_1time_var_axsym},
    // basic stencil operator kernel with openMP static,1 schedule, intel's SIMD directive, and cache blocking in the Y direction

    /* 4 */    {"star", 4, 1, STAR, VARIABLE_COEFFICIENT_AXSYM,
        iso_ref_8space_1time_var_axsymi,
        stat_sched_iso_ref_8space_1time_var_axsym,
        swd_iso_ref_8space_1time_var_axsym},
    // basic stencil operator kernel with openMP static,1 schedule, intel's SIMD directive, and cache blocking in the Y direction

    /* 5 */    {"star", 1, 1, STAR, VARIABLE_COEFFICIENT_NOSYM,
        iso_ref_2space_1time_var_nosymi,
        stat_sched_iso_ref_2space_1time_var_nosym,
        swd_iso_ref_2space_1time_var_nosym},
    // basic stencil operator kernel with openMP static,1 schedule, intel's SIMD directive, and cache blocking in the Y direction

    /* end */    {0, 0, 0, 0, 0, 0, 0, 0},
};

#endif /* WRAPPERS_H_ */
