#include "stencils.h"

/*
 * ISO stencil 8th-order-in-space-2nd-order-in-time with constant coefficient
 */
#ifdef FUNC_BODY
#undef FUNC_BODY
#endif
#define FUNC_BODY() { \
ux[i] = ((real_t) (2.0))*vx[i] - ux[i] + ROC2(i,j,k)*(coef[0]*vx[i] \
+coef[1]*(vx[i+1]+vx[i-1]) \
+coef[1]*(vx[i+nnx]+vx[i-nnx]) \
+coef[1]*(vx[+nnxy+i]+vx[-nnxy+i]) \
+coef[2]*(vx[i+2]+vx[i-2]) \
+coef[2]*(vx[i+2*nnx]+vx[i-2*nnx]) \
+coef[2]*(vx[+2*nnxy+i]+vx[-2*nnxy+i]) \
+coef[3]*(vx[i+3]+vx[i-3]) \
+coef[3]*(vx[i+3*nnx]+vx[i-3*nnx]) \
+coef[3]*(vx[+3*nnxy+i]+vx[-3*nnxy+i]) \
+coef[4]*(vx[i+4]+vx[i-4]) \
+coef[4]*(vx[i+4*nnx]+vx[i-4*nnx]) \
+coef[4]*(vx[+4*nnxy+i]+vx[-4*nnxy+i]) ); \
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME iso_ref

#include "stencils_list.h"


/*
 * ISO stencil 2nd-order-in-space-1st-order-in-time with constant coefficient
 */
#ifdef FUNC_BODY
#undef FUNC_BODY
#endif
#define FUNC_BODY() { \
ux[i] = coef[0]*vx[i] \
+coef[1]*(vx[i+1]+vx[i-1]) \
+coef[1]*(vx[i+nnx]+vx[i-nnx]) \
+coef[1]*(vx[-nnxy+i]+vx[+nnxy+i]); \
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME iso_ref_2space_1time

#include "stencils_list.h"


/*
 * ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient
 */
#ifdef FUNC_BODY
#undef FUNC_BODY
#endif
#define FUNC_BODY() { \
ux[i] = COEF(0,i,j,k)*vx[i] \
+COEF(1,i,j,k)*(vx[i+1]+vx[i-1]) \
+COEF(1,i,j,k)*(vx[i+nnx]+vx[i-nnx]) \
+COEF(1,i,j,k)*(vx[+nnxy+i]+vx[-nnxy+i]); \
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME iso_ref_2space_1time_var

#include "stencils_list.h"


/*
 * ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient axis symmetry
 */
#ifdef FUNC_BODY
#undef FUNC_BODY
#endif
#define FUNC_BODY() { \
ux[i] = COEF(0,i,j,k)*vx[i] \
+COEF(1,i,j,k)*(vx[i+1]+vx[i-1]) \
+COEF(2,i,j,k)*(vx[i+nnx]+vx[i-nnx]) \
+COEF(3,i,j,k)*(vx[+nnxy+i]+vx[-nnxy+i]); \
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME iso_ref_2space_1time_var_axsym

#include "stencils_list.h"


/*
 * ISO stencil 8th-order-in-space-1st-order-in-time with variable axis symmetric coefficient
 */
#ifdef FUNC_BODY
#undef FUNC_BODY
#endif
#define FUNC_BODY() { \
ux[i] = COEF(0 ,i,j,k)*vx[i] \
+COEF(1 ,i,j,k)*(vx[i+1]+vx[i-1]) \
+COEF(2 ,i,j,k)*(vx[i+nnx]+vx[i-nnx]) \
+COEF(3 ,i,j,k)*(vx[+nnxy+i]+vx[-nnxy+i]) \
+COEF(4 ,i,j,k)*(vx[i+2]+vx[i-2]) \
+COEF(5 ,i,j,k)*(vx[i+2*nnx]+vx[i-2*nnx]) \
+COEF(6 ,i,j,k)*(vx[+2*nnxy+i]+vx[-2*nnxy+i]) \
+COEF(7 ,i,j,k)*(vx[i+3]+vx[i-3]) \
+COEF(8 ,i,j,k)*(vx[i+3*nnx]+vx[i-3*nnx]) \
+COEF(9 ,i,j,k)*(vx[+3*nnxy+i]+vx[-3*nnxy+i]) \
+COEF(10,i,j,k)*(vx[i+4]+vx[i-4]) \
+COEF(11,i,j,k)*(vx[i+4*nnx]+vx[i-4*nnx]) \
+COEF(12,i,j,k)*(vx[+4*nnxy+i]+vx[-4*nnxy+i]); \
}

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME iso_ref_8space_1time_var_axsym

#include "stencils_list.h"


/*
 * ISO stencil 2nd-order-in-space-1st-order-in-time with variable coefficient no symmetry
 */
#ifdef FUNC_BODY
#undef FUNC_BODY
#endif
#define FUNC_BODY() { \
ux[i] = COEF(0,i,j,k)*vx[i] \
+COEF(1,i,j,k)*vx[i-1] \
+COEF(2,i,j,k)*vx[i+1] \
+COEF(3,i,j,k)*vx[i-nnx] \
+COEF(4,i,j,k)*vx[i+nnx] \
+COEF(5,i,j,k)*vx[-nnxy+i] \
+COEF(6,i,j,k)*vx[+nnxy+i]; \
  }

#ifdef FUNC_NAME
#undef FUNC_NAME
#endif
#define FUNC_NAME iso_ref_2space_1time_var_nosym

#include "stencils_list.h"



#include "solar_kernels.h"

/*
 * Wrap the created functions in array
 */

struct StencilInfo stencil_info_list[] = {
    // parameters: name, semi-bandwidth, order in time, domain replicas, shape, coefficients variability
    {"star", 4, 2, 3,  STAR, CONSTANT_COEFFICIENT, REGULAR},
    {"star", 1, 1, 2,  STAR, CONSTANT_COEFFICIENT, REGULAR},
    {"star", 1, 1, 4,  STAR, VARIABLE_COEFFICIENT, REGULAR},
    {"star", 1, 1, 6,  STAR, VARIABLE_COEFFICIENT_AXSYM, REGULAR},
    {"star", 4, 1, 15, STAR, VARIABLE_COEFFICIENT_AXSYM, REGULAR},
    {"star", 1, 1, 9 , STAR, VARIABLE_COEFFICIENT_NOSYM, REGULAR},
    {"star", 1, 1, 40, STAR, SOLAR_COEFFICIENT, SOLAR},
    {0, 0, 0, 0, 0, 0, 0},
};

spt_blk_func_t spt_blk_func_list[] = {
    iso_ref,
    iso_ref_2space_1time,
    iso_ref_2space_1time_var,
    iso_ref_2space_1time_var_axsym,
    iso_ref_8space_1time_var_axsym,
    iso_ref_2space_1time_var_nosym,
    solar,};

spt_blk_func_t stat_sched_func_list[] = {
    stat_sched_iso_ref,
    stat_sched_iso_ref_2space_1time,
    stat_sched_iso_ref_2space_1time_var,
    stat_sched_iso_ref_2space_1time_var_axsym,
    stat_sched_iso_ref_8space_1time_var_axsym,
    stat_sched_iso_ref_2space_1time_var_nosym,
    stat_sched_solar,};

mwd_func_t swd_func_list[] = {
    swd_iso_ref,
    swd_iso_ref_2space_1time,
    swd_iso_ref_2space_1time_var,
    swd_iso_ref_2space_1time_var_axsym,
    swd_iso_ref_8space_1time_var_axsym,
    swd_iso_ref_2space_1time_var_nosym,
    swd_solar,};


mwd_func_t mwd_func_list[] = {  /* 0 */
    mwd_iso_ref,
    mwd_iso_ref_2space_1time,
    mwd_iso_ref_2space_1time_var,
    mwd_iso_ref_2space_1time_var_axsym,
    mwd_iso_ref_8space_1time_var_axsym,
    mwd_iso_ref_2space_1time_var_nosym,
    mwd_solar,};

mwd_func_t femwd_func_list[] = { /* 1 */
    femwd_iso_ref,
    femwd_iso_ref_2space_1time,
    femwd_iso_ref_2space_1time_var,
    femwd_iso_ref_2space_1time_var_axsym,
    femwd_iso_ref_8space_1time_var_axsym,
    femwd_iso_ref_2space_1time_var_nosym,
    femwd_solar,};

mwd_func_t rsmwd_func_list[] = { /* 2 */
    rsmwd_iso_ref,
    rsmwd_iso_ref_2space_1time,
    rsmwd_iso_ref_2space_1time_var,
    rsmwd_iso_ref_2space_1time_var_axsym,
    rsmwd_iso_ref_8space_1time_var_axsym,
    rsmwd_iso_ref_2space_1time_var_nosym,
    rsmwd_solar,};

mwd_func_t rsfemwd_func_list[] = { /* 3 */
    rsfemwd_iso_ref,
    rsfemwd_iso_ref_2space_1time,
    rsfemwd_iso_ref_2space_1time_var,
    rsfemwd_iso_ref_2space_1time_var_axsym,
    rsfemwd_iso_ref_8space_1time_var_axsym,
    rsfemwd_iso_ref_2space_1time_var_nosym,
    not_supported_mwd,};



const char *MWD_name[] = {"Wavefront", 
                           "Fixed execution wavefronts", 
                           "Relaxed synchronization wavefront", 
                           "Relaxed synchronization wavefront with fixed execution", 
                           0,};

mwd_func_t *mwd_list[] = {mwd_func_list, femwd_func_list, rsmwd_func_list, rsfemwd_func_list, };


