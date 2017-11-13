/** 
 * @copyright (c) 2014- King Abdullah University of Science and
 *                      Technology (KAUST). All rights reserved.
 **/
 

/**
 * @file src/utils.h 

 * GIRIH is a high performance stencil framework using wavefront 
 * 	diamond tiling.
 * GIRIH is provided by KAUST.
 *
 * @version 1.0.0
 * @author Tareq Malas
 * @date 2017-11-13
 **/

#ifndef UTILS_H_
#define UTILS_H_

#include "data_structures.h"

void init(Parameters *);
void set_centered_source(Parameters *);
void init_coeff(Parameters *);
void arrays_allocate(Parameters *p);
void reset_timers(Profile * p);
void print_3Darray(char *, real_t * restrict , int , int , int , int );
void copy_params_struct(Parameters a, Parameters * b);
void domain_data_fill(Parameters * p);
void mpi_halo_finalize(Parameters *p);
void arrays_free(Parameters *p);

#include <stdio.h>
//extern void domain_decompose(Parameters *p);
extern void mpi_halo_init(Parameters *);
extern void domain_decompose(Parameters *p);

extern int get_ntg(Parameters);
extern struct time_stepper TSList[];
extern const char *MWD_name[];

extern void dynamic_intra_diamond_ts(Parameters *p);
extern void intra_diamond_info_init(Parameters *p);

extern clu_func_t clu_func_list[];
extern struct StencilInfo stencil_info_list[];
extern spt_blk_func_t spt_blk_func_list[];
extern spt_blk_func_t stat_sched_func_list[];
extern mwd_func_t swd_func_list[];

extern mwd_func_t mwd_func_list[];
extern mwd_func_t femwd_func_list[];
extern mwd_func_t *mwd_list[];

extern void swd_iso_ref_split KERNEL_MWD_SIG;
extern void mwd_iso_ref_split KERNEL_MWD_SIG;
extern void iso_ref_split KERNEL_SIG;


#endif /* UTILS_H_ */
