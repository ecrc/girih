#ifndef UTILS_H_
#define UTILS_H_

#include "data_structures.h"

void init(Parameters *);
void set_centered_source(Parameters *);
void init_coeff(Parameters *);
void arrays_allocate(Parameters *p);
void reset_timers(Profile * p);
void print_3Darray(char *, FLOAT_PRECISION * restrict , int , int , int , int );
void copy_params_struct(Parameters a, Parameters * b);
void domain_data_fill(Parameters * p);
void mpi_halo_finalize(Parameters *p);


//extern void domain_decompose(Parameters *p);
extern void mpi_halo_init(Parameters *);

// 1-stride of ISO stencils
extern void iso_ref_8space_2time_stride STRIDE1_SIG;
extern void iso_ref_2space_1time_stride STRIDE1_SIG;
extern void iso_ref_2space_1time_var_stride STRIDE1_SIG;
extern void iso_ref_2space_1time_var_axsym_stride STRIDE1_SIG;
extern void domain_decompose(Parameters *p);


extern struct time_stepper TSList[];
extern struct Kernel KernelList[];

extern void intra_diamond_1wf_comp(Parameters *p, int yb_r, int ye_r, int b_inc, int e_inc);
extern void intra_diamond_all_mwf_comp(Parameters *p, int yb_r, int ye_r, int b_inc, int e_inc, int tid);
extern void dynamic_intra_diamond_ts(Parameters *p);

extern mwd_func_t mwd_func_list[];
extern mwd_func_t femwd_func_list[];
extern void iso_ref_all_wf_split KERNEL_MWD_SIG;


#endif /* UTILS_H_ */