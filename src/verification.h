#include "data_structures.h"

#define U(i,j,k)         (u[((k)*(nny)+(j))*(nnx)+(i)])
#define V(i,j,k)         (v[((k)*(nny)+(j))*(nnx)+(i)])
#define ROC2(i,j,k)   (roc2[((k)*(nny)+(j))*(nnx)+(i)])
#define TARGET_DOMAIN(i,j,k)   (target_domain[((k)*(ny)+(j))*(nx)+(i)])
#define SNAPSHOT_ERROR(i,j,k)     (snapshot_error[((k)*(nny)+(j))*(nnx)+(i)])

#define COEF(m,i,j,k) (coef[((k)*(nny)+(j))*(nnx)+(i)+(((unsigned long) (ln_domain))*(m))])

void verify(Parameters *);

void aggregate_MPI_subdomains(Parameters vp, FLOAT_PRECISION * restrict aggr_array);
void verification_printing(Parameters vp);
FLOAT_PRECISION *restrict aggregate_subdomains(Parameters vp);
void compare_results(FLOAT_PRECISION *restrict u, FLOAT_PRECISION *restrict target_domain, int alignment, int nx, int ny, int nz, int NHALO, Parameters p);
void verify_serial_generic(FLOAT_PRECISION * , Parameters);

void std_kernel_8space_2time( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2);
void std_kernel_2space_1time( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2);
void std_kernel_2space_1time_var( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2);
void std_kernel_2space_1time_var_axsym( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2);
void std_kernel_8space_1time_var_axsym( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2);
void std_kernel_2space_1time_var_nosym( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2);
void solar_kernel( const int shape[3],
    const FLOAT_PRECISION * restrict coef, FLOAT_PRECISION * restrict u,
    const FLOAT_PRECISION * restrict v, const FLOAT_PRECISION * restrict roc2);

extern void copy_params_struct(Parameters a, Parameters * b);
extern void print_param(Parameters p);
extern void arrays_free(Parameters *);
extern void print_3Darray(char *, FLOAT_PRECISION * restrict , int , int , int , int );
extern void print_3Darray_solar(char *, FLOAT_PRECISION * restrict , int , int , int , int );
extern void mpi_halo_finalize(Parameters *);
extern void init_coeff(Parameters *);
extern void arrays_allocate(Parameters *p);
extern void domain_data_fill(Parameters * p);

extern void mpi_halo_init(Parameters *);
extern void check_merr(int e);

extern struct time_stepper TSList[];
