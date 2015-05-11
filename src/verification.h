#include "data_structures.h"

#define U(i,j,k)         (u[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])
#define V(i,j,k)         (v[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])
#define ROC2(i,j,k)   (roc2[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])
#define TARGET_DOMAIN(i,j,k)   (target_domain[((1ULL)*((k)*(ny)+(j))*(nx)+(i))])
#define SNAPSHOT_ERROR(i,j,k)     (snapshot_error[((1ULL)*((k)*(nny)+(j))*(nnx)+(i))])

#define COEF(m,i,j,k) (coef[((k)*(nny)+(j))*(nnx)+(i)+((ln_domain)*(m))])

void verify(Parameters *);

void aggregate_MPI_subdomains(Parameters vp, real_t * restrict aggr_array);
void verification_printing(Parameters vp);
real_t *restrict aggregate_subdomains(Parameters vp);
void compare_results(real_t *restrict u, real_t *restrict target_domain, int alignment, int nx, int ny, int nz, int NHALO, Parameters p);
void verify_serial_generic(real_t * , Parameters);

void std_kernel_8space_2time( const int shape[3],
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2);
void std_kernel_2space_1time( const int shape[3],
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2);
void std_kernel_2space_1time_var( const int shape[3],
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2);
void std_kernel_2space_1time_var_axsym( const int shape[3],
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2);
void std_kernel_8space_1time_var_axsym( const int shape[3],
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2);
void std_kernel_2space_1time_var_nosym( const int shape[3],
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2);
void solar_kernel( const int shape[3],
    const real_t * restrict coef, real_t * restrict u,
    const real_t * restrict v, const real_t * restrict roc2);

extern void copy_params_struct(Parameters a, Parameters * b);
extern void print_param(Parameters p);
extern void arrays_free(Parameters *);
extern void print_3Darray(char *, real_t * restrict , int , int , int , int );
extern void print_3Darray_solar(char *, real_t * restrict , int , int , int , int );
extern void mpi_halo_finalize(Parameters *);
extern void init_coeff(Parameters *);
extern void arrays_allocate(Parameters *p);
extern void domain_data_fill(Parameters * p);

extern void mpi_halo_init(Parameters *);
extern void check_merr(int e);

extern struct time_stepper TSList[];
