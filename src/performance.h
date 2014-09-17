#include "data_structures.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mpi.h"


extern struct time_stepper TSList[];

#define PRINT_TIME_DETAILS (1)

extern void mpi_halo_init(Parameters *);
extern void print_param(Parameters);
extern void arrays_allocate(Parameters *p);
extern void init_coeff(Parameters *);
extern void domain_data_fill(Parameters * p);
extern void reset_timers(Profile *);
extern void reset_wf_timers(Parameters *);
extern void performance_results(Parameters *p, double t, double t_max, double t_min, double, double);
extern void mpi_halo_finalize(Parameters *);
extern void arrays_free(Parameters *p);
