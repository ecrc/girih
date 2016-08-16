#include "performance.h"

// from GNU.org
int compare_doubles (const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}

void performance_test(Parameters * p){
  int i;
  int tests_remain;
  double t_start, t_end, t=0.0;
  double t_min = 1000000.0, t_max= -1.0, t_med, tpercycle, tpertest;
  double  t_ts_main_max= -1.0, t_ts_main_min = 1000000.0;
  double t2;
  double *ttests = (double*) malloc(p->n_tests*sizeof(double));

  // Create the required MPI data structures
  mpi_halo_init(p);

  time_t now;
  if (p->mpi_rank == 0) {
    time(&now);
    printf("Started on %s", ctime(&now));
    // print the used parameters in the experiment if required
    if(p->verbose ==1) print_param(*p);
  }

  // allocate and initialize the required arrays for the performance experiments
  arrays_allocate(p);
  init_coeff(p);
  domain_data_fill(p);
  if(p->source_point_enabled == 1)
    for(i=0; i<p->nt; i++) p->source[i] = (real_t) i;

  // run the performance experiments of the target kernel
  if (p->mpi_rank == 0) {
    printf("\n******************************************************\n");
    printf("Performance results\n");
    printf("******************************************************\n");
  }
  tests_remain = p->n_tests;

  if(p->target_ts != 2) {
    #pragma omp parallel
    {
      MARKER_START("calc");
    }
  }
  while(tests_remain--) {
    reset_timers(&(p->prof));
    reset_wf_timers(p);
    MPI_Barrier(MPI_COMM_WORLD);
    t_start = MPI_Wtime();
    TSList[p->target_ts].func(p);
    t2 = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    t_end = MPI_Wtime();
    tpertest = t_end - t_start;
    tpercycle = (t_end - t_start) / p->nt;
    ttests[tests_remain] = tpertest;

    p->prof.total = t_end - t_start;
    p->prof.wait += (t_end - t2);
    p->prof.communicate += p->prof.wait;
    p->prof.others = p->prof.total - p->prof.communicate - p->prof.compute;

#if PRINT_TIME_DETAILS
    if (p->mpi_rank == 0) {
      if (tests_remain > 0)
        printf("Rank 0 TEST#%02d time: %e\n",(p->n_tests - tests_remain),tpertest);
      else if (tests_remain == 0) {
        printf("Rank 0 TEST#%02d time: %e\n",(p->n_tests - tests_remain),tpertest);
        printf("******************************************************\n");
      }
    }
#endif
    t += tpertest/p->n_tests;
    if (t_min > tpercycle) t_min = tpercycle;
    if (t_max < tpercycle) t_max = tpercycle;
    if (t_ts_main_min > p->prof.ts_main) t_ts_main_min = p->prof.ts_main;
    if (t_ts_main_max < p->prof.ts_main) t_ts_main_max = p->prof.ts_main;
  }
  if(p->target_ts != 2) {
    #pragma omp parallel
    {
      MARKER_STOP("calc");
    }
  }

  // compute the tests median
  qsort(ttests, p->n_tests, sizeof(double), compare_doubles);  
  t_med = ttests[p->n_tests/2];

  // collect and print the performance results
  performance_results(p, t, t_max, t_min, t_med, t_ts_main_max, t_ts_main_min);

  // clean up
  if (p->mpi_rank == 0) {
    time(&now);
    printf("COMPLETED SUCCESSFULLY on %s", ctime(&now));
  }
  mpi_halo_finalize(p);
  arrays_free(p);
  free(ttests);

}
