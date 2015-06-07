#include "data_structures.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


typedef struct Tune_Params{
  int thread_group_size, th_z, th_y, th_x, diamond_width, num_wavefronts;
}Tune_Params;

int get_ntg(Parameters p){
  return (int) ceil(1.0*p.num_threads/p.stencil_ctx.thread_group_size);
}


void cpu_bind_init(Parameters *p){

  // Source for finding number of CPUs: https://software.intel.com/en-us/blogs/2013/10/31/applying-intel-threading-building-blocks-observers-for-thread-affinity-on-intel
  cpu_set_t *mask;
  int ncpus;
  for ( ncpus = sizeof(cpu_set_t)/8; ncpus < 16*1024; ncpus <<= 1 ) {
    mask = CPU_ALLOC( ncpus );
    if ( !mask ) break;
    const size_t size = CPU_ALLOC_SIZE( ncpus );
    CPU_ZERO_S( size, mask );
    const int err = sched_getaffinity( 0, size, mask );
    if ( !err ) break;
    CPU_FREE( mask );
    mask = NULL;
    if ( errno != EINVAL )  break;
  }
  if ( !mask )
    printf("Warning: Failed to obtain process affinity mask. Thread affinitization is disabled.\n");

  p->stencil_ctx.setsize = CPU_ALLOC_SIZE(ncpus);
  p->stencil_ctx.bind_masks = (cpu_set_t**) malloc(p->num_threads*sizeof(cpu_set_t*));

  int i, ib, idx=0;
  ib=0;
#if __MIC__
  ib = 1;
#endif
  for(i=ib; i<p->num_threads*p->th_stride/p->th_block+ib;i++){
    if((i-ib)%p->th_stride < p->th_block){
      p->stencil_ctx.bind_masks[idx] = CPU_ALLOC( ncpus );
      CPU_ZERO_S(p->stencil_ctx.setsize, p->stencil_ctx.bind_masks[idx]);
      CPU_SET_S(i,p->stencil_ctx.setsize, p->stencil_ctx.bind_masks[idx]);
      idx++;
    }
  }

  int max_threads = omp_get_max_threads();
  int *phys_cpu = (int*) malloc(max_threads*sizeof(int));

  omp_set_nested(1);
  // Set the affinity to reduce the cost of first run
  int num_thread_groups = get_ntg(*p);
#pragma omp parallel num_threads(num_thread_groups) PROC_BIND(spread)
  {
    int mtid = omp_get_thread_num();
  #pragma omp parallel shared(mtid)  num_threads(p->stencil_ctx.thread_group_size) PROC_BIND(master)
    {
      int tid = omp_get_thread_num();
      int gtid = tid + mtid * p->stencil_ctx.thread_group_size;

      int err = sched_setaffinity(0, p->stencil_ctx.setsize, p->stencil_ctx.bind_masks[gtid]);
      if(err==-1) printf("WARNING: Could not set CPU Affinity of thread:%d error:%d\n", gtid, err);
      phys_cpu[gtid] = sched_getcpu();
    }
  }

  printf("Using OS threads binding (OS_tid->Phys_id):");
  for(i=0;i<p->num_threads;i++){
    printf("  %d->%d", i, phys_cpu[i]);
  }
  printf("\n");

  free(phys_cpu);
}


uint64_t get_mwf_size(Parameters p, int t_dim){
  uint64_t diam_width, diam_height, wf_updates, wf_elements, lnx, t_order, total_points;

  t_order = p.stencil.time_order;
  diam_width = (t_dim+1)*2*p.stencil.r;
  int nwf = p.stencil_ctx.num_wf;
  diam_height = t_dim*2*p.stencil.r + nwf;

  lnx = p.ldomain_shape[0];

  wf_updates = (t_dim+1)*(t_dim+1)*2 * p.stencil.r; // Y-T projection
  wf_elements = (wf_updates - diam_width) * p.stencil.r + diam_width + diam_width*(nwf-1);

  switch(p.stencil.coeff){
  case CONSTANT_COEFFICIENT:
    total_points = ( (t_order+1)             *wf_elements + (diam_width + diam_height )*2*p.stencil.r) * lnx * sizeof(real_t);
    break;

  case VARIABLE_COEFFICIENT:
    total_points = ( (t_order+1 + (1+p.stencil.r) )*wf_elements + (diam_width + diam_height )*2*p.stencil.r) * lnx * sizeof(real_t);
    break;

  case VARIABLE_COEFFICIENT_AXSYM:
    total_points = ( (t_order+1 + (1+3*p.stencil.r) )*wf_elements + (diam_width + diam_height )*2*p.stencil.r) * lnx * sizeof(real_t);
    break;

  case VARIABLE_COEFFICIENT_NOSYM:
    total_points = ( (t_order+1 + (1+6*p.stencil.r) )*wf_elements + (diam_width + diam_height )*2*p.stencil.r) * lnx * sizeof(real_t);
    break;

  case SOLAR_COEFFICIENT:
    total_points = ( (p.stencil.nd)*wf_elements + (diam_width + diam_height )*12*p.stencil.r) * lnx * 2*sizeof(real_t);
    break;


  default:
    printf("ERROR: unknown type of stencil\n");
    exit(1);
    break;
  }
  //printf("npx:%d  nx:%d  lnx:%lu  updates:%lu  elements:%lu  total:%lu\n", p.t.shape[0],  p.ldomain_shape[0], lnx, wf_updates, wf_elements, total_points);
  return total_points;
}
double run_tuning_test(Parameters *tp){
  double t = 0.0;
  int reps = 1, orig_reps;
  float obt_perf;

  do{
    tp->nt = reps*(tp->t_dim+1)*2 + 2;
    tp->prof.ts_main = 0;
    dynamic_intra_diamond_ts(tp);
    t = tp->prof.ts_main;
    orig_reps =reps;
    reps *= (int) ceil(3.0/t);
  } while(t < 2.0);

  uint64_t lups = (tp->ln_stencils*tp->nt - tp->idiamond_pro_epi_logue_updates);
  obt_perf =  lups/tp->prof.ts_main;

  printf("%06.2f]  time:%es  lups:%lu  cache block size:%lukiB  reps:%d\n",
      obt_perf/(1e6), tp->prof.ts_main, lups, get_mwf_size(*tp, tp->t_dim)*get_ntg(*tp)/1024, orig_reps);

  return obt_perf;
}

void get_feasible_diam_widthes(Parameters p, int *diam_widthes, int *num_diam_widthes){
  int n_diam_widthes, max_t_dim, i, y_len, lt_dim, diam_width, wf_len, ntg, wf_size, num_wf, idx;
  int cache_size_cond, int_diam_cond, wf_len_cond, cuncurrency_cond, diam_concurrency;

  y_len = p.stencil_shape[1]/p.t.shape[1];
  num_wf = (p.stencil_ctx.num_wf>0? p.stencil_ctx.num_wf: 1);
  ntg = get_ntg(p);
  n_diam_widthes=0;
  for(i=y_len/p.stencil.r; i>=4; i-=4){
    lt_dim = i/2 - 1;
    diam_width = i*p.stencil.r;
    wf_len = lt_dim*2*p.stencil.r+1 + num_wf-1;

    wf_size = get_mwf_size(p, lt_dim);
    cache_size_cond = wf_size*ntg < (1ULL*MAX_CACHE_SIZE*1024);

    cuncurrency_cond = (y_len/diam_width) >= ntg;
    int_diam_cond = y_len%diam_width == 0;
    wf_len_cond = wf_len <= p.stencil_shape[2];

//    printf("i:%d, diam_width %d,  cuncurrency_cond %d, cache_size_cond %d, int_diam_cond %d, wf_len_cond %d, cache_blk_size: %lu kB\n",
//        i, diam_width, cuncurrency_cond, cache_size_cond, int_diam_cond, wf_len_cond, wf_size*ntg/1024);
    if( (int_diam_cond == 1) && (wf_len_cond == 1) && (cuncurrency_cond == 1)  && (cache_size_cond == 1) ){ // consider limitation in z and concurrency
      n_diam_widthes++;
    }
  }

  // check if no feasible diam width is found
  if(n_diam_widthes == 0){
    diam_widthes=NULL;
    *num_diam_widthes=0;
    return;
  }

  // allocate and fill the diam widthes array
  diam_widthes = (int*) malloc(n_diam_widthes*sizeof(int));
  idx=0;
  for(i=y_len/p.stencil.r; i>=4; i-=4){
    lt_dim = i/2 - 1;
    diam_width = i*p.stencil.r;
    wf_len = lt_dim*2*p.stencil.r+1 + num_wf-1;

    wf_size = get_mwf_size(p, lt_dim);
    cache_size_cond = wf_size*ntg < (1ULL*MAX_CACHE_SIZE*1024);

    cuncurrency_cond = (y_len/diam_width) >= ntg;
    int_diam_cond = y_len%diam_width == 0;
    wf_len_cond = wf_len <= p.stencil_shape[2];
    if( (int_diam_cond == 1) && (wf_len_cond == 1) && (cuncurrency_cond == 1)  && (cache_size_cond == 1) ){ // consider limitation in z and concurrency
      diam_widthes[idx++] = lt_dim;
    }
  }
  *num_diam_widthes = n_diam_widthes;

}

void auto_tune_diam_nwf(Parameters *op){
  Parameters tp;
  copy_params_struct(*op, &tp);
  tp.stencil_shape[0] = tp.stencil_shape[0]/tp.t.shape[0];
  tp.stencil_shape[1] = tp.stencil_shape[1]/tp.t.shape[1];
  tp.stencil_shape[2] = tp.stencil_shape[2]/tp.t.shape[2];
  tp.t.shape[0]=1;
  tp.t.shape[1]=1;
  tp.t.shape[2]=1;
  tp.mpi_size = 1;

  if( (tp.wavefront == 1) && (tp.stencil_ctx.thread_group_size != 1) ) // multi-thread group
    tp.wavefront = -1;

  double exp_perf, cur_perf, prev_nwf_perf, prev_diam_perf;
  int i, lt_dim, max_t_dim, prev_max_nwf, diam_width, prev_t_dim;
  uint64_t wf_size, ntg;
  int cache_size_cond, int_diam_cond, wf_len_cond, cuncurrency_cond, diam_concurrency;
  int diam_height;
  int thz = tp.stencil_ctx.th_z;
  ntg = get_ntg(tp);

  tp.stencil_ctx.num_wf = thz; // set the number of wavefronts to the minimum possible value
  if(tp.mwd_type == 3) tp.stencil_ctx.num_wf = thz*tp.stencil.r;


  // find max possible diamond width and allocate memory accordingly, if not pre set
  if(op->t_dim == -1){
    max_t_dim = -1;
    for(i=tp.stencil_shape[1]/tp.stencil.r; i>=4; i-=4){
      lt_dim = i/2 - 1;
      diam_width = i*tp.stencil.r;
      diam_height = lt_dim*2*tp.stencil.r+1 + tp.stencil_ctx.num_wf-1;

      wf_size = get_mwf_size(tp, lt_dim);
      cache_size_cond = wf_size*ntg < (uint64_t) (MAX_CACHE_SIZE*1024);

      cuncurrency_cond = (tp.stencil_shape[1]/diam_width) >= ntg;
      int_diam_cond = tp.stencil_shape[1]%diam_width == 0;
      wf_len_cond = diam_height <= tp.stencil_shape[2];

  //    printf("i:%d, diam_width %d,  cuncurrency_cond %d, cache_size_cond %d, int_diam_cond %d, wf_len_cond %d, cache_blk_size: %lu kB\n",
  //        i, diam_width, cuncurrency_cond, cache_size_cond, int_diam_cond, wf_len_cond, wf_size*ntg/1024);
      if( (int_diam_cond == 1) && (wf_len_cond == 1) && (cuncurrency_cond == 1)  && (cache_size_cond == 1) ){ // consider limitation in z and concurrency
        tp.t_dim = lt_dim;
        max_t_dim = lt_dim;
  //      tp.stencil_shape[1] = i*tp.stencil.r;
        break;
      }
    }

    if (max_t_dim == -1){
      op->t_dim = -1;
      return;
    }

    printf("[AUTO TUNE] max. diamond width: %d\n", (max_t_dim+1)*2*tp.stencil.r);
  }

  // initialize the data of the tuning experiments
  init(&tp);
  arrays_allocate(&tp);
  init_coeff(&tp);
  domain_data_fill(&tp);

  // Allocate the wavefront profiling timers
  int num_thread_groups = get_ntg(tp);
  tp.stencil_ctx.t_wait        = (double *) malloc(sizeof(double)*tp.num_threads);
  tp.stencil_ctx.t_wf_main     = (double *) malloc(sizeof(double)*num_thread_groups);
  tp.stencil_ctx.t_wf_comm = (double *) malloc(sizeof(double)*num_thread_groups);
  tp.stencil_ctx.t_wf_prologue = (double *) malloc(sizeof(double)*num_thread_groups);
  tp.stencil_ctx.t_wf_epilogue = (double *) malloc(sizeof(double)*num_thread_groups);
  tp.stencil_ctx.wf_num_resolved_diamonds = (double *) malloc(sizeof(double)*num_thread_groups);
  tp.stencil_ctx.t_group_wait = (double *) malloc(sizeof(double)*num_thread_groups);

  if(op->t_dim == -1){ // tune both diamond widht and number of frontlines
    // brute force search for best diamond/nwf starting from small to large
    // test diamond sizes from smallest to largest
    prev_diam_perf = -1;
    prev_t_dim = 0;
    for(i=4; i<=(max_t_dim+1)*2; i+=4){ // loop over diamond sizes
      tp.t_dim = i/2 - 1;
      diam_width = i*tp.stencil.r;
      wf_size = get_mwf_size(tp, tp.t_dim);
      diam_concurrency = tp.stencil_shape[1]/diam_width;

      cache_size_cond = wf_size*ntg > (uint64_t) (tp.cache_size*1024);
      cuncurrency_cond = diam_concurrency >= ntg;
      int_diam_cond = tp.stencil_shape[1]%diam_width == 0;
  //    printf("i:%d, diam_width %d,  cuncurrency_cond %d, cache_size_cond %d, int_diam_cond %d, wf_len_cond %d, cache_blk_size: %lu kB\n",
  //        i, diam_width, cuncurrency_cond, cache_size_cond, int_diam_cond, wf_len_cond, wf_size*ntg/1024);

      if( (int_diam_cond == 1) && (cuncurrency_cond == 1) && (cache_size_cond == 1) ){ // check diamond size validity
        printf("[AUTO TUNE] Diamond width:%02d  [wavefronts #: pefromance (MLUPS/s)]\n", (tp.t_dim+1)*2*tp.stencil.r);
        // loop over increasing number of wavefronts per update
        prev_nwf_perf = -1;
        tp.stencil_ctx.num_wf = thz; // start with smallest possible number of updates
        if(tp.mwd_type == 3) tp.stencil_ctx.num_wf = thz*tp.stencil.r;
        tp.idiamond_pro_epi_logue_updates = (uint64_t) (tp.stencil_shape[0] * tp.stencil_shape[2]) * (uint64_t) (2*diam_concurrency) * ((tp.t_dim+1)*(tp.t_dim+1) + (tp.t_dim+1))*tp.stencil.r;

        while(1){
          wf_size = get_mwf_size(tp, tp.t_dim);
          cache_size_cond = wf_size*ntg > (uint64_t) (tp.cache_size*1024);
          diam_height = tp.t_dim*2*tp.stencil.r+1 +tp.stencil_ctx.num_wf-1;
          wf_len_cond = diam_height <= tp.stencil_shape[2];

          if( (wf_len_cond==1) && (cache_size_cond == 1) ){
            printf("[AUTO TUNE]     [%03d: ",tp.stencil_ctx.num_wf);
            exp_perf = run_tuning_test(&tp);

  //          printf("tgs:%d  nwf:%d  perf:%6.2f  prev_nwf_perf:%6.2f\n", tgs, tp.stencil_ctx.num_wf, exp_perf/1024/1024, prev_nwf_perf/1024/1024);

            // termination criteria for the nwf
            if (exp_perf < prev_nwf_perf){
              tp.stencil_ctx.num_wf -= thz;
              break;
            }
            else{
              prev_nwf_perf = exp_perf;
              tp.stencil_ctx.num_wf += thz;
            }

          }
          else{ // invalid wavefront length
            tp.stencil_ctx.num_wf -= thz;
            break;
          }

        }
        if(tp.stencil_ctx.num_wf < thz){
          tp.stencil_ctx.num_wf = thz;
          //printf("ERROR: Invalid Wavefronts #\n");
          //exit(1);
        }

        // termination criteria for diamond size
        if (prev_nwf_perf < prev_diam_perf){
          tp.t_dim = prev_t_dim; // revert to previous diamond size
          tp.stencil_ctx.num_wf = prev_max_nwf;
          break;
        }
        else{
          prev_diam_perf = prev_nwf_perf;
          prev_t_dim = tp.t_dim;
          prev_max_nwf = tp.stencil_ctx.num_wf;
        }

      }
    }
    if (op->t_dim < 1) op->t_dim = 1;

    op->t_dim = tp.t_dim;
    op->stencil_ctx.num_wf = tp.stencil_ctx.num_wf;


  } else { // diamond width is provided but not the number of frontlines

    diam_width =  ((tp.t_dim+1)*2)*tp.stencil.r;
    diam_concurrency = tp.stencil_shape[1]/diam_width;
    prev_nwf_perf = -1;
    tp.idiamond_pro_epi_logue_updates = (uint64_t) (tp.stencil_shape[0] * tp.stencil_shape[2]) * (uint64_t) (2*diam_concurrency) * ((tp.t_dim+1)*(tp.t_dim+1) + (tp.t_dim+1))*tp.stencil.r;

    while(1){
      diam_height = tp.t_dim*2*tp.stencil.r+1 +tp.stencil_ctx.num_wf-1;
      wf_len_cond = diam_height <= tp.stencil_shape[2];

      if(wf_len_cond==1){
        printf("[AUTO TUNE]     [%03d: ",tp.stencil_ctx.num_wf);
        exp_perf = run_tuning_test(&tp);

//          printf("tgs:%d  nwf:%d  perf:%6.2f  prev_nwf_perf:%6.2f\n", tgs, tp.stencil_ctx.num_wf, exp_perf/1024/1024, prev_nwf_perf/1024/1024);

        // termination criteria for the nwf
        if (exp_perf < prev_nwf_perf){
          tp.stencil_ctx.num_wf -= thz;
          break;
        }
        else{
          prev_nwf_perf = exp_perf;
          tp.stencil_ctx.num_wf += thz;
        }

      }
      else{ // invalid wavefront length
        tp.stencil_ctx.num_wf -= thz;
        break;
      }

    }
    if(tp.stencil_ctx.num_wf < thz){
      tp.stencil_ctx.num_wf = thz;
    }
    if(tp.mwd_type == 3){
      if(tp.stencil_ctx.num_wf < thz*tp.stencil.r){
        tp.stencil_ctx.num_wf = thz*tp.stencil.r;
      }
    }
    op->stencil_ctx.num_wf = tp.stencil_ctx.num_wf;
  }

  free(tp.stencil_ctx.t_wait);
  free(tp.stencil_ctx.t_wf_main);
  free(tp.stencil_ctx.t_wf_comm);
  free(tp.stencil_ctx.t_wf_prologue);
  free(tp.stencil_ctx.t_wf_epilogue);
  free(tp.stencil_ctx.wf_num_resolved_diamonds);
  free(tp.stencil_ctx.t_group_wait);
  arrays_free(&tp);
  mpi_halo_finalize(&tp);
}

int* get_num_factors(int val, int *n_factors){
  int *fact_l;
  int n_fact=0, i;
  for(i=1; i<=val; i++){
    if(val%i==0) n_fact++;
  }
  fact_l = (int*) malloc(n_fact*sizeof(int));
  int idx=0;
  for(i=1; i<=val; i++){
    if(val%i==0){
      fact_l[idx] = i;
      idx++;
    }
  }
  *n_factors = n_fact;
  return fact_l;
}
void get_tgs_tune_params_lists(Parameters *p, Tune_Params **ret_tune_cases_l, int *num_tune_cases, int **ret_tgs_l, int *num_tgs){
  int i, x, y, z, max_tgs, tgsi, n_tgs, n_thz, n_thy, n_thx, n_tgs_factors, idx, n_tune_cases;
  int *tgs_l, *thz_l, *thy_l, *thx_l, *tgs_factors_l;
  Tune_Params *tune_cases_l;

  max_tgs = p->num_threads;
  if(p->num_threads > MAX_THREAD_GROUP_SIZE){
    max_tgs = MAX_THREAD_GROUP_SIZE;
    printf("INFO: Using maximum configured thread group size %d\n", MAX_THREAD_GROUP_SIZE);
  }

  // Get possible thread group sizes
  if(p->stencil_ctx.thread_group_size >= 1){ // thread group size passed by the user
    tgs_l = (int*) malloc(sizeof(int));
    tgs_l[0] = p->stencil_ctx.thread_group_size;
    n_tgs=1;

    tgs_factors_l = get_num_factors(p->stencil_ctx.thread_group_size, &n_tgs_factors);
  } else { // not passed by the user nor 1WD
    tgs_l = get_num_factors(max_tgs, &n_tgs);
    tgs_factors_l = tgs_l;
    n_tgs_factors = n_tgs;
  }

  // Get the other intra-tile params list if not set by the user
  if(p->stencil_ctx.th_z == -1){
    thz_l = tgs_factors_l;
    n_thz = n_tgs_factors;
  } else{
    thz_l = (int*) malloc(sizeof(int));
    thz_l[0] = p->stencil_ctx.th_z;
    n_thz = 1;
  }
  if(p->stencil_ctx.th_y == -1){
    thy_l = tgs_factors_l;
    n_thy = n_tgs_factors;
  } else{
    if(p->stencil_ctx.th_y>2) RAISE_ERROR("Threads along y-axis must be 1 or 2")
    thy_l = (int*) malloc(sizeof(int));
    thy_l[0] = p->stencil_ctx.th_y;
    n_thy = 1;
  }
  if(p->stencil_ctx.th_x == -1){
    thx_l = tgs_factors_l;
    n_thx = n_tgs_factors;
  } else{
    if(p->stencil_ctx.th_x>MAX_X_THREADS) RAISE_ERROR("Threads along x-axis must be less than 4")
    thx_l = (int*) malloc(sizeof(int));
    thx_l[0] = p->stencil_ctx.th_x;
    n_thx = 1;
  }

  // Get all feasible intra-tile parallelizm combinations
  n_tune_cases=0;
  for(tgsi=0; tgsi<n_tgs;tgsi++){
    for(z=0; z<n_thz; z++){
      for(y=0; y<n_thy; y++){
        for(x=0; x<n_thx; x++){
          if(thx_l[x]<=MAX_X_THREADS  && thy_l[y]<=2)
            if(thx_l[x]*thy_l[y]*thz_l[z] == tgs_l[tgsi]) n_tune_cases++;
        }
      }
    }
  }
  if(n_tune_cases == 0) RAISE_ERROR("Invalid thread group parallelism dimensions")
  tune_cases_l = (Tune_Params *) malloc(n_tune_cases*sizeof(Tune_Params));
  idx=0;
  for(tgsi=0; tgsi<n_tgs;tgsi++){
    for(z=0; z<n_thz; z++){
      for(y=0; y<n_thy; y++){
        for(x=0; x<n_thx; x++){
          if(thx_l[x]<=MAX_X_THREADS  && thy_l[y]<=2){
            if(thx_l[x]*thy_l[y]*thz_l[z] == tgs_l[tgsi]){
              tune_cases_l[idx].thread_group_size = tgs_l[tgsi];
              tune_cases_l[idx].th_z = thz_l[z];
              tune_cases_l[idx].th_y = thy_l[y];
              tune_cases_l[idx].th_x = thx_l[x];
              idx++;
            }
          }
        }
      }
    }
  }
  *num_tune_cases = n_tune_cases;
  *ret_tune_cases_l = tune_cases_l;
  *num_tgs = n_tgs;
  *ret_tgs_l = tgs_l;
}


void auto_tune_params(Parameters *p){
  int i, n_tune_cases, n_tgs, tune_case;
  Tune_Params *tune_cases_l=NULL;
  int *tgs_l=NULL;
  get_tgs_tune_params_lists(p, &tune_cases_l, &n_tune_cases, &tgs_l, &n_tgs);

  printf("Auto-tuning for the following cases (thread group size, th_z, th_y, th_x): ");
  for(i=0; i<n_tune_cases; i++)
    printf("(%d, %d, %d, %d) ", tune_cases_l[i].thread_group_size, tune_cases_l[i].th_z, tune_cases_l[i].th_y, tune_cases_l[i].th_x);
  printf("\n");


//  exit(0);
  // Get the maximum possible diamond width in the current configuration to initialize autotuning

  for(tune_case=0; tune_case<n_tune_cases; tune_case++){
    // set the current test cases
    p->stencil_ctx.thread_group_size = tune_cases_l[tune_case].thread_group_size;
    p->stencil_ctx.th_x = tune_cases_l[tune_case].th_x;
    p->stencil_ctx.th_y = tune_cases_l[tune_case].th_y;
    p->stencil_ctx.th_z = tune_cases_l[tune_case].th_z;


    // Initialize thread binings
    if(p->stencil_ctx.use_manual_cpu_bind==1){
      cpu_bind_init(p);
    }

    // auto-tune for diamond width and number of frontlines
    p->wf_larger_blk_size = 0;
    p->larger_t_dim = 0;
    if( ((p->t_dim == -1) || (p->stencil_ctx.num_wf==-1) )){
      if(p->mpi_rank == 0){
        auto_tune_diam_nwf(p);
      }
      if(p->t_dim==-1)
        RAISE_ERROR("EXITING: No feasible diamond width for this problem configurations")
    }

    // Free the CPU binding array after each test
    if(p->stencil_ctx.use_manual_cpu_bind==1){
      for(i=0; i<p->num_threads;i++){
        CPU_FREE(p->stencil_ctx.bind_masks[i]);
      }
      free(p->stencil_ctx.bind_masks);
    }

  }

  // broad cast the autotuning params
  if(p->mpi_size > 1){
    MPI_Bcast(&(p->t_dim), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(p->stencil_ctx.num_wf), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(p->stencil_ctx.thread_group_size), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(p->stencil_ctx.th_x), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(p->stencil_ctx.th_y), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(p->stencil_ctx.th_z), 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  p->wf_blk_size = get_mwf_size(*p, p->t_dim);

}


void intra_diamond_info_init(Parameters *p){
  int i, in_dimension=0;
  int nt2, remain, min_z;
  int nt = p->nt;
  int diam_concurrency, num_thread_groups;
  uint64_t diam_width;

  if(p->in_auto_tuning==0){
    p->in_auto_tuning = 1;
    auto_tune_params(p);
    p->in_auto_tuning = 0;
  }

  // Initialize thread binings
  if((p->stencil_ctx.use_manual_cpu_bind==1) && (p->in_auto_tuning==0)){
    cpu_bind_init(p);
  }

  diam_width = (p->t_dim+1) * 2 * p->stencil.r;
  diam_concurrency = (p->stencil_shape[1]/p->t.shape[1]) / diam_width;
  num_thread_groups = get_ntg(*p);
  if(num_thread_groups > diam_concurrency)
    if(p->mpi_rank ==0){
      printf("###ERROR: the number of thread groups exceed the available concurrency. Consider using %d thread groups or less\n", ((diam_concurrency>1)?diam_concurrency-1:1));
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }

  // check for thread assignment validity
  if(p->stencil_ctx.thread_group_size > p->num_threads){
    if(p->mpi_rank ==0){
      printf("###WARNING: Requested thread group size is larger the total available threads \n");
    }
  }

  // check thread group size validity
  if(p->stencil_ctx.thread_group_size != p->stencil_ctx.th_x * p->stencil_ctx.th_y*
                                         p->stencil_ctx.th_z * p->stencil_ctx.th_c){
     if(p->mpi_rank ==0){
      fprintf(stderr, "###ERROR: Thread group size must be consistent with parallelizm in all dimensions\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }
  }

  // check number of threads along x-axis vailidity
  if(p->stencil_ctx.th_x > p->lstencil_shape[0]){
     if(p->mpi_rank ==0){
      fprintf(stderr, "###ERROR: no sufficient concurrency along the x-axis\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }
  }
  // check if thread group sizes are equal
  if(p->num_threads%p->stencil_ctx.thread_group_size != 0){
    if(p->mpi_rank ==0){
      fprintf(stderr, "###ERROR: threads number must be multiples of thread group size\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }
  }

  // Check for size validity in the direction of the wavefront
  if( (p->wavefront == 1) && (p->stencil_ctx.thread_group_size == 1) ){ // single-thread group
    min_z = (p->t_dim*2)*p->stencil.r+1;
    if(p->stencil_shape[2] < min_z){
      if(p->mpi_rank ==0) fprintf(stderr,"ERROR: The single core wavefront requires a minimum size of %d at the Z direction in the current configurations\n", min_z);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }
  }

  if( (p->stencil_ctx.num_wf%p->stencil_ctx.th_z != 0) && (p->stencil_ctx.thread_group_size != 1) ){
    if(p->mpi_rank ==0) fprintf(stderr,"ERROR: num_wavefronts must be multiples of thread groups size\n");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
  }


  if(p->mwd_type == 3){ // relaxed synchronization
    if(p->stencil_ctx.num_wf/p->stencil_ctx.thread_group_size < p->stencil.r){
      if(p->mpi_rank ==0) fprintf(stderr,"ERROR: number of frontlines per thread must be greater or equal to the stencil radius in the relaxed synchronization implementation\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }
  }

  if( (p->wavefront == 1) && (p->stencil_ctx.thread_group_size != 1) ) // multi-thread group
    p->wavefront = -1;

    // Check for size validity in the direction of the wavefront
    if( (p->wavefront != 0) && (p->stencil_ctx.thread_group_size != 1) ){ // multi-thread group

      if(p->stencil.type==REGULAR)
        min_z = (p->t_dim*2)*p->stencil.r+1 + p->stencil_ctx.num_wf -1;
      else if (p->stencil.type==SOLAR)
        min_z = (p->t_dim*2+1)*p->stencil.r+1 + p->stencil_ctx.num_wf -1;

      if(p->stencil_shape[2] < min_z){
        if(p->mpi_rank ==0) fprintf(stderr,"ERROR: The multi-core wavefront requires a minimum size of %d at the Z direction in the current configurations\n", min_z);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(1);
      }
    }


    // Allocate the wavefront profiling timers
    num_thread_groups = get_ntg(*p);
    p->stencil_ctx.t_wait        = (double *) malloc(sizeof(double)*p->num_threads);
    p->stencil_ctx.t_wf_main     = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_wf_comm = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_wf_prologue = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_wf_epilogue = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.wf_num_resolved_diamonds = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_group_wait = (double *) malloc(sizeof(double)*num_thread_groups);

    // number of Stencil updates in the prologue and epilogue (trapezoid and inverted trapezoid)
//    p->idiamond_pro_epi_logue_updates = (uint64_t) (p->stencil_shape[0]/p->t.shape[0] * p->stencil_shape[2]/p->t.shape[2]) * (uint64_t) (2*diam_concurrency) * ((diam_width*diam_width)/4 - (diam_width/2));

    if(p->stencil.type == REGULAR){
      p->idiamond_pro_epi_logue_updates = (uint64_t) (p->stencil_shape[0]/p->t.shape[0] * p->stencil_shape[2]/p->t.shape[2]) * (uint64_t) (2*diam_concurrency) * ((p->t_dim+1)*(p->t_dim+1) + (p->t_dim+1))*p->stencil.r;
    }else if(p->stencil.type == SOLAR){
      p->idiamond_pro_epi_logue_updates = (uint64_t) (p->stencil_shape[0]/p->t.shape[0] * p->stencil_shape[2]/p->t.shape[2] * diam_concurrency)
                                         *((p->t_dim+1)*(p->t_dim+1)*2)*p->stencil.r;
    }

    if(p->source_point_enabled == 1){
      p->source_point_enabled = 0;
      if(p->mpi_rank ==0){
        printf("###INFO: Source point update disabled. Intra-diamond method does not support source point updates\n");
      }
    }
    if(p->t_dim < 1){
      if(p->mpi_rank == 0) fprintf(stderr,"ERROR: Diamond method does not support unrolling in time less than 1\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }
    if(p->t_dim%2 == 0){
      fprintf(stderr,"ERROR: diamond method does not supports even time unrolling\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }
    if( ( p->t.shape[0]>1 ) || (p->t.shape[2]>1) ){
      fprintf(stderr,"ERROR: Intra-diamond method supports domain decomposition across the Y direction only\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }

    // round the number of time steps to the nearest valid number
    if(p->stencil.type == REGULAR){
      remain = (nt-2)%((p->t_dim+1)*2);
    }else if(p->stencil.type == SOLAR){
      remain = (nt)%((p->t_dim+1)*2);
    }
    if(remain != 0){
      nt2 = nt + (p->t_dim+1)*2 - remain;
      if(nt2 != nt){
        if( (p->mpi_rank ==0) && (p->verbose ==1) ){
          printf("###INFO: Modified nt from %03d to %03d for the intra-diamond method to work properly\n", nt ,nt2);
        }
        p->nt= nt2;
      }
    }



  // count the topology dimensions containing the source point
  for(i=0; i<3; i++)
    if(((p->source_pt[i]-p->stencil.r) >= p->gb[i]) && ((p->source_pt[i]-p->stencil.r) <= p->ge[i])) in_dimension++;

  if(in_dimension == 3){
    p->has_source=1;
    for(i=0; i<3; i++) p->lsource_pt[i] = p->source_pt[i] - p->gb[i];
  } else{
    p->has_source=0;
    for(i=0; i<3; i++) p->lsource_pt[i] = -1;
  }

  // Update problem size information
  int t_dim = p->t_dim;
  if(p->t.shape[1] > 1)
    p->ldomain_shape[1] += (t_dim+1)*p->stencil.r;

  p->is_last = 0;
  if(p->t.rank_coords[1] == (p->t.shape[1]-1)) {
    p->is_last = 1;
    if(p->t.shape[1] > 1)
      p->ldomain_shape[1] += 2*p->stencil.r; //consider the interior boundary layers in the last MPI rank
//      if(p->lsource_pt[1] >= p->lstencil_shape[1]+p->stencil.r) p->lsource_pt[1] += 2*p->stencil.r; // consider the interior boundary region
  }

  if (p->lstencil_shape[1] < p->stencil.r*(t_dim+1)*2){
    fprintf(stderr,"ERROR: Intra-diamond method requires the sub-domain size to fit at least one diamond: %d elements in Y [stencil_radius*2*(time_unrolls+1)]. Given %d elements\n"
        ,p->stencil.r*(t_dim+1)*2, p->lstencil_shape[1]);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
  }
  if (floor(p->lstencil_shape[1] / (p->stencil.r*(t_dim+1)*2.0)) != p->lstencil_shape[1] / (p->stencil.r*(t_dim+1)*2.0)){
    if (p->mpi_rank ==0) fprintf(stderr,"ERROR: Intra-diamond method requires the sub-domain size to be multiples of the diamond width: %d elements [stencil_radius*2*(time_unrolls+1)]\n"
        ,p->stencil.r*(t_dim+1)*2);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
  }
}

