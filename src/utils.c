#include <string.h>
#include <getopt.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "utils.h"

void check_merr(int e) {
  switch(e){
    case 22: // EINVAL
      printf("alignment error\n");
    break;

    case 12: //ENOMEM
      printf("no sufficient memory\n");
    break;
  }

}

// set the default values of the kernel parameters
void param_default(Parameters *p) {
  p->stencil_shape[0] = 256;
  p->stencil_shape[1] = 64;
  p->stencil_shape[2] = 64;

  p->alignment = 16;
  p->target_ts = 0; //Naive TS
  p->target_kernel = 0; //Basic ISO stencil kernel
  p->stencil.r = stencil_info_list[p->target_kernel].r;
  p->n_tests = 3;
  p->nt = 100;
  p->verify = 0;
  p->source=NULL;
  p->verbose = 1;
  p->debug = 0;
  p->source_point_enabled = 0;
  p->array_padding = 1;
  p->mwd_type = 0;
  p->use_omp_stat_sched = 0;

  // diamond method
  p->t_dim = -1;
//  p->nb_diamond_chunk=1;
  p->halo_concat = 1;

  p->cache_size = 8192;

  // get the number of available OpenMP threads
  p->num_threads = 1;
#if defined(_OPENMP)
  p->num_threads = omp_get_max_threads();
#endif
  // default thread group size for methods using single thread group
  p->stencil_ctx.thread_group_size = p->num_threads;
  p->wavefront = 1; // default to using wavefront in the tile

  p->stencil_ctx.num_wf = -1;

  p->stencil_ctx.enable_likwid_m = 1;

  p->stencil_ctx.bs_x = 1e6;

  // Topology parameters
  p->t.is_periodic[0]=0;
  p->t.is_periodic[1]=0;
  p->t.is_periodic[2]=0;
  p->t.shape[0] = 1;
  p->t.shape[1] = 1;
  p->t.shape[2] = 1;

  p->h[0].is_contiguous = 1;
  p->h[1].is_contiguous = 1;
  p->h[2].is_contiguous = 1;

  p->h[0].size = 0;
  p->h[1].size = 0;
  p->h[2].size = 0;

  // Initialize the default stencil coefficients values
  FLOAT_PRECISION coef[] = {-0.28472, 0.16000, -0.02000, 0.00254,
      -0.00018, -0.18472, 0.19, -0.0500, 0.00554, -0.0009, 0.00354};
  int i;
  for(i=0; i<11; i++) p->g_coef[i] = (FLOAT_PRECISION) coef[i];

  // profiling information
  reset_timers(&(p->prof));
}
void reset_timers(Profile * p){
  p->compute = 0.;
  p->communicate = 0.;
  p->send_recv = 0.;
  p->wait = 0.;
  p->total = 0.;
  p->others = 0.;
  p->ts_main = 0.;
  p->ts_others = 0.;
}

void reset_wf_timers(Parameters * p){
  int i;
  int num_thread_groups = (int) ceil(1.0*p->num_threads/p->stencil_ctx.thread_group_size);

  // reset if the wavefront profiling is allocated
  if( (p->wavefront != 0) && (p->target_ts == 2) ) {

    for(i=0; i<p->num_threads; i++) p->stencil_ctx.t_wait[i] = 0.0;
    for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_wf_main[i] = 0.0;
    for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_wf_comm[i] = 0.0;
    for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_wf_prologue[i] = 0.0;
    for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_wf_epilogue[i] = 0.0;
    for(i=0; i<num_thread_groups; i++) p->stencil_ctx.wf_num_resolved_diamonds[i] = 0.0;
    for(i=0; i<num_thread_groups; i++) p->stencil_ctx.t_group_wait[i] = 0.0;
  }
}

void arrays_allocate(Parameters *p) {
  int male0, male1, male2;

  male1 = posix_memalign((void **)&(p->U1), p->alignment, sizeof(FLOAT_PRECISION)*p->ln_domain); check_merr(male1);
  male1 = posix_memalign((void **)&(p->U2), p->alignment, sizeof(FLOAT_PRECISION)*p->ln_domain); check_merr(male1);

  if(p->stencil.time_order == 2)
    male1 = posix_memalign((void **)&(p->U3), p->alignment, sizeof(FLOAT_PRECISION)*p->ln_domain); check_merr(male1);

  if(p->source_point_enabled==1)
    male0 = posix_memalign((void **)&(p->source), p->alignment, sizeof(FLOAT_PRECISION)*p->nt); check_merr(male0);

  // allocate the size of the coefficients matrix according to the stencil type
  unsigned long coef_size;
  switch(p->stencil.coeff){
  case CONSTANT_COEFFICIENT:
    coef_size = 10;
    break;

  case VARIABLE_COEFFICIENT:
    coef_size = p->ln_domain*(1 + p->stencil.r);
    break;

  case VARIABLE_COEFFICIENT_AXSYM:
    coef_size = p->ln_domain*(1 + 3*p->stencil.r);
    break;

  case VARIABLE_COEFFICIENT_NOSYM:
    coef_size = p->ln_domain*(1 + 6*p->stencil.r);
    break;

  default:
    printf("ERROR: unknown type of stencil\n");
    exit(1);
    break;
  }

  male2 = posix_memalign((void **)&(p->coef), p->alignment, sizeof(FLOAT_PRECISION)*coef_size); check_merr(male2);

  if (p->verbose == 1){
    if( (p->mpi_rank ==0)  || (p->mpi_rank == p->mpi_size-1))
    printf("[rank=%d] alloc. dom(err=%d):%fGiB coef(err=%d):%fGiB total:%fGiB\n", p->mpi_rank,
            male1, sizeof(FLOAT_PRECISION)*p->ln_domain*1.0/(1024*1024*1024),
            male2, sizeof(FLOAT_PRECISION)*coef_size*1.0/(1024*1024*1024),
            sizeof(FLOAT_PRECISION)*(coef_size+ 2*p->ln_domain)*1.0/(1024*1024*1024));
  }
}

void arrays_free(Parameters *p) {
  free(p->coef);
  free(p->U1);
  free(p->U2);

  if(p->stencil.time_order == 2)
      free(p->U3);

  if(p->source_point_enabled==1)
    free(p->source);
}

void set_centered_source(Parameters *p) {
  p->source_pt[0] = (p->stencil_shape[0]+2*p->stencil.r)/2 -1;
  p->source_pt[1] = (p->stencil_shape[1]+2*p->stencil.r)/2 -1;
  p->source_pt[2] = (p->stencil_shape[2]+2*p->stencil.r)/2 -1;
}

unsigned long get_mwf_size(Parameters *p, int t_dim){
  unsigned long diam_width, diam_height, wf_updates, wf_elements, lnx, t_order, total_points;

  t_order = p->stencil.time_order;
  diam_width = (t_dim+1)*2*p->stencil.r;
  int nwf = p->stencil_ctx.num_wf;
  diam_height = t_dim*2*p->stencil.r + nwf;

  int bsl = p->stencil_ctx.bs_x + t_dim*p->stencil.r;
  lnx = (bsl < p->ldomain_shape[0] ? bsl : p->ldomain_shape[0]);

  wf_updates = (t_dim+1)*(t_dim+1)*2 * p->stencil.r; // Y-T projection
  wf_elements = (wf_updates - diam_width) * p->stencil.r + diam_width + diam_width*(nwf-1);

  switch(p->stencil.coeff){
  case CONSTANT_COEFFICIENT:
    total_points = ( (t_order+1)             *wf_elements + (diam_width + diam_height )*2*p->stencil.r) * lnx * sizeof(FLOAT_PRECISION);
    break;

  case VARIABLE_COEFFICIENT:
    total_points = ( (t_order+1 + (1+p->stencil.r) )*wf_elements + (diam_width + diam_height )*2*p->stencil.r) * lnx * sizeof(FLOAT_PRECISION);
    break;

  case VARIABLE_COEFFICIENT_AXSYM:
    total_points = ( (t_order+1 + (1+3*p->stencil.r) )*wf_elements + (diam_width + diam_height )*2*p->stencil.r) * lnx * sizeof(FLOAT_PRECISION);
    break;

  case VARIABLE_COEFFICIENT_NOSYM:
    total_points = ( (t_order+1 + (1+6*p->stencil.r) )*wf_elements + (diam_width + diam_height )*2*p->stencil.r) * lnx * sizeof(FLOAT_PRECISION);
    break;

  default:
    printf("ERROR: unknown type of stencil\n");
    exit(1);
    break;
  }
  //printf("npx:%d  nx:%d  lnx:%lu  updates:%lu  elements:%lu  total:%lu\n", p->t.shape[0],  p->ldomain_shape[0], lnx, wf_updates, wf_elements, total_points);
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

  unsigned long lups = (tp->ln_stencils*tp->nt - tp->idiamond_pro_epi_logue_updates);
  obt_perf =  lups/tp->prof.ts_main;

  printf("%06.2f]  time:%es  lups:%lu  cache block size:%lukiB  reps:%d\n",
      obt_perf/(1e6), tp->prof.ts_main, lups, get_mwf_size(tp, tp->t_dim)*(tp->num_threads/tp->stencil_ctx.thread_group_size)/1024, orig_reps);

  return obt_perf;
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
  unsigned long wf_size, ntg;
  int tgs = tp.stencil_ctx.thread_group_size;
  int cache_size_cond, int_diam_cond, wf_len_cond, cuncurrency_cond, diam_concurrency;
  int diam_height;
  ntg = tp.num_threads / tgs;

  tp.stencil_ctx.num_wf = tgs; // set the number of wavefronts to the minimum possible value
  if((tp.mwd_type == 2) | (tp.mwd_type == 3)) tp.stencil_ctx.num_wf = tgs*tp.stencil.r;


  // find max possible diamond width and allocate memory accordingly, if not pre set
  if(op->t_dim == -1){
    max_t_dim = -1;
    for(i=tp.stencil_shape[1]/tp.stencil.r; i>=4; i-=4){
      lt_dim = i/2 - 1;
      diam_width = i*tp.stencil.r;
      diam_height = lt_dim*2*tp.stencil.r+1 + tp.stencil_ctx.num_wf-1;

      wf_size = get_mwf_size(&tp, lt_dim);
      cache_size_cond = wf_size*ntg < (unsigned long) (MAX_CACHE_SIZE*1024);

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
  int num_thread_groups = (int) ceil(1.0*tp.num_threads/tp.stencil_ctx.thread_group_size);
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
      wf_size = get_mwf_size(&tp, tp.t_dim);
      diam_concurrency = tp.stencil_shape[1]/diam_width;

      cache_size_cond = wf_size*ntg > (unsigned long) (tp.cache_size*1024);
      cuncurrency_cond = diam_concurrency >= ntg;
      int_diam_cond = tp.stencil_shape[1]%diam_width == 0;
  //    printf("i:%d, diam_width %d,  cuncurrency_cond %d, cache_size_cond %d, int_diam_cond %d, wf_len_cond %d, cache_blk_size: %lu kB\n",
  //        i, diam_width, cuncurrency_cond, cache_size_cond, int_diam_cond, wf_len_cond, wf_size*ntg/1024);

      if( (int_diam_cond == 1) && (cuncurrency_cond == 1) && (cache_size_cond == 1) ){ // check diamond size validity
        printf("[AUTO TUNE] Diamond width:%02d  [wavefronts #: pefromance (MLUPS/s)]\n", (tp.t_dim+1)*2*tp.stencil.r);
        // loop over increasing number of wavefronts per update
        prev_nwf_perf = -1;
        tp.stencil_ctx.num_wf = tgs; // start with smallest possible number of updates
        if( (tp.mwd_type == 2) | (tp.mwd_type == 3) ) tp.stencil_ctx.num_wf = tgs*tp.stencil.r;
        tp.idiamond_pro_epi_logue_updates = (unsigned long) (tp.stencil_shape[0] * tp.stencil_shape[2]) * (unsigned long) (2*diam_concurrency) * ((tp.t_dim+1)*(tp.t_dim+1) + (tp.t_dim+1))*tp.stencil.r;

        while(1){
          wf_size = get_mwf_size(&tp, tp.t_dim);
          cache_size_cond = wf_size*ntg > (unsigned long) (tp.cache_size*1024);
          diam_height = tp.t_dim*2*tp.stencil.r+1 +tp.stencil_ctx.num_wf-1;
          wf_len_cond = diam_height <= tp.stencil_shape[2];

          if( (wf_len_cond==1) && (cache_size_cond == 1) ){
            printf("[AUTO TUNE]     [%03d: ",tp.stencil_ctx.num_wf);  
            exp_perf = run_tuning_test(&tp);

  //          printf("tgs:%d  nwf:%d  perf:%6.2f  prev_nwf_perf:%6.2f\n", tgs, tp.stencil_ctx.num_wf, exp_perf/1024/1024, prev_nwf_perf/1024/1024);

            // termination criteria for the nwf
            if (exp_perf < prev_nwf_perf){
              tp.stencil_ctx.num_wf -= tgs;
              break;
            }
            else{
              prev_nwf_perf = exp_perf;
              tp.stencil_ctx.num_wf += tgs;
            }

          }
          else{ // invalid wavefront length
            tp.stencil_ctx.num_wf -= tgs;
            break;
          }

        }
        if(tp.stencil_ctx.num_wf < tgs){
          tp.stencil_ctx.num_wf = tgs;
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
    tp.idiamond_pro_epi_logue_updates = (unsigned long) (tp.stencil_shape[0] * tp.stencil_shape[2]) * (unsigned long) (2*diam_concurrency) * ((tp.t_dim+1)*(tp.t_dim+1) + (tp.t_dim+1))*tp.stencil.r;

    while(1){
      diam_height = tp.t_dim*2*tp.stencil.r+1 +tp.stencil_ctx.num_wf-1;
      wf_len_cond = diam_height <= tp.stencil_shape[2];

      if(wf_len_cond==1){
        printf("[AUTO TUNE]     [%03d: ",tp.stencil_ctx.num_wf);  
        exp_perf = run_tuning_test(&tp);

//          printf("tgs:%d  nwf:%d  perf:%6.2f  prev_nwf_perf:%6.2f\n", tgs, tp.stencil_ctx.num_wf, exp_perf/1024/1024, prev_nwf_perf/1024/1024);

        // termination criteria for the nwf
        if (exp_perf < prev_nwf_perf){
          tp.stencil_ctx.num_wf -= tgs;
          break;
        }
        else{
          prev_nwf_perf = exp_perf;
          tp.stencil_ctx.num_wf += tgs;
        }

      }
      else{ // invalid wavefront length
        tp.stencil_ctx.num_wf -= tgs;
        break;
      }

    }
    if(tp.stencil_ctx.num_wf < tgs){
      tp.stencil_ctx.num_wf = tgs;
    }
    if( (tp.mwd_type == 2) | (tp.mwd_type == 3) ){ 
      if(tp.stencil_ctx.num_wf < tgs*tp.stencil.r){
        tp.stencil_ctx.num_wf = tgs*tp.stencil.r;
      }
    }
    op->stencil_ctx.num_wf = tp.stencil_ctx.num_wf;
  }

  //simple tuning of blocking in X
  diam_concurrency = tp.stencil_shape[1]/((tp.t_dim+1)*2*tp.stencil.r);
  tp.idiamond_pro_epi_logue_updates = (unsigned long) (tp.stencil_shape[0] * tp.stencil_shape[2]) * (unsigned long) (2*diam_concurrency) * ((tp.t_dim+1)*(tp.t_dim+1) + (tp.t_dim+1))*tp.stencil.r;
  int prev_bs_x, base_bs_x;
  double div=1.0;
  if (tp.stencil_ctx.bs_x > tp.stencil_shape[0]/tp.t.shape[0]){//no blocking set by user
    tp.stencil_ctx.bs_x = tp.stencil_shape[0]/tp.t.shape[0];
    base_bs_x = tp.stencil_ctx.bs_x;
    printf("[AUTO TUNE]     Diamond width:%02d, wavefronts #:%d  [Blk. in X: pefromance (MLUPS/s)]\n", (tp.t_dim+1)*2*tp.stencil.r, tp.stencil_ctx.num_wf);
    printf("[AUTO TUNE]          [%03d: ", tp.stencil_ctx.bs_x);  
    prev_nwf_perf = run_tuning_test(&tp);
    while(1){
      div += 1;
      prev_bs_x = tp.stencil_ctx.bs_x; 
      tp.stencil_ctx.bs_x = ceil( ((double)(base_bs_x)) /div);
      printf("[AUTO TUNE]          [%03d: ", tp.stencil_ctx.bs_x);  
      exp_perf = run_tuning_test(&tp);

      if(exp_perf < prev_nwf_perf){
        op->stencil_ctx.bs_x = prev_bs_x;
        break;
      }
      prev_nwf_perf = exp_perf;
    }
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

void set_kernels(Parameters *p){

#if USE_SPLIT_STRIDE // separate central line update
  p->stencil.spt_blk_func = iso_ref_split; // Note keeping the stat. sched. same

  if(p->stencil_ctx.thread_group_size == 1){ //1WD use implementation
    p->stencil.mwd_func = swd_iso_ref_split;
  } else {
    p->stencil.mwd_func = mwd_iso_ref_split;
  }

  // set central line updates kernels
  p->stencil_ctx.clu_func = clu_func_list[p->target_kernel];

#else

  p->stencil.spt_blk_func = spt_blk_func_list[p->target_kernel];

  // use static openmp schedule if set for spatial blocking time steppers
  if( (p->target_ts==0 || p->target_ts==1) && (p->use_omp_stat_sched==1) ){      
    p->stencil.spt_blk_func = stat_sched_func_list[p->target_kernel];
  }

  p->stencil.mwd_func = mwd_list[p->mwd_type][p->target_kernel];



  if(p->stencil_ctx.thread_group_size == 1){ //1WD use implementation
    p->stencil.mwd_func = swd_func_list[p->target_kernel];
  }
#endif

  p->stencil.name = stencil_info_list[p->target_kernel].name;
  p->stencil.coeff = stencil_info_list[p->target_kernel].coeff;
  p->stencil.r = stencil_info_list[p->target_kernel].r;
  p->stencil.time_order = stencil_info_list[p->target_kernel].time_order;
  p->stencil.shape = stencil_info_list[p->target_kernel].shape;
  p->stencil.stat_sched_func = stat_sched_func_list[p->target_kernel];

}

void standard_info_init(Parameters *p){
  int i, in_dimension=0;
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
}
void intra_diamond_info_init(Parameters *p){
  int i, in_dimension=0;
  int nt2, remain, min_z;
  int nt = p->nt;
  int diam_concurrency, num_thread_groups;
  unsigned long diam_width;


    p->wf_larger_blk_size = 0;
    p->larger_t_dim = 0;

    if( (p->stencil_ctx.thread_group_size !=-1) && ((p->t_dim == -1) || (p->stencil_ctx.num_wf==-1) )){ // thread group size is defined but not the diamond size
      if(p->mpi_rank == 0){
        p->stencil_ctx.enable_likwid_m = 0;
        auto_tune_diam_nwf(p);
        p->stencil_ctx.enable_likwid_m = 1;
      }
      if(p->mpi_size > 1){ // broad cast the autotuning params
        MPI_Bcast(&(p->t_dim), 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&(p->stencil_ctx.num_wf), 1, MPI_INT, 0, MPI_COMM_WORLD);
      }
      if(p->t_dim == -1){
        if(p->mpi_rank == 0) 
          printf("EXITING: No feasible diamond width for this problem configurations\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(0);
      }
    }
    p->wf_blk_size = get_mwf_size(p, p->t_dim);


    diam_width = (p->t_dim+1) * 2 * p->stencil.r;
    diam_concurrency = (p->stencil_shape[1]/p->t.shape[1]) / diam_width;

    if (p->stencil_ctx.thread_group_size ==-1) { // setup default thread group information
      if(p->num_threads == 1)
        p->stencil_ctx.thread_group_size = 1;
      else
        p->stencil_ctx.thread_group_size = ( (diam_concurrency<2) ? p->num_threads : (p->num_threads/2) );
    }
    else { // thread group size is set by the user
      num_thread_groups = (int) ceil(1.0*p->num_threads/p->stencil_ctx.thread_group_size);
      //printf("threads:%d  num_thread_groups:%d thread_group_size:%d width:%d stencils:%d concurrency:%d\n",p->num_threads, num_thread_groups, p->stencil_ctx.thread_group_size, diam_width, p->stencil_shape[1]/p->t.shape[1], diam_concurrency);
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
    }

    // check for block size in X validity
    if(p->stencil_ctx.bs_x%p->stencil.r != 0) {
      if(p->mpi_rank ==0){
        fprintf(stderr, "###ERROR: Block size in in X must be multiples of the stencil radius\n");
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(1);
      }
    }

    // check if thread group sizes are equal
    if(p->num_threads%p->stencil_ctx.thread_group_size != 0){
      if(p->mpi_rank ==0){
        fprintf(stderr, "###ERROR: thread groups must be equal in size when using pre-computed wavefront assignment\n");
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

    if( (p->stencil_ctx.num_wf%p->stencil_ctx.thread_group_size != 0) && (p->stencil_ctx.thread_group_size != 1) ){
      if(p->mpi_rank ==0) fprintf(stderr,"ERROR: num_wavefronts must be multiples of thread groups size\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }


  if (p->mwd_type == 1){ // fixed execution to data
    if( (p->stencil_ctx.num_wf%p->stencil_ctx.thread_group_size != 0) && (p->stencil_ctx.thread_group_size != 1) ){
      if(p->mpi_rank ==0) fprintf(stderr,"ERROR: number of wavefronts must be multiples of thread group size\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(1);
    }
  }

  if( (p->mwd_type == 2) | (p->mwd_type == 3) ){ // relaxed synchronization
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

      min_z = (p->t_dim*2)*p->stencil.r+1 + p->stencil_ctx.num_wf -1;
      if(p->stencil_shape[2] < min_z){
        if(p->mpi_rank ==0) fprintf(stderr,"ERROR: The multi-core wavefront requires a minimum size of %d at the Z direction in the current configurations\n", min_z);

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
        exit(1);
      }
    }


    // Allocate the wavefront profiling timers
    num_thread_groups = (int) ceil(1.0*p->num_threads/p->stencil_ctx.thread_group_size);
    p->stencil_ctx.t_wait        = (double *) malloc(sizeof(double)*p->num_threads);
    p->stencil_ctx.t_wf_main     = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_wf_comm = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_wf_prologue = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_wf_epilogue = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.wf_num_resolved_diamonds = (double *) malloc(sizeof(double)*num_thread_groups);
    p->stencil_ctx.t_group_wait = (double *) malloc(sizeof(double)*num_thread_groups);
    
    // number of Stencil updates in the prologue and epilogue (trapezoid and inverted trapezoid)
//    p->idiamond_pro_epi_logue_updates = (unsigned long) (p->stencil_shape[0]/p->t.shape[0] * p->stencil_shape[2]/p->t.shape[2]) * (unsigned long) (2*diam_concurrency) * ((diam_width*diam_width)/4 - (diam_width/2));
    p->idiamond_pro_epi_logue_updates = (unsigned long) (p->stencil_shape[0]/p->t.shape[0] * p->stencil_shape[2]/p->t.shape[2]) * (unsigned long) (2*diam_concurrency) * ((p->t_dim+1)*(p->t_dim+1) + (p->t_dim+1))*p->stencil.r;

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
    remain = (nt-2)%((p->t_dim+1)*2);
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
void init(Parameters *p) {
  int q, r, i;

  set_kernels(p);
  p->n_stencils = p->stencil_shape[0] * p->stencil_shape[1] * p->stencil_shape[2];


  if(p->stencil.r > 10){
    if(p->mpi_rank == 0) fprintf(stderr,"ERROR: Stencil operators with radius greater than 10 are not supported\n");
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    exit(1);
  }

  if( (p->mpi_size == 1) && (p->halo_concat ==1) ){
    p->halo_concat = 0;
    if(p->verbose==1) printf("###INFO: Halo concatenation disabled. It does not make sense in single process run\n");
  }

  p->source_point_enabled = 0;

  // compute the global boundaries of each subdomain
  for(i=0; i<3; i++) {
    if(p->t.shape[i] > 1) {
      q = (int)(p->stencil_shape[i]/p->t.shape[i]);
      r = p->stencil_shape[i] % p->t.shape[i];
      if(p->t.rank_coords[i] < r) {
        p->lstencil_shape[i] = q+1;
        p->gb[i] = p->t.rank_coords[i] * (q+1);
      }else {
        p->lstencil_shape[i] = q;
        p->gb[i] = r * (q+1) +
            (p->t.rank_coords[i] - r) * q;
      }
    } else { // i.e. one rank in this dimension
      p->lstencil_shape[i] = p->stencil_shape[i];
      p->gb[i] = 0;
    }
    p->ge[i] = p->gb[i] + p->lstencil_shape[i] - 1;
  }

  int padding_comp, padding_size = 0;
  if (p->array_padding == 1) {
    padding_comp = (p->lstencil_shape[0]+2*p->stencil.r)%p->alignment;
    if (padding_comp != 0) padding_size = p->alignment - padding_comp;
  }
  p->ldomain_shape[0] = p->lstencil_shape[0]+2*p->stencil.r + padding_size;
  p->ldomain_shape[1] = p->lstencil_shape[1]+2*p->stencil.r;
  p->ldomain_shape[2] = p->lstencil_shape[2]+2*p->stencil.r;

  // calculate the block size in Y to satisfy the layer condition at the spatially blocked code
  if (p->cache_size >0){

    int num_thread_groups = (int) ceil(1.0*p->num_threads / p->stencil_ctx.thread_group_size);

    if(p->use_omp_stat_sched==0){
      p->stencil_ctx.bs_y = (p->cache_size*1024)/((num_thread_groups* (p->stencil_ctx.thread_group_size+(2*p->stencil.r)))*p->ldomain_shape[0]*sizeof(FLOAT_PRECISION));
    } else {// tailored for the Xeon Phi
      p->stencil_ctx.bs_y = (p->cache_size*1024)/(p->stencil_ctx.thread_group_size*(1+2*p->stencil.r)*p->ldomain_shape[0]*sizeof(FLOAT_PRECISION));
    }
    // set minimum block size if cache is not sufficient
    if(p->stencil_ctx.bs_y == 0) p->stencil_ctx.bs_y=1;

  } else {
    p->stencil_ctx.bs_y = 1000000; // make the block larger than the domain
  }

  // set the local source point and other information
  switch(p->target_ts){
  case 0: //set source information for standard methods
  case 1:
    standard_info_init(p);
    break;
  case 2: // intra-diamond methods
    intra_diamond_info_init(p);
    break;
  }

  p->ln_domain = p->ldomain_shape[0] * p->ldomain_shape[1]* p->ldomain_shape[2];
  p->ln_stencils = p->lstencil_shape[0] * p->lstencil_shape[1] * p->lstencil_shape[2];


  if(p->debug ==1){
    MPI_Barrier(MPI_COMM_WORLD);
    if(p->mpi_rank == 0) {
      printf("\n******************************************************\n"); fflush(stdout);
      printf("DEBUG domain decomposition information BEGIN\n"); fflush(stdout);
      printf("******************************************************\n"); fflush(stdout); sleep(1);
    }
    int j,k;
    for(j=0; j<p->mpi_size; j++){
      if(j == p->mpi_rank){

        printf("[%02d]:top(%02d,%02d,%02d)\n", p->mpi_rank, p->t.rank_coords[0], p->t.rank_coords[1], p->t.rank_coords[2]); fflush(stdout);
        printf("ln_domain:%06lu  lnstencil:%06lu\n", p->ln_domain, p->ln_stencils); fflush(stdout);
        printf("  Local stencil Shape:(%03d,%03d,%03d)\n", p->lstencil_shape[0], p->lstencil_shape[1], p->lstencil_shape[2]); fflush(stdout);
        printf("  Local domain Shape: (%03d,%03d,%03d)\n", p->ldomain_shape[0], p->ldomain_shape[1], p->ldomain_shape[2]); fflush(stdout);
        printf("  Local begin:        (%03d,%03d,%03d)\n", p->gb[0], p->gb[1], p->gb[2]); fflush(stdout);
        printf("  Local end:          (%03d,%03d,%03d)\n", p->ge[0], p->ge[1], p->ge[2]); fflush(stdout);
        printf("  Local source point: (%03d,%03d,%03d)\n", p->lsource_pt[0]-p->stencil.r, p->lsource_pt[1]-p->stencil.r, p->lsource_pt[2]-p->stencil.r); fflush(stdout);

        printf("\n"); fflush(stdout);
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if(p->mpi_rank == 0) {
      printf("******************************************************\n"); fflush(stdout);
      printf("DEBUG domain decomposition information END\n"); fflush(stdout);
      printf("******************************************************\n\n"); fflush(stdout);
    }
  }
}


void init_coeff(Parameters * p) {
  int i, k, ax;

  switch(p->stencil.coeff){
  case CONSTANT_COEFFICIENT:
    for(i=0;i<p->stencil.r+1;i++)
        p->coef[i] = p->g_coef[i];
    break;

  case VARIABLE_COEFFICIENT:
    for(k=0; k <= p->stencil.r; k++){
      for(i=0; i<p->ln_domain; i++){
        p->coef[i + k*p->ln_domain] = p->g_coef[k];
      }
    }
    break;

  case VARIABLE_COEFFICIENT_AXSYM:
    // central point coeff
    for(i=0; i<p->ln_domain; i++){
      p->coef[i] = p->g_coef[0];
    }
    for(k=0; k < p->stencil.r; k++){
      for(ax=0; ax<3; ax++){
        for(i=0; i<p->ln_domain; i++){
          p->coef[i + p->ln_domain + 3*k*p->ln_domain + ax*p->ln_domain] = p->g_coef[k+1];
        }
      }
    }
    break;

  case VARIABLE_COEFFICIENT_NOSYM:
    // central point coeff
    for(i=0; i<p->ln_domain; i++){
      p->coef[i] = p->g_coef[0];
    }
    for(k=0; k < p->stencil.r; k++){
      for(ax=0; ax<3; ax++){
        for(i=0; i<p->ln_domain; i++){
          p->coef[i + p->ln_domain + 6*k*p->ln_domain +  2*ax   *p->ln_domain] = p->g_coef[k+1];
          p->coef[i + p->ln_domain + 6*k*p->ln_domain + (2*ax+1)*p->ln_domain] = p->g_coef[k+1];
        }
      }
    }
    break;

  default:
    printf("ERROR: unknown type of stencil\n");
    exit(1);
    break;
  }
}

void copyv(int *a, int * b, int n) {
  int i;
  for(i=0;i<n;i++) b[i] = a[i];
}
void copy_halo_struct(Halo a, Halo *b) {
  copyv(a.shape, b->shape, 3);
  copyv(a.recv_b, b->recv_b, 3);
  copyv(a.recv_e, b->recv_e, 3);
  copyv(a.send_b, b->send_b, 3);
  copyv(a.send_e, b->send_e, 3);
  b->recv_hb = a.recv_hb;
  b->recv_he = a.recv_he;
  b->send_hb = a.send_hb;
  b->send_he = a.send_he;
}
void copy_params_struct(Parameters a, Parameters * b) {

  b->alignment =a.alignment;
  b->verbose =a.verbose;
  b->debug =a.debug;
  copyv(a.stencil_shape, b->stencil_shape, 3);
  b->n_stencils = a.n_stencils;
  b->ln_domain = a.ln_domain;
  b->target_ts = a.target_ts;
  b->target_kernel = a.target_kernel;
  b->source_point_enabled = a.source_point_enabled;
  b->mpi_rank = a.mpi_rank;
  b->mpi_size = a.mpi_size;
  b->n_tests = a.n_tests;
  b->nt = a.nt;
  b->verify = a.verify;
  b->num_threads = a.num_threads;
  b->array_padding = a.array_padding;
  b->mwd_type = a.mwd_type;
  b->use_omp_stat_sched = a.use_omp_stat_sched;
  b->stencil_ctx.bs_y = a.stencil_ctx.bs_y;
  b->stencil_ctx.bs_x = a.stencil_ctx.bs_x;
  b->stencil_ctx.thread_group_size = a.stencil_ctx.thread_group_size;
  b->stencil_ctx.clu_func = a.stencil_ctx.clu_func;
  b->stencil_ctx.num_wf = a.stencil_ctx.num_wf;

  
  copyv(a.source_pt, b->source_pt, 3);
  copyv(a.lstencil_shape, b->lstencil_shape, 3);
  copyv(a.ldomain_shape, b->ldomain_shape, 3);
  copyv(a.gb, b->gb, 3);
  copyv(a.ge, b->ge, 3);
  copyv(a.lsource_pt, b->lsource_pt, 3);
  b->has_source = a.has_source;
  b->halo_concat = a.halo_concat;
  b->wavefront = a.wavefront;
  b->idiamond_pro_epi_logue_updates = a.idiamond_pro_epi_logue_updates;
  b->t_dim = a.t_dim;
  b->is_last = a.is_last;

  b->cache_size =   a.cache_size;

  //   * restrict U1, * restrict U2, * restrict U3, * restrict coef, * restrict source;
  //
  //  copy_halo_struct(a.h[0],&(b->h[0]));
  //  copy_halo_struct(a.h[1],&(b->h[1]));
  //  copy_halo_struct(a.h[2],&(b->h[2]));
  b->h[0].is_contiguous = a.h[0].is_contiguous;
  b->h[1].is_contiguous = a.h[1].is_contiguous;
  b->h[2].is_contiguous = a.h[2].is_contiguous;

  b->h[0].size = a.h[0].size;
  b->h[1].size = a.h[1].size;
  b->h[2].size = a.h[2].size;

  b->t.right = a.t.right;
  b->t.left = a.t.left;
  b->t.up = a.t.up;
  b->t.down = a.t.down;
  b->t.front = a.t.front;
  b->t.back = a.t.back;
  copyv(a.t.shape, b->t.shape, 3);
  copyv(a.t.is_periodic, b->t.is_periodic, 3);
  copyv(a.t.rank_coords, b->t.rank_coords, 3);
  b->t.cart_comm = a.t.cart_comm;

  b->stencil.name = a.stencil.name;
  b->stencil.r = a.stencil.r;
  b->stencil.time_order = a.stencil.time_order;
  b->stencil.shape = a.stencil.shape;
  b->stencil.coeff = a.stencil.coeff;
  b->stencil.spt_blk_func  = a.stencil.spt_blk_func;
  b->stencil.stat_sched_func = a.stencil.stat_sched_func;
  b->stencil.mwd_func = a.stencil.mwd_func;

}

#define pU1(i,j,k)          (p->U1[((k)*(p->ldomain_shape[1])+(j))*(p->ldomain_shape[0])+(i)])
#define pU2(i,j,k)          (p->U2[((k)*(p->ldomain_shape[1])+(j))*(p->ldomain_shape[0])+(i)])
#define pU3(i,j,k)          (p->U3[((k)*(p->ldomain_shape[1])+(j))*(p->ldomain_shape[0])+(i)])
void domain_data_fill(Parameters * p){
  int i,j,k, gi, gj, gk;
  FLOAT_PRECISION r;
  int xb, xe, yb, ye, zb, ze;
  for(i=0; i<p->ln_domain;i++){
    p->U1[i] = 0.0;
    p->U2[i] = 0.0;

    if(p->stencil.time_order == 2)
      p->U3[i] = 0.0;
  }
  xb = 0;
  yb = 0;
  zb = 0;
  xe = p->lstencil_shape[0]+2*p->stencil.r;
  ye = p->lstencil_shape[1]+2*p->stencil.r;
  ze = p->lstencil_shape[2]+2*p->stencil.r;
  if(p->t.rank_coords[0] == 0) xb += p->stencil.r;
  if(p->t.rank_coords[1] == 0) yb += p->stencil.r;
  if(p->t.rank_coords[2] == 0) zb += p->stencil.r;
  if(p->t.rank_coords[0] == p->t.shape[0]-1) xe -= p->stencil.r;
  if(p->t.rank_coords[1] == p->t.shape[1]-1) ye -= p->stencil.r;
  if(p->t.rank_coords[2] == p->t.shape[2]-1) ze -= p->stencil.r;
  // fill the local stencil subdomain according to the global location and pad the boundary with zeroes
  for(k=zb; k<ze; k++){
    for(j=yb; j<ye; j++){
      for(i=xb; i<xe; i++){
        gi = i + p->gb[0];
        gj = j + p->gb[1];
        gk = k + p->gb[2];

        r = 1.0/3 * (1.0*gi/p->stencil_shape[0] + 1.0*gj/p->stencil_shape[1]  + 1.0*gk/p->stencil_shape[2]);

        pU1(i, j, k) = r*1.845703;
        pU2(i, j, k) = r*1.845703;
        if(p->stencil.time_order == 2)
          pU3(i, j, k) = r*1.845703;
      }
    }
  }
  // fill the extra allocation of diamond tiles
  int gb_1;
  if( (p->target_ts == 2) && (p->stencil.time_order == 2) && p->mpi_size > 1){
    // last process extends to the beginning of the domain
    if(p->t.rank_coords[1] == p->t.shape[1]-1) {
      yb = ye + 2*p->stencil.r;
      ye = yb + (p->t_dim+1+1)*p->stencil.r;
      gb_1 = -yb+p->stencil.r;
    } else {
      yb = ye;
      ye += (p->t_dim+1)*p->stencil.r;
      gb_1 = p->gb[1];
    }
    for(k=zb; k<ze; k++){
      for(j=yb; j<ye; j++){
        for(i=xb; i<xe; i++){
          gi = i + p->gb[0];
          gj = j + gb_1;
          gk = k + p->gb[2];

          r = 1.0/3 * (1.0*gi/p->stencil_shape[0] + 1.0*gj/p->stencil_shape[1]  + 1.0*gk/p->stencil_shape[2]);

          pU1(i, j, k) = r*1.845703;
          pU2(i, j, k) = r*1.845703;
          if(p->stencil.time_order == 2)
            pU3(i, j, k) = r*1.845703;
        }
      }
    }
  }

  // set source points at the first and last YZ plains
  if(p->t.rank_coords[0] == 0){
    for(k=0; k<p->ldomain_shape[2]; k++){
      for(j=0; j<p->ldomain_shape[1]; j++){
        pU1(0, j, k) += BOUNDARY_SRC_VAL;
        pU2(0, j, k) += BOUNDARY_SRC_VAL;
      }
    }
  }
  if(p->t.rank_coords[0] == p->t.shape[0]-1){
    for(k=0; k<p->ldomain_shape[2]; k++){
      for(j=0; j<p->ldomain_shape[1]; j++){
        pU1(p->lstencil_shape[0]+2*p->stencil.r-1, j, k) += BOUNDARY_SRC_VAL;
        pU2(p->lstencil_shape[0]+2*p->stencil.r-1, j, k) += BOUNDARY_SRC_VAL;
      }
    }
  }
}

void performance_results(Parameters *p, double t, double t_max, double t_min, double t_ts_main_max, double t_ts_main_min){
  double max_comm, max_comp, max_wait, max_others, max_total;
  double min_comm, min_comp, min_wait, min_others, min_total;
  double mean_comm, mean_comp, mean_wait, mean_others, mean_total;


  int i;
  unsigned long total_stencils;
  int num_thread_groups = (int) ceil(1.0*p->num_threads/p->stencil_ctx.thread_group_size);

  // look for NANs and zero results
  int k, zeroes_p, n_zeroes=0, nans=0;
  for (k=0; k<p->ln_domain; k++){
    // check for nan and -inf/inf
    nans += (p->U1[k] * 0) != 0;
    // check for zeroes
#if DP
      if(fabsf(p->U1[k]) < 1e-6) n_zeroes++;
#else
      if(fabs(p->U1[k]) < 1e-6) n_zeroes++;
#endif
  }
  zeroes_p = 100*n_zeroes/p->ln_domain;
  if(zeroes_p > 90){
    printf("\n******************************************************\n");
    printf("##WARNING[rank:%d]: %d%% of the sub domain contains zeroes. This might result in inaccurate performance results\n", p->mpi_rank, zeroes_p);
    printf("******************************************************\n\n");
  }
  if(nans > 0){
    printf("\n******************************************************\n");
    printf("##WARNING[rank:%d]: %d nan and/or -inf/inf values in the final sub domain solution. This might result in inaccurate performance results\n", p->mpi_rank, nans);
    printf("******************************************************\n\n");
  }


  // Collect timing statistics
  ierr = MPI_Reduce(&(p->prof.compute), &max_comp, 1, MPI_DOUBLE, MPI_MAX, 0, p->t.cart_comm); CHKERR(ierr);
  ierr = MPI_Reduce(&(p->prof.communicate), &max_comm, 1, MPI_DOUBLE, MPI_MAX, 0, p->t.cart_comm); CHKERR(ierr);
  ierr = MPI_Reduce(&(p->prof.wait), &max_wait, 1, MPI_DOUBLE, MPI_MAX, 0, p->t.cart_comm); CHKERR(ierr);
  ierr = MPI_Reduce(&(p->prof.others), &max_others, 1, MPI_DOUBLE, MPI_MAX, 0, p->t.cart_comm); CHKERR(ierr);
  ierr = MPI_Reduce(&(p->prof.total), &max_total, 1, MPI_DOUBLE, MPI_MAX, 0, p->t.cart_comm); CHKERR(ierr);

  ierr = MPI_Reduce(&(p->prof.compute), &min_comp, 1, MPI_DOUBLE, MPI_MIN, 0, p->t.cart_comm); CHKERR(ierr);
  ierr = MPI_Reduce(&(p->prof.communicate), &min_comm, 1, MPI_DOUBLE, MPI_MIN, 0, p->t.cart_comm); CHKERR(ierr);
  ierr = MPI_Reduce(&(p->prof.wait), &min_wait, 1, MPI_DOUBLE, MPI_MIN, 0, p->t.cart_comm); CHKERR(ierr);
  ierr = MPI_Reduce(&(p->prof.others), &min_others, 1, MPI_DOUBLE, MPI_MIN, 0, p->t.cart_comm); CHKERR(ierr);
  ierr = MPI_Reduce(&(p->prof.total), &min_total, 1, MPI_DOUBLE, MPI_MIN, 0, p->t.cart_comm); CHKERR(ierr);

  ierr = MPI_Reduce(&(p->prof.compute), &mean_comp, 1, MPI_DOUBLE, MPI_SUM, 0, p->t.cart_comm); CHKERR(ierr); mean_comp/=p->mpi_size;
  ierr = MPI_Reduce(&(p->prof.communicate), &mean_comm, 1, MPI_DOUBLE, MPI_SUM, 0, p->t.cart_comm); CHKERR(ierr); mean_comm/=p->mpi_size;
  ierr = MPI_Reduce(&(p->prof.wait), &mean_wait, 1, MPI_DOUBLE, MPI_SUM, 0, p->t.cart_comm); CHKERR(ierr); mean_wait/=p->mpi_size;
  ierr = MPI_Reduce(&(p->prof.others), &mean_others, 1, MPI_DOUBLE, MPI_SUM, 0, p->t.cart_comm); CHKERR(ierr); mean_others/=p->mpi_size;
  ierr = MPI_Reduce(&(p->prof.total), &mean_total, 1, MPI_DOUBLE, MPI_SUM, 0, p->t.cart_comm); CHKERR(ierr); mean_total/=p->mpi_size;


  ////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////
  // print the performance results
  if (p->mpi_rank == 0) {
    printf("Total memory allocation per MPI rank: %lu MiB\n", sizeof(FLOAT_PRECISION)*p->ln_domain*3/1024/1024);
    printf("Total time(s): %e\n", t*p->n_tests);
    printf("time/test(s): %e\n", t);

    if(p->target_ts != 2){
      printf("\nRANK0 MStencil/s  MIN: %f  \n", p->ln_stencils/(1e6*t_max));
      printf("RANK0 MStencil/s  AVG: %f  \n", p->ln_stencils/(1e6*(t/p->nt) ));
      printf("RANK0 MStencil/s  MAX: %f  \n", p->ln_stencils/(1e6*t_min));

      printf("\n******************************************************\n");
      printf("RANK0 Total: %f (s) -%06.2f%%\n", p->prof.total,p->prof.total/p->prof.total*100);
      printf("RANK0 Computation: %f (s) - %05.2f%%\n", p->prof.compute, p->prof.compute/p->prof.total*100);
      printf("RANK0 Communication: %f (s) - %05.2f%%\n", p->prof.communicate,p->prof.communicate/p->prof.total*100);
      printf("RANK0 Waiting: %f (s) - %05.2f%%\n", p->prof.wait,p->prof.wait/p->prof.total*100);
      printf("RANK0 Other: %f (s) - %05.2f%%\n", p->prof.others,p->prof.others/p->prof.total*100);
      printf("\n******************************************************\n");
      printf("MEAN Total: %f (s) -%06.2f%%\n", mean_total,mean_total/mean_total*100);
      printf("MEAN Computation: %f (s) - %05.2f%%\n", mean_comp, mean_comp/mean_total*100);
      printf("MEAN Communication: %f (s) - %05.2f%%\n", mean_comm,mean_comm/mean_total*100);
      printf("MEAN Waiting: %f (s) - %05.2f%%\n", mean_wait,mean_wait/mean_total*100);
      printf("MEAN Other: %f (s) - %05.2f%%\n", mean_others,mean_others/mean_total*100);

      printf("\n******************************************************\n");
      printf("MAX Total: %f (s)\n", max_total);
      printf("MAX Computation: %f (s)\n", max_comp);
      printf("MAX Communication: %f (s)\n", max_comm);
      printf("MAX Waiting: %f (s)\n", max_wait);
      printf("MAX Other: %f (s)\n", max_others);

      printf("\n******************************************************\n");
      printf("MIN Total: %f (s)\n", min_total);
      printf("MIN Computation: %f (s)\n", min_comp);
      printf("MIN Communication: %f (s)\n", min_comm);
      printf("MIN Waiting: %f (s)\n", min_wait);
      printf("MIN Other: %f (s)\n", min_others);
    }

    if( (p->wavefront !=0) && (p->target_ts == 2) ) {
//      printf("\nTotal RANK0 MStencil/s  MIN: %f  \n", p->ln_stencils/(1e6*t_max));
//      printf("Total RANK0 MStencil/s  MAX: %f  \n", p->ln_stencils/(1e6*t_min));

      printf("******************************************************\n");
      total_stencils = ((unsigned long) p->ln_stencils * (unsigned long) p->nt - p->idiamond_pro_epi_logue_updates)/(1e6);
      printf("MWD main-loop RANK0 MStencil/s  MIN: %f\n", total_stencils*1.0/(t_ts_main_max));
      printf("MWD main-loop RANK0 MStencil/s  MAX: %f\n", total_stencils*1.0/(t_ts_main_min));

      printf("******************************************************\n");
      printf("%-27s %f (s) - %05.2f%%\n", "RANK0 ts main loop:",
          p->prof.ts_main,p->prof.ts_main/p->prof.total*100);
      printf("%-27s %f (s) - %05.2f%%\n", "RANK0 ts prologue/epilogue:",
          p->prof.ts_others, p->prof.ts_others/p->prof.total*100);
      printf("%-27s %f (s) - %05.2f%%\n","RANK0 ts others:",
          p->prof.total-(p->prof.ts_main+p->prof.ts_others), (p->prof.total-(p->prof.ts_main+p->prof.ts_others))/p->prof.total*100);
      printf("******************************************************\n");
      printf("%-27s", "Metric \\ core:");
      for(i=0; i<p->num_threads; i++) printf("  core %02d     ", i);
      printf("\n");

      printf("%-27s", "Wavefront barrier wait [s]:");
      for(i=0; i<p->num_threads; i++) printf("  %e", p->stencil_ctx.t_wait[i]);
      printf("\n");
      printf("%-27s", "Wavefront barrier wait [%]:");
      for(i=0; i<p->num_threads; i++) printf("  %05.2f       ", p->stencil_ctx.t_wait[i]/p->prof.ts_main*100);
      printf("\n\n");

      printf("%-27s", "Metric \\ thread group:");
      for(i=0; i<num_thread_groups; i++) printf("  group %02d    ", i);
      printf("\n");

      printf("%-27s", "Wavefront steady state [s]:");
      for(i=0; i<num_thread_groups; i++) printf("  %e", p->stencil_ctx.t_wf_main[i]);
      printf("\n");
      printf("%-27s", "Wavefront steady state [%]:");
      for(i=0; i<num_thread_groups; i++) printf("  %05.2f       ", p->stencil_ctx.t_wf_main[i]/p->prof.ts_main*100);
      printf("\n");

      printf("%-27s", "Wavefront startup/end [s]:");
      for(i=0; i<num_thread_groups; i++) printf("  %e", p->stencil_ctx.t_wf_prologue[i] + p->stencil_ctx.t_wf_epilogue[i]);
      printf("\n");
      printf("%-27s", "Wavefront startup/end [%]:");
      for(i=0; i<num_thread_groups; i++) printf("  %05.2f       ",
          (p->stencil_ctx.t_wf_prologue[i] + p->stencil_ctx.t_wf_epilogue[i])/p->prof.ts_main*100);
      printf("\n");

      printf("%-27s", "Wavefront communication [s]:");
      for(i=0; i<num_thread_groups; i++) printf("  %e", p->stencil_ctx.t_wf_comm[i]);
      printf("\n");
      printf("%-27s", "Wavefront communication [%]:");
      for(i=0; i<num_thread_groups; i++) printf("  %05.2f       ", p->stencil_ctx.t_wf_comm[i]/p->prof.ts_main*100);
      printf("\n");

      printf("%-27s", "Wavefront others [s]:");
      for(i=0; i<num_thread_groups; i++) printf("  %e",
          p->prof.ts_main -(p->stencil_ctx.t_wf_main[i] + p->stencil_ctx.t_wf_prologue[i] + p->stencil_ctx.t_wf_comm[i] + p->stencil_ctx.t_wf_epilogue[i] + p->stencil_ctx.t_group_wait[i]));
      printf("\n");
      printf("%-27s", "Wavefront others [%]:");
      for(i=0; i<num_thread_groups; i++) printf("  %05.2f       ",
          (p->prof.ts_main - (p->stencil_ctx.t_wf_main[i] + p->stencil_ctx.t_wf_prologue[i] + p->stencil_ctx.t_wf_comm[i] + p->stencil_ctx.t_wf_epilogue[i]))/p->prof.ts_main*100);
      printf("\n");

      printf("%-27s", "Group spin-wait [s]:");
      for(i=0; i<num_thread_groups; i++) printf("  %e", p->stencil_ctx.t_group_wait[i]);
      printf("\n");
      printf("%-27s", "Group spinn-wait [%]:");
      for(i=0; i<num_thread_groups; i++) printf("  %05.2f       ", p->stencil_ctx.t_group_wait[i]/p->prof.ts_main*100);
      printf("\n");
 
      printf("%-27s", "Resolved diamonds:");
      for(i=0; i<num_thread_groups; i++) printf("  %e", p->stencil_ctx.wf_num_resolved_diamonds[i]);

    }
    printf("\n******************************************************\n");

  }
}

void print_param(Parameters p) {

  char *coeff_type, *precision;
  int wf_halo = p.stencil.r;
  int diam_height;

  diam_height = (p.t_dim*2)*p.stencil.r + p.stencil_ctx.num_wf;

  switch (p.stencil.coeff){
  case CONSTANT_COEFFICIENT:
    coeff_type = "constant";
    break;
  case VARIABLE_COEFFICIENT:
    coeff_type = "variable";
    break;
  case VARIABLE_COEFFICIENT_AXSYM:
    coeff_type = "variable axis-symmetric";
    break;
  case VARIABLE_COEFFICIENT_NOSYM:
    coeff_type = "variable no-symmetry";
    break;
  }

  precision = ((sizeof(FLOAT_PRECISION)==4)?"SP":"DP");

  printf("\n******************************************************\n");
  printf("Parameters settings\n");
  printf("******************************************************\n");
  printf("Time stepper name: %s\n", TSList[p.target_ts].name);
  printf("Stencil Kernel name: %s\n", p.stencil.name);
  printf("Stencil Kernel semi-bandwidth: %d\n", p.stencil.r);
  printf("Stencil Kernel coefficients: %s\n", coeff_type);
  printf("Precision: %s\n", precision);
  printf("Global domain    size:%lu    nx:%d    ny:%d    nz:%d\n", p.n_stencils, p.stencil_shape[0],p.stencil_shape[1],p.stencil_shape[2]);
  printf("Rank 0 domain    size:%lu    nx:%d    ny:%d    nz:%d\n", p.ln_stencils, p.lstencil_shape[0],p.lstencil_shape[1],p.lstencil_shape[2]);
  printf("Number of time steps: %d\n", p.nt);
  printf("Alignment size: %d Bytes\n", p.alignment);
  printf("Number of tests: %d\n", p.n_tests);
  printf("Verify:   %d\n", p.verify);
  printf("Source point enabled: %d\n", p.source_point_enabled);
  printf("Time unroll:   %d\n", p.t_dim);
  printf("Using separate call to central line update: %d\n", USE_SPLIT_STRIDE);
  printf("Halo concatenation: %d\n", p.halo_concat);
  // Print kernel specific parameters
  switch(p.target_ts){
  case 0:
  case 1:
    if(p.h[2].is_contiguous ==1) printf("MPI datatype is contiguous across the Z direction\n");
    printf("Block size in Y: %d\n", p.stencil_ctx.bs_y);
    printf("OpenMP schedule: ");
    if(p.use_omp_stat_sched==1) printf("static\n"); else printf("static1\n");
    break;
  case 2: // dynamic scheduling intra diamond methods
    printf("Block size in X: %d\n", p.stencil_ctx.bs_x);
    printf("Enable wavefronts: %d\n", p.wavefront!=0);
    if(p.stencil_ctx.thread_group_size!=1) printf("Wavefront parallel strategy: %s\n", MWD_name[p.mwd_type]);
    printf("Intra-diamond width:   %d\n", (p.t_dim+1)*2*p.stencil.r);
    printf("Wavefront width:  %d\n", diam_height);
    printf("Cache block size/wf (kiB): %lu\n", p.wf_blk_size/1024);
    printf("Total cache block size (kiB): %lu\n", (p.num_threads/p.stencil_ctx.thread_group_size) * p.wf_blk_size/1024);
    printf("Next larger cache block size/wf (kiB): %lu (diam_width=%d)\n", p.wf_larger_blk_size/1024, (p.larger_t_dim+1)*2*p.stencil.r);
    printf("Intra-diamond prologue/epilogue MStencils: %lu\n", p.idiamond_pro_epi_logue_updates/(1000*1000));
    printf("Multi-wavefront updates: %d\n", p.stencil_ctx.num_wf);
    printf("Thread group size: %d\n", p.stencil_ctx.thread_group_size);
    break;
  }

  printf("OpenMP Threads: %d\n", p.num_threads);
  printf("Assumed usable cache size: %dKiB\n",   p.cache_size);
  printf("MPI size: %d\n", p.mpi_size);
  printf("Processors topology (npx, npy, npz): %02d,%02d,%02d\n", p.t.shape[0], p.t.shape[1], p.t.shape[2]);
  printf("******************************************************\n");
}

void list_kernels(Parameters *p){
  int i;
  char *coeff_type;

  if(p->mpi_rank==0) {
    printf("Available time steppers:\n#    Name\n");
    i = 0;
    while(1){
      if (TSList[i].name == 0) break;
      printf("%02d   %s\n",i, TSList[i].name);
      i++;
    }

    printf("\nAvailable stencil kernels:\n");
    i = 0;
    while(1){
      switch (stencil_info_list[i].coeff){
      case CONSTANT_COEFFICIENT:
        coeff_type = "constant";
        break;
      case VARIABLE_COEFFICIENT:
        coeff_type = "variable";
        break;
      case VARIABLE_COEFFICIENT_AXSYM:
        coeff_type = "variable axis-symmetric";
        break;
      case VARIABLE_COEFFICIENT_NOSYM:
        coeff_type = "variable no-symmetry";
        break;
      }
      if (stencil_info_list[i].name == 0) break;
      printf("%02d  stencil_op:%s  time-order:%d  radius:%d  coeff:%s\n",
          i, stencil_info_list[i].name, stencil_info_list[i].time_order, stencil_info_list[i].r,coeff_type);
      i++;
    }

    printf("\nAvailable MWD implementations:\n#    Name\n");
    i = 0;
    while(1){
      if (MWD_name[i] == 0) break;
      printf("%02d   %s\n",i, MWD_name[i]);
      i++;
    }

  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(0);
}

void print_help(Parameters *p){
  if(p->mpi_rank == 0){
    printf(
        "Note: default values are available at 'param_default()' in src/utils.c\n"
        "Usage:\n\n"

        "  --help\n"
        "       Show available options\n"
        "  --list\n"
        "       List the available time stepper kernels and their number.\n"
        "       Description of the time steppers is available at wrappers.h\n"
        "  --verify <bool>\n"
        "       Verify the correctness of the selected time stepper using\n"
        "       the provided domain and topology information\n"
        "       (This option disables performance measurement experiments)\n"
        "\nGeneral experiment parameters:\n"
        "  --target-ts <integer>\n"
        "       Select the desired time stepper kernel\n"
        "       (can acquire the available time stepper through --list option)\n"
//        "  --stencil_radius <integer>\n"
//        "       Select the radius of the stencil operator (maximum is 10)\n"
//        "  --const-coeff-stencil <bool>  (default 1)"
//        "       Determine whether the stencil coefficients are constant (1) or variable (2)\n"
        "  --target-kernel <integer>\n"
        "       Select the desired stencil kernel\n"
        "       (can acquire the available kernels through --list option)\n"
        "  --nx <integer>\n"
        "       Global domain size in the x direction\n"
        "  --ny <integer>\n"
        "       Global domain size in the y direction\n"
        "  --nz <integer>\n"
        "       Global domain size in the z direction\n"
        "  --nt <integer>\n"
        "       The desired number of time steps\n"
        "  --npx <integer>\n"
        "       Number of MPI processes across the x direction\n"
        "  --npy <integer>\n"
        "       Number of MPI processes across the y direction\n"
        "  --npz <integer>\n"
        "       Number of MPI processes across the z direction\n"
        "  --n-tests <integer>\n"
        "       Specify how many times the time stepper is run at the\n"
        "       performance measurement experiments\n"
        "  --alignment <integer>\n"
        "       Specify the memory alignment value (in bytes) of the\n"
        "       allocated domain arrays\n"
//        "  --disable-source-point\n"
//        "       disables the source point update in the solution domain\n"

        "\nDisplay options:\n"
        "  --verbose <bool>\n"
        "       Show the used configuration information\n"
        "  --debug <bool>\n"
        "       Print detailed debugging information\n"

        "\nSpecialized arguments:\n"
        "  --t-dim <integer>    (specific to the Diamond methods)\n"
        "       Specifies the value of unrolling in time\n"
        "  --z-mpi-contig <bool>    (Specific to standard methods: 0-1)\n"
        "       Uses contiguous MPI datatype in the z direction of the MPI\n"
        "       topology at the standard methods\n"
        "  --halo-concatenate <integer>  (experimental feature)\n"
        "       Explicitly concatenate halo information before send and unpack them\n"
        "       at receiving end using multi-threading.\n"
        "       (Works for decomposition across X and Y only)\n"
        "  --thread-group-size <integer>   (specific to method 9)\n"
        "       Set the thread group size for methods supporting multiple thread groups\n"
        "  --cache-size <integer>\n"
        "       The usable last level cache size for spatial blocking\n"
        "  --wavefront <bool>\n"
        "       Indicate whether to use wavefront in the tile. 1 for yes and 0 for no (default 1)\n"
        "  --num-wavefronts <int>\n"
        "       Set the number of wavefronts updated per wavefront iteration (default 1)\n"
        "  --bsx <int>\n"
        "       Set block size in X for MWD (default 1e6)\n"
        "  --mwd-type <int>\n"
        "       Select one of the MWD implementations from the one available at --list option\n"
        "  --use-omp-stat-sched\n"
        "       Use OpenMP static schedule instead of static,1 at the spatial blocking time steppers\n"




//        "  --target-parallel-wavefront <integer>\n"
//        "       Indicate specify multi-threaded wavefront parallelization strategy (specific to ts 9)\n"
    );
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  exit(0);
}
void parse_args (int argc, char** argv, Parameters * p)
{ // for more details see http://libslack.org/manpages/getopt.3.html
  int c;
  int thread_group_size = -1;

  while (1)
  {
    int option_index = 0;
    static struct option long_options[] =
    {
        {"nz", 1, 0, 0},
        {"ny", 1, 0, 0},
        {"nx", 1, 0, 0},
        {"nt", 1, 0, 0},
        {"alignment", 1, 0, 0},
        {"verbose", 1, 0, 0},
        {"debug", 1, 0, 0},
        {"target-ts", 1, 0, 0},
        {"target-kernel", 1, 0, 0},
        {"n-tests", 1, 0, 0},
        {"verify", 1, 0, 0},
        {"list", 0, 0, 0},
        {"help", 0, 0, 0},
        {"npx", 1, 0, 0},
        {"npy", 1, 0, 0},
        {"npz", 1, 0, 0},
        {"t-dim", 1, 0, 0},
        {"z-mpi-contig", 1, 0, 0},
        {"disable-source-point", 0, 0, 0},
        {"halo-concatenate", 1, 0, 0},
        {"thread-group-size", 1, 0, 0},
        {"cache-size", 1, 0, 0},
        {"wavefront", 1, 0, 0},
        {"num-wavefronts", 1, 0, 0},
        {"pad-array", 0, 0, 0},
        {"bsx", 1, 0, 0},
        {"mwd-type", 1, 0, 0},
        {"use-omp-stat-sched", 0, 0, 0},
//        {"target-parallel-wavefront", 1, 0, 0},
        {0, 0, 0, 0}
    };

    c = getopt_long (argc, argv, "",
        long_options, &option_index);
    if (c == -1)
      break;

    switch (c)
    {
    case 0:
      if(strcmp(long_options[option_index].name, "nz") == 0) p->stencil_shape[2] = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "ny") == 0) p->stencil_shape[1] = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "nx") == 0) p->stencil_shape[0] = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "nt") == 0) p->nt = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "npx") == 0) p->t.shape[0] = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "npy") == 0) p->t.shape[1] = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "npz") == 0) p->t.shape[2] = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "alignment") == 0) p->alignment = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "verbose") == 0) p->verbose = atoi(optarg)!=0;
      else if(strcmp(long_options[option_index].name, "target-ts") == 0) p->target_ts = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "target-kernel") == 0) p->target_kernel = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "n-tests") == 0) p->n_tests = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "verify") == 0) p->verify = atoi(optarg)!=0;
      else if(strcmp(long_options[option_index].name, "debug") == 0) p->debug = atoi(optarg)!=0;
      else if(strcmp(long_options[option_index].name, "t-dim") == 0) p->t_dim = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "z-mpi-contig") == 0) p->h[2].is_contiguous = atoi(optarg)!=0;
      else if(strcmp(long_options[option_index].name, "list") == 0) list_kernels(p);
      else if(strcmp(long_options[option_index].name, "help") == 0) print_help(p);
      else if(strcmp(long_options[option_index].name, "disable-source-point") == 0) p->source_point_enabled=0;
      else if(strcmp(long_options[option_index].name, "halo-concatenate") == 0) p->halo_concat = atoi(optarg)!=0;
      else if(strcmp(long_options[option_index].name, "thread-group-size") == 0) thread_group_size = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "cache-size") == 0) p->cache_size = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "wavefront") == 0) p->wavefront = atoi(optarg)!=0;
      else if(strcmp(long_options[option_index].name, "num-wavefronts") == 0) p->stencil_ctx.num_wf = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "pad-array") == 0) p->array_padding = 1;
      else if(strcmp(long_options[option_index].name, "bsx") == 0) p->stencil_ctx.bs_x = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "mwd-type") == 0) p->mwd_type = atoi(optarg);
      else if(strcmp(long_options[option_index].name, "use-omp-stat-sched") == 0) p->use_omp_stat_sched = 1;
//      else if(strcmp(long_options[option_index].name, "target-parallel-wavefront") == 0) p->target_parallel_wavefront = atoi(optarg);
     break;

    default:
      if(p->mpi_rank == 0){
        fprintf(stderr, "Invalid arguments\n\n");
      }
      print_help(p);
      break;
    }
  }

  if (optind < argc)
  {
    if(p->mpi_rank == 0){
      fprintf(stderr, "Invalid arguments\n\n");
    }
    print_help(p);
  }

  set_centered_source(p);

  // allow thread group size change only for methods supporting multiple thread groups
  if(p->target_ts == 2) {
    p->stencil_ctx.thread_group_size = thread_group_size;
    p->use_omp_stat_sched = 0;
  }
}

void print_3Darray(char *filename, FLOAT_PRECISION * restrict array, int nx, int ny, int nz, int halo) {
  int i,j,k;
  FILE *fp;
  if(!(fp = fopen(filename, "w")))
    printf("ERROR: cannot open file for writing\n");

  for (i=halo; i<nz-halo; i++) {
    // new page
    fprintf(fp, "\n***************** slice # %d *************\n",i+1-halo);
    for (j=halo; j<ny-halo; j++) {
      for (k=halo; k<nx-halo; k++) {
        fprintf(fp, "%+.4e ", array[(i*ny+j)*nx+k]);
      }
      // new line
      fprintf(fp, "\n");
    }
  }

  fclose(fp);
}
