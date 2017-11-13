/** 
 * @copyright (c) 2014- King Abdullah University of Science and
 *                      Technology (KAUST). All rights reserved.
 **/
 

/**
 * @file src/kernels/solar_kernels.h 

 * GIRIH is a high performance stencil framework using wavefront 
 * 	diamond tiling.
 * GIRIH is provided by KAUST.
 *
 * @version 1.0.0
 * @author Tareq Malas
 * @date 2017-11-13
 **/

#include "solar_spt_blk.ic"
#include "solar_spt_blk_stat_sched.ic"

#include "solar_naive.ic" // For 1WD implementation
#include "solar_par_components.ic" // OpenMP parallel componets
#include "solar_1wf.ic"
//#include "solar_mwf.ic"
#include "solar_femwf.ic"
//#include "solar_rsmwf.ic"
//#include "solar_rsfemwf.ic"

