/** 
 * @copyright (c) 2014- King Abdullah University of Science and
 *                      Technology (KAUST). All rights reserved.
 **/
 

/**
 * @file src/wrappers.h 

 * GIRIH is a high performance stencil framework using wavefront 
 * 	diamond tiling.
 * GIRIH is provided by KAUST.
 *
 * @version 1.0.0
 * @author Tareq Malas
 * @date 2017-11-13
 **/

#ifndef WRAPPERS_H_
#define WRAPPERS_H_

#include "data_structures.h"

extern void naive_nonblocking_ts(Parameters *p);
extern void halo_first_ts(Parameters *p);
extern void dynamic_intra_diamond_ts(Parameters *p);

struct time_stepper TSList[] = {
/* 0 */    {"Spatial Blocking", naive_nonblocking_ts},
    // Naive Halo exchange implementation with non-locking MPI communication
/* 1 */    {"Halo-first", halo_first_ts},
    // Compute the sides then communicate the halo while computing the middle
/* 2 */    {"Diamond", dynamic_intra_diamond_ts},
    // Dynamic scheduling Intra-Diamond method implementation with no communication computation overlapping
/* end */    {0, 0},
};

#endif /* WRAPPERS_H_ */
