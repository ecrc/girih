/** 
 * @copyright (c) 2014- King Abdullah University of Science and
 *                      Technology (KAUST). All rights reserved.
 **/
 

/**
 * @file src/driver.h 

 * GIRIH is a high performance stencil framework using wavefront 
 * 	diamond tiling.
 * GIRIH is provided by KAUST.
 *
 * @version 1.0.0
 * @author Tareq Malas
 * @date 2017-11-13
 **/

#ifndef DRIVER_H_
#define DRIVER_H_

#define _XOPEN_SOURCE 600       /* make stdlib.h provide posix_memalign() in strict C99 mode */

#include "data_structures.h"
#include <stdio.h>
#include "mpi.h"
#include "wrappers.h"

extern void performance_test(Parameters * p);
extern void param_default(Parameters *);
extern void parse_args (int argc, char** argv, Parameters *);
extern void mpi_topology_init(Parameters *);
extern void init(Parameters *);
extern void verify(Parameters *);

#endif /* DRIVER_H_ */
