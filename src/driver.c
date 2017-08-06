#include "driver.h"
#include <stdlib.h> //@KADIR for exit();

int main(int argc, char** argv)
{

  LIKWID_MARKER_INIT;

  int provided, flag, claimed;
  MPI_Init_thread( &argc, &argv, MPI_THREAD_MULTIPLE, &provided );
//  MPI_Is_thread_main( &flag );
//  if (!flag)
//      printf( "This thread called init_thread but Is_thread_main gave false\n" );fflush(stdout);
  MPI_Query_thread( &claimed );
  if (claimed != provided)
      printf( "Query thread gave thread level %d but Init_thread gave %d\n", claimed, provided );fflush(stdout);

//  MPI_Init(&argc,&argv);
  int mpi_rank, mpi_size;
  MPI_Comm_rank (MPI_COMM_WORLD, &(mpi_rank));
  MPI_Comm_size (MPI_COMM_WORLD, &(mpi_size));

  // configure the experiment's parameters
  Parameters p;
  p.mpi_rank = mpi_rank;
  p.mpi_size = mpi_size;
  param_default(&p);

  parse_args(argc, argv, &p);

#pragma omp parallel num_threads(p.num_threads)
  {
    LIKWID_MARKER_THREADINIT;
  }


  // Simple error checking
  if (p.t.shape[0]*p.t.shape[1]*p.t.shape[2] != p.mpi_size) {
    if(p.mpi_rank==0) fprintf(stderr,"ERROR: requested MPI topology shape does not match the available processes count: \n\tRequested:%03d \n\tAvailable:%03d\n",
        p.t.shape[0]*p.t.shape[1]*p.t.shape[2], p.mpi_size);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 1;
  }

  // Create the MPI topology
  mpi_topology_init(&p);
  //  MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  //printf("%s %d\tKADIR: ENABLE AUTO TUNING\n", __FILE__, __LINE__);
  // initialize time-stepper specific requirements
  init(&p);

  p.source_point_enabled = 1; //@KADIR DIAMOND DISABLES SOURCE_POINT_ENABLED

  // Verify the time stepper kernel if required
  if (p.verify !=0) {
    verify(&p);
  } else { // do performance tests
    performance_test(&p);
  }

  MPI_Finalize();

  LIKWID_MARKER_CLOSE;

  return 0;
}
