#!/bin/bash
#SBATCH -J test           # Job name
#SBATCH -o test.o%j       # Name of stdout output file(%j expands to jobId)
#SBATCH -e test.o%j       # Name of stderr output file(%j expands to jobId)
#SBATCH -p development    # Submit to the 'normal' or 'development' queue
#SBATCH -N 1              # Total number of nodes requested (16 cores/node)
#SBATCH -n 2              # Total number of mpi tasks requested
#SBATCH -t 00:05:00       # Run time (hh:mm:ss)

# Set the number of threads per task(Default=1)
export OMP_NUM_THREADS=8

# Run the hybrid application 
ibrun $HOME/mwd_kernels/build/mwd_kernel --npz 2 --nt 100 --nx 64 --ny 64 --nz 64 --target-ts 1
