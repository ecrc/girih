#!/bin/bash

#SBATCH -J intra_diamond           # Job name
#SBATCH -o intra_diamond.o%j       # Name of stdout output file(%j expands to jobId)
#SBATCH -e intra_diamond.o%j       # Name of stderr output file(%j expands to jobId)
#SBATCH -p development    # Submit to the 'normal' or 'development' queue
#SBATCH -N 1              # Total number of nodes
#SBATCH -n 2              # Total number of mpi tasks requested
#SBATCH -t 00:15:00       # Run time (hh:mm:ss)

# Set the number of threads per task(Default=1)
export OMP_NUM_THREADS=8

# Run the hybrid application 
ibrun tacc_affinity $HOME/mwd_kernels/build/mwd_kernel --n-tests 2 --target-ts 6 --npx 1 --npy 2 --npz 1 --nx 512 --ny 4096 --nz 512 --nt 802 --t-dim 1 --halo-concatenate 1
