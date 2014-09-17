#!/bin/bash

#SBATCH -J halo_first           # Job name
#SBATCH -o halo_first.o%j       # Name of stdout output file(%j expands to jobId)
#SBATCH -e halo_first.o%j       # Name of stderr output file(%j expands to jobId)
#SBATCH -p development    # Submit to the 'normal' or 'development' queue
#SBATCH -N 16             # Total number of nodes
#SBATCH -n 32             # Total number of mpi tasks requested
#SBATCH -t 00:15:00       # Run time (hh:mm:ss)

# Set the number of threads per task(Default=1)
export OMP_NUM_THREADS=7

export MPICH_ASYNC_PROGRESS=1 

# Run the hybrid application 
ibrun tacc_affinity $HOME/mwd_kernels/build/mwd_kernel --disable-source-point --target-ts 2 --npx 1 --npy 1 --npz 32 --nx 512 --ny 512 --nz 8192 --nt 802 --z-mpi-contig 1
