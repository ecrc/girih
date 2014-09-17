#!/bin/bash
#SBATCH -J naive           # Job name
#SBATCH -o naive.o%j       # Name of stdout output file(%j expands to jobId)
#SBATCH -e naive.o%j       # Name of stderr output file(%j expands to jobId)
#SBATCH -p development    # Submit to the 'normal' or 'development' queue
#SBATCH -N 1              # Total number of nodes
#SBATCH -n 2              # Total number of mpi tasks requested
#SBATCH -t 00:15:00       # Run time (hh:mm:ss)

# Set the number of threads per task(Default=1)
export OMP_NUM_THREADS=8
export OFFLOAD_REPORT=2

# Run the hybrid application 
ibrun tacc_affinity $HOME/mwd_kernels/build/mwd_kernel --alignment 64 --n-tests 2 --disable-source-point --target-ts 1 --npx 1 --npy 2 --npz 1 --nx 512 --ny 512 --nz 1024 --nt 4
