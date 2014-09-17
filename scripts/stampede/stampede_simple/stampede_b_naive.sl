#!/bin/bash
#SBATCH -J blocking_naive           # Job name
#SBATCH -o blocking_naive.o%j       # Name of stdout output file(%j expands to jobId)
#SBATCH -e blocking_naive.o%j       # Name of stderr output file(%j expands to jobId)
#SBATCH -p development    # Submit to the 'normal' or 'development' queue
#SBATCH -N 1              # Total number of nodes
#SBATCH -n 2              # Total number of mpi tasks requested
#SBATCH -t 00:59:00       # Run time (hh:mm:ss)
#SBATCH -A TG-ASC130040

# Set the number of threads per task(Default=1)
export OMP_NUM_THREADS=8

# Run the hybrid application
ibrun -o 0 -n 1 tacc_affinity $HOME/mwd_kernels/build/mwd_kernel --target-ts 0 --n-tests 2 --disable-source-point --halo-concatenate 0 --npx 1 --npy 1 --npz 1 --nx 512 --ny 512 --nz 512 --nt 802
