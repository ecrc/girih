#!/bin/bash
#SBATCH --account=k1137
#SBATCH --job-name=bench-girih
#SBATCH --output=bench-girih-%j.out
#SBATCH --error=bench-girih-%j.err
#SBATCH --nodes=1
#SBATCH --time=00:30:00

module swap PrgEnv-cray PrgEnv-intel
#module load boost/1.58

#export LD_LIBRARY_PATH=/scratch/x_etiennv":"$LD_LIBRARY_PATH

#export KMP_AFFINITY=scatter,1,0,granularity=fine 

#export OMP_NUM_THREADS=16
#srun -o 0 --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=3
#2 --ntasks-per-socket=16 --ntasks-per-core=1 --cpu_bind=cores /project/k1137/sof
#t/bin/gptexe2 --xml=modelling.xml --cbx=14 --cby=2 --nt_user=10 

#cd /scratch/x_tonelltl/runs/prod1

#srun -o 0 --ntasks=256 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=32 --ntasks-per-socket=16 --ntasks-per-core=1 --cpu_bind=cores /scratch/x_etiennv/gptexe2 --xml=modelling.xml --cbx=14 --cby=2  --tmax=24

#srun -o 0 --ntasks=256 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=32 --ntasks-per-socket=16 --ntasks-per-core=1\
# --cpu_bind=cores /scratch/x_etiennv/gptexe2 --xml=modelling.xml --cbx=14 --cby=2

#OMP_NUM_THREADS=32 srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --nx 512 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2
# --thread-group-size 8
#OMP_NUM_THREADS=32 srun --ntasks=2 --cpus-per-task=16 --ntasks-per-node=2 ./build/mwd_kernel --npx 1 --npy 2 --npz 1 --nx 512 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2
#OMP_NUM_THREADS=32 srun --ntasks=2 --cpus-per-task=16 --hint=nomultithread --ntasks-per-node=2 ./build/mwd_kernel --npx 1 --npy 2 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 4
export OMP_NUM_THREADS=32

gs=501; nt=2000
srun ./build/mwd_kernel --nx $gs  --ny $gs --nz $gs --nt $nt --target-kernel 0 --mwd-type 0 --target-ts 0 --verify 1 $source 
exit 0
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
