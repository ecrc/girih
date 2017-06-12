#!/bin/bash
#SBATCH --account=k1137
#SBATCH --job-name=girih
#SBATCH --output=/project/k1137/akbudak/girih/out/%j
#SBATCH --error=/project/k1137/akbudak/girih/out/%j
#SBATCH --nodes=1
#SBATCH --time=12:00:00

hn=`hostname`
echo $hn
sha=0
if [[ "$hn" == nid0* ]]; then 
    echo shaheen; 
    sha=1
fi
if [ $sha -eq 1 ]; then
    module list
    module swap PrgEnv-cray PrgEnv-intel
    module list
    dir=/lustre/project/k1137/akbudak/girih
    #dir=/lustre/project/k1137/akbudak/ogirih  #FIXME
    sruncmd="srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 " 
else
    dir=/home/akbudak/girih
    sruncmd=""
fi


gs=512; nt=2101;verify=0;ntests=1  #exawave
#gs=512; nt=506; num_threads="1 2 4 8 12 16 20 24 28 32";verify=0;ntests=2   #eage paper
#gs=512; nt=506; num_threads="20 24 28 32";verify=0;ntests=2   #eage paper
nexp=22
#shaheen     0 1 2 3 4 5  6  7  8  9  10 11 12 13 14  15 16 17 18 19 20 21 22
num_threads=(1 2 4 8 8 16 16 16 20 20 24 24 24 24 24  28 28 28 28 32 32 32 32)
tgss=(       1 2 2 2 4 2  4  8  2  4  2  4  6  8  12  2  4  7  14 2  4  8  16)

nexp=17
#shihab      0  1  2  3  4  5  6  7  8  9  10 11 12 13 14  15 16 17 18 19 20 21 22
num_threads=(1  3  3  6  6  9  9  18 18 18 18 36 36 36  36 36 36 36            )
tgss=(       1  1  3  2  3  1  3  2  3  6  9  2  3  4   6  9  12 18            )

nexp=20
#jasmine     0  1  2  3  4  5  6  7  8  9  10 11 12 13 14  15 16 17 18 19 20 21 22
num_threads=(1  4  6  6  8  8  12 12 12 16 16 16 20  20 20 20 24 24 24 24 24            )
tgss=(       1  2  2  3  2  4  3  4  6  2  4  8  2   4  5  10 2  4  6  8  12            )
#for i in `seq 0 $nexp`;do
for i in 20;do
    nthread=${num_threads[i]}
    tgs=${tgss[i]}
    export OMP_NUM_THREADS=$nthread
    #DIAMOND
    #cmd=$dir"/build/mwd_kernel --nx $gs --ny $gs --nz $gs --nt $nt --mwd-type 1 --target-ts 2 --verify $verify --npx 1 --npy 1 --npz 1 --threads $nthread  --n-tests $ntests --thread-group-size $tgs --target-kernel 0   --thx 1 --thy 2 --thz 2 --thc 1"
    cmd=$dir"/build/mwd_kernel --nx $gs --ny $gs --nz $gs --nt $nt --mwd-type 1 --target-ts 2 --verify $verify --npx 1 --npy 1 --npz 1 --threads $nthread  --n-tests $ntests --thread-group-size $tgs "
    echo $sruncmd $cmd
    $sruncmd $cmd 
    #DIAMOND
    #srun  --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx $gs --ny $gs --nz $gs --nt $nt --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
done
exit 0
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
