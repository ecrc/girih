#!/bin/bash
#SBATCH --account=k1137
#SBATCH --job-name=girih
#SBATCH --output=/project/k1137/akbudak/girih/out/%j
#SBATCH --error=/project/k1137/akbudak/girih/out/%j
#SBATCH --nodes=1
##SBATCH --time=12:00:00
#SBATCH --time=00:05:00 #TODO
#SBATCH --partition=debug

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
    sruncmd="srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 " # -vvvvvvv " 
else
    dir=/home/akbudak/girih
    sruncmd=""
    if [[ "$hn" == "uwork" ]]; then
        dir=/home/kadir/girih
    fi
fi


gs=512; nt=2101;verify=0;ntests=1  #exawave
gs=512; nt=506;verify=0;ntests=1  #exawave
gs=80; nt=250;verify=1;ntests=1  #debug
#gs=512; nt=506; num_threads="1 2 4 8 12 16 20 24 28 32";verify=0;ntests=1   #eage paper
#gs=512; nt=506; num_threads="20 24 28 32";verify=0;ntests=2   #eage paper
#shaheen     0 1 2 3 4 5  6  7  8  9  10 11 12 13 14  15 16 17 18 19 20 21 22
num_threads=(1 2 4 8 8 16 16 16 20 20 24 24 24 24 24  28 28 28 28 32 32 32 32)
tgss=(       1 2 2 2 4 2  4  8  2  4  2  4  6  8  12  2  4  7  14 2  4  8  16)
nexp=22

#shihab      0  1  2  3  4  5  6  7  8  9  10 11 12 13 14  15 16 17 18 19 20 21 22
#num_threads=(1  3  3  6  6  9  9  18 18 18 18 36 36 36  36 36 36 36            )
#tgss=(       1  1  3  2  3  1  3  2  3  6  9  2  3  4   6  9  12 18            )
#nexp=17

#jasmine     0  1  2  3  4  5  6  7  8  9  10 11 12 13 14  15 16 17 18 19 20 21 22
#num_threads=(1  4  6  6  8  8  12 12 12 16 16 16 20  20 20 20 24 24 24 24 24            )
#tgss=(       1  2  2  3  2  4  3  4  6  2  4  8  2   4  5  10 2  4  6  8  12            )
#nexp=20
#for i in `seq 1 $nexp`;do
for i in 4;do
    nthread=${num_threads[i]}
    tgs=${tgss[i]}
    export OMP_NUM_THREADS=$nthread
    #export KMP_AFFINITY=verbose
    call_combined_function="--call_combined_function"
    #DIAMOND
    cmd=$dir"/build/mwd_kernel --nx $gs --ny $gs --nz $gs --nt $nt --mwd-type 1 --target-ts 2 --verify $verify --npx 1 --npy 1 --npz 1 --threads $nthread  --n-tests $ntests --thread-group-size $tgs --target-kernel 0   --thx 1 --thy 2 --thz 2 --thc 1 $call_combined_function"
    #cmd=$dir"/build/mwd_kernel --nx $gs --ny $gs --nz $gs --nt $nt --mwd-type 1 --target-ts 2 --verify $verify --npx 1 --npy 1 --npz 1 --threads $nthread  --n-tests $ntests --thread-group-size $tgs "
    echo $sruncmd $cmd;$sruncmd $cmd 

    #cmd=$dir"/build/mwd_kernel --nx $gs --ny $gs --nz $gs --nt $nt --mwd-type 1 --target-ts 2 --verify $verify --npx 1 --npy 1 --npz 1  --n-tests $ntests --thread-group-size $tgs --target-kernel 0   --thx 1 --thy 2 --thz 2 --thc 1"
    #echo $sruncmd $cmd;$sruncmd $cmd 

    #export KMP_AFFINITY=verbose,"granularity=fine,proclist=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31],explicit"
    #echo $sruncmd $cmd;$sruncmd $cmd 

    #export KMP_AFFINITY=verbose,"granularity=fine,proclist=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31],explicit"
    #echo $sruncmd $cmd;$sruncmd $cmd 
    
    #export KMP_AFFINITY=verbose,"granularity=fine,proclist=[0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62],explicit"
    #sruncmd="srun --ntasks=1 --cpus-per-task=64 --ntasks-per-node=1 " # -vvvvvvv " 
    #echo $sruncmd $cmd;$sruncmd $cmd 

    #DIAMOND
    #srun  --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx $gs --ny $gs --nz $gs --nt $nt --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
done
exit 0
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
srun --ntasks=1 --cpus-per-task=32 --hint=nomultithread --ntasks-per-node=1 ./build/mwd_kernel --npx 1 --npy 1 --npz 1 --nx 513 --ny 512 --nz 512 --nt 506 --target-kernel 0 --mwd-type 1 --target-ts 2 --thread-group-size 8 --thx 1 --thy 1 --thz 8 
