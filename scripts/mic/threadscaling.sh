#!/bin/sh 
export LIKWID_LIB=-L/home/hpc/unrz/unrz280/local.mic/lib && export LIKWID_INC=-I/home/hpc/unrz/unrz280/local.mic/include
export PATH=$PATH:~/local.mic/bin
res_dir=$HOME/girih/results/mic_thread_scaling/
mkdir -p $res_dir

run_exp() {

  res_file="${res_dir}kernel${k}_tgs${tgs}_ts${ts}_th${t}_grp_${group}_${affinity}_affinity.txt"

  thx=`expr $t - 1`

  export OMP_NUM_THREADS=$t; export KMP_AFFINITY=${affinity}; likwid-perfctr -m -c 0-$thx -g $group -- /home/hpc/unrz/unrz280/girih/build_dp/mwd_kernel --n-tests 2 --nx $N --ny $N --nz $N --target-ts $ts --target-kernel $k --nt $nt --use-omp-stat-sched --cache-size 4096 --thread-group-size $tgs --cache-size $cs | tee -a ${res_file} 

}


k=1
N=480
ts=0
nt=100
tgs=0
affinity=balanced

for group in MEM1 MEM2 MEM3 MEM4 MEM5 CPI
do

  t=1
  te=60
  ccs=256
  while [ $t -le $te ]
  do
    cs=`expr $t \* $ccs`
    run_exp
    echo $t
    t=`expr $t + 1`
  done

  t=61
  te=120
  ccs=128
  while [ $t -le $te ]
  do
    cs=`expr $t \* $ccs`
    run_exp
    echo $t
    t=`expr $t + 1`
  done


done
