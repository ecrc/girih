#!/bin/sh 

. /home/hpc/unrz/unrz280/girih/scripts/mic/run_test.sh

#k=1
#N=768
k=5
N=456


ts=0

# 1 thread per core tests
res_dir=$HOME/girih/results/mic_thread_scaling/
mkdir -p $res_dir

for blk_size in 1 2 4
do
  for group in MEM1 MEM2 MEM3 MEM4 MEM5 CPI
  do

    t=`expr $blk_size \* 6`
    te=`expr $blk_size \* 60`
    ccs=`expr 256 / $blk_size`

    while [ $t -le $te ]
    do
      if [ $ts -eq 0 ]
      then
        cs=`expr $t \* $ccs`
        tgs=0
        mwd_t=0
      else
        cs=64
        tgs=$blk_size
        mwd_t=2
      fi

      affinity="E:N:${t}:${blk_size}:4"
      nt=`expr $t \* 1`
      if [ $nt -gt 50 ]
      then
        nt=50
      elif [ $nt -lt 4 ]
      then 
        nt=4
      fi 

      run_test $t $affinity $group $N $ts $k $nt $tgs $cs $mwd_t

      t=`expr $t + $blk_size \* 6`

    done
  done
done



