run_exp() {


  t=$1
  top=$2
  aff=$3
  group=$4
  N=$5
  ts=$6
  k=$7
  nt=$8
  tgs=$9
  cs=$10

  res_file="${res_dir}kernel${k}_tgs${tgs}_ts${ts}_th${t}_thpgrp${top}_grp_${group}_${aff}_affinity.txt"

  thx=`expr $t - 1`

  echo export OMP_NUM_THREADS=$t ';' export KMP_PLACE_THREADS=$top ';' export KMP_AFFINITY=${aff},verbose ';' likwid-perfctr -m -c 0-$thx -g $group -- /home/hpc/unrz/unrz280/girih/build_dp/mwd_kernel --n-tests 4 --nx $N --ny $N --nz $N --target-ts $ts --target-kernel $k --nt $nt --use-omp-stat-sched --thread-group-size $tgs --cache-size $cs | tee -a ${res_file} 

}

run_exp_no_prof() {


  t=$1
  top=$2
  aff=$3
  N=$4
  ts=$5
  k=$6
  nt=$7
  tgs=$8
  cs=$9

  res_file="${res_dir}kernel${k}_tgs${tgs}_ts${ts}_th${t}_thpgrp${top}_${aff}_affinity.txt"

  export OMP_NUM_THREADS=$t ; export KMP_PLACE_THREADS=$top ; export KMP_AFFINITY=${aff},verbose ;  /home/hpc/unrz/unrz280/girih/build_dp/mwd_kernel --n-tests 4 --nx $N --ny $N --nz $N --target-ts $ts --target-kernel $k --nt $nt --use-omp-stat-sched --thread-group-size $tgs --cache-size $cs | tee -a ${res_file} 

}


