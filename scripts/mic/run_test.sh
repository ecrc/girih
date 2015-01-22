run_test() {


  t=$1
  aff=$2
  group=$3
  N=$4
  ts=$5
  k=$6
  nt=$7
  tgs=$8
  cs=$9
  mwdt=$10

  res_file="${res_dir}kernel${k}_tgs${tgs}_ts${ts}_th${t}_grp_${group}_affinity${aff}_mwdtype${mwdt}_N${N}.txt"

  thx=`expr $t - 1`

  export OMP_NUM_THREADS=$t ;likwid-perfctr -m -C $aff -g $group -- /home/hpc/unrz/unrz280/girih/build_dp/mwd_kernel --n-tests 4 --nx $N --ny $N --nz $N --target-ts $ts --target-kernel $k --nt $nt --use-omp-stat-sched --thread-group-size $tgs --cache-size $cs --mwd-type $mwdt | tee -a ${res_file} 

#  echo export OMP_NUM_THREADS=$t ';' export KMP_PLACE_THREADS=$top ';' export KMP_AFFINITY=${aff},verbose ';' likwid-perfctr -m -c 0-$thx -g $group -- /home/hpc/unrz/unrz280/girih/build_dp/mwd_kernel --n-tests 4 --nx $N --ny $N --nz $N --target-ts $ts --target-kernel $k --nt $nt --use-omp-stat-sched --thread-group-size $tgs --cache-size $cs | tee -a ${res_file} 

}

