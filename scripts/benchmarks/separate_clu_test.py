#!/usr/bin/env python

def submit_experiment(kernel, ts, nx, ny, nz, nt, is_dp, outfile, target_dir, tgs, cs, exe_cmd):
  import os
  import subprocess
  from string import Template
  from utils import ensure_dir    

  job_template=Template(
"""$exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --nt $nt --target-kernel $kernel --cache-size $cs --thread-group-size $tgs | tee $outfile""")

  target_dir = os.path.join(os.path.abspath("."),target_dir)
  ensure_dir(target_dir)
  outpath = os.path.join(target_dir,outfile)

  if(is_dp==1):
    exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
  else:
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")

  job_cmd = job_template.substitute(cs=cs, tgs=tgs, nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel, ts=ts, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

  job_cmd = exe_cmd + job_cmd
 
  print job_cmd
  sts = subprocess.call(job_cmd, shell=True)



def naive(target_dir, exp_name, cs, exe_cmd): 
  is_dp = 1
  ts = 0
  tgs = 0
  points = list(range(32, 1025, 32))
  kernels_limits = [0, 1025, 0, 0, 521, 769]

  for kernel in [1, 4, 5]:
    for N in points:
      if N < kernels_limits[kernel]:
        outfile=(exp_name + 'isDP%d_ts%d_kernel%d_tgs%d_N%d.txt' % (is_dp, ts, kernel, tgs,  N))
        nx = N
        ny = N
        nz = N
        nt = max(int(10 * 4e9/(nx*ny*nz*3)), 50)
        submit_experiment(kernel, ts=ts, nx=nx, ny=ny, nz=nz, nt=nt, is_dp=is_dp, 
                               outfile=outfile, target_dir=target_dir, tgs=tgs, cs=cs, exe_cmd=exe_cmd)


def mwd(target_dir, exp_name, cs, exe_cmd, tgs, kernels): 
  is_dp = 1
  ts = 2
  kr = [4,1,1,1,4,1]
  points = list(range(32, 1025, 64))
  kernels_limits = [0, 1025, 0, 0, 521, 769]

  for kernel in kernels:
    for N in points:
      if N < kernels_limits[kernel]:
        if N > (10/tgs)*4*kr[kernel]:
          outfile=(exp_name + 'isDP%d_ts%d_kernel%d_tgs%d_N%d.txt' % (is_dp, ts, kernel, tgs,  N))
          nx = N
          ny = N
          nz = N
          nt = max(int(10 * 4e9/(nx*ny*nz*3)), 50)
          submit_experiment(kernel, ts=ts, nx=nx, ny=ny, nz=nz, nt=nt, is_dp=is_dp, 
                               outfile=outfile, target_dir=target_dir, tgs=tgs, cs=cs, exe_cmd=exe_cmd)


def main():

  from utils import create_project_tarball 

  import socket
  hostname = socket.gethostname()

  exp_name = "emmy_"
  run_cmd = "export OMP_NUM_THREADS=10; likwid-perfctr -m -C S0:0-9 -s 0x03 -g MEM "
  cs = 8192
 
 
  exp_name = exp_name + "test_no_separate_clu_at_%s" % (hostname)  

  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)

  target_dir='results/' + exp_name 

  mwd(target_dir, exp_name[:-9], 2048, run_cmd, 1, [1, 4, 5]) 
  mwd(target_dir, exp_name[:-9], 2048, run_cmd, 2, [1, 5]) 
  mwd(target_dir, exp_name[:-9], 2048, run_cmd, 5, [4]) 
  naive(target_dir, exp_name[:-9], cs, run_cmd)



if __name__ == "__main__":
  main()
