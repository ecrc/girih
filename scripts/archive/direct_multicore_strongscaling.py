#!/usr/bin/env python

def submit_experiment(kernel, ts, tb, nx, ny, nz, nt, is_dp, outfile, target_dir, tgs, fuse_wf, cs, exe_cmd):
  import os
  import subprocess
  from string import Template
  from utils import ensure_dir    

  job_template=Template(
"""$exec_path --n-tests 3 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --nt $nt --target-kernel $kernel --t-dim $tb --cache-size $cs --thread-group-size $tgs --fuse-wavefronts $fuse_wf | tee $outfile""")


  target_dir = os.path.join(os.path.abspath("."),target_dir)
  ensure_dir(target_dir)
  outpath = os.path.join(target_dir,outfile)

  if(is_dp==1):
    exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
  else:
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


  job_cmd = job_template.substitute(cs=cs, tgs=tgs, fuse_wf=fuse_wf, tb=tb, nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel, ts=ts, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

  job_cmd = exe_cmd + job_cmd
 
  print job_cmd
  sts = subprocess.call(job_cmd, shell=True)


def naive_nx512(target_dir, exp_name, cs, exe_cmd): 
  is_dp = 1
  tb = 0
  ts = 0
  nx = 512
  kernel = 2
  tgs = 0
  fuse_wf = 0
  for N in range(32,528, 16):
    outfile=(exp_name + 'isDP%d_ts%d_kernel%d_cs%d_tgs%d_fuse_wf%d_N%d_TB%d' % (is_dp, ts, kernel, cs, tgs, fuse_wf,  N, tb))
    ny = N
    nz = N
    nt = int(10 * 4e9/(nx*ny*nz*3))
    submit_experiment(kernel, tb=tb, ts=ts, nx=nx, ny=ny, nz=nz, nt=nt, is_dp=is_dp, 
                           outfile=outfile, target_dir=target_dir, tgs=tgs, fuse_wf=fuse_wf, cs=cs, exe_cmd=exe_cmd)

def swd_nx512(target_dir, exp_name, cs, exe_cmd): 
  is_dp = 1
  tb = -1
  ts = 9
  nx = 512
  kernel = 2
  tgs = 1
  fuse_wf = 0
  for N in range(48,528, 16):
    outfile=(exp_name + 'isDP%d_ts%d_kernel%d_cs%d_tgs%d_fuse_wf%d_N%d_TB%d' % (is_dp, ts, kernel, cs, tgs, fuse_wf,  N, tb))
    ny = N
    nz = N
    nt = int(10 * 4e9/(nx*ny*nz*3))
    submit_experiment(kernel, tb=tb, ts=ts, nx=nx, ny=ny, nz=nz, nt=nt, is_dp=is_dp, 
                           outfile=outfile, target_dir=target_dir, tgs=tgs, fuse_wf=fuse_wf, cs=cs, exe_cmd=exe_cmd)

def m2wd_nx512(target_dir, exp_name, cs, exe_cmd): 
  is_dp = 1
  tb = -1
  ts = 9
  nx = 512
  kernel = 2
  tgs = 2
  fuse_wf = 0
  for N in range(32,528, 16):
    outfile=(exp_name + 'isDP%d_ts%d_kernel%d_cs%d_tgs%d_fuse_wf%d_N%d_TB%d' % (is_dp, ts, kernel, cs, tgs, fuse_wf,  N, tb))
    ny = N
    nz = N
    nt = int(10 * 4e9/(nx*ny*

nz*3))
    submit_experiment(kernel, tb=tb, ts=ts, nx=nx, ny=ny, nz=nz, nt=nt, is_dp=is_dp, 
                           outfile=outfile, target_dir=target_dir, tgs=tgs, fuse_wf=fuse_wf, cs=cs, exe_cmd=exe_cmd)


def mwd_nx512(target_dir, exp_name, cs, exe_cmd, tgs): 
  is_dp = 1
  tb = -1
  ts = 9
  nx = 512
  kernel = 2
  fuse_wf = 0
  for N in range(32,528, 16):
    outfile=(exp_name + 'isDP%d_ts%d_kernel%d_cs%d_tgs%d_fuse_wf%d_N%d_TB%d' % (is_dp, ts, kernel, cs, tgs, fuse_wf,  N, tb))
    ny = N
    nz = N
    nt = int(10 * 4e9/(nx*ny*

nz*3))
    submit_experiment(kernel, tb=tb, ts=ts, nx=nx, ny=ny, nz=nz, nt=nt, is_dp=is_dp, 
                           outfile=outfile, target_dir=target_dir, tgs=tgs, fuse_wf=fuse_wf, cs=cs, exe_cmd=exe_cmd)





def main():

  from utils import create_project_tarball 

  import socket
  hostname = socket.gethostname()

  #exp_name = "emmy_"
  #run_cmd = "export OMP_NUM_THREADS=10; likwid-perfctr -c S0:0-9 -s 0x03 -g MEM "
  #cs = 11000
 
  exp_name = "stampede_"
  run_cmd = "export OMP_NUM_THREADS=8; numactl -N 0 "
  cs = 8192
  
  exp_name = exp_name + "strongscaling_mwd_naive_nx512_at_%s" % (hostname)  

  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)

  target_dir='results/' + exp_name 



 # naive_nx512(target_dir, exp_name, cs, run_cmd)
 # swd_nx512(target_dir, exp_name, cs, run_cmd)
 # m2wd_nx512(target_dir, exp_name, cs, run_cmd)
  mwd_nx512(target_dir, exp_name, cs, run_cmd, 4)
  mwd_nx512(target_dir, exp_name, cs, run_cmd, 8)


if __name__ == "__main__":
  main()
