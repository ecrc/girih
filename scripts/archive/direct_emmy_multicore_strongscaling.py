#!/usr/bin/env python

def submit_emmy_experiment(kernel, ts, tb, perf_ctr, nx, ny, nz, nt, is_dp, outfile, target_dir, tgs, fuse_wf, cs):
  import os
  import subprocess
  from string import Template
  from utils import ensure_dir    

  job_template=Template(
"""export OMP_NUM_THREADS=10; likwid-perfctr -c S0:0-9 -s 0x03 -g $perf_ctr $exec_path --n-tests 3 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --nt $nt --target-kernel $kernel --t-dim $tb --cache-size $cs --thread-group-size $tgs --fuse-wavefronts $fuse_wf | tee $outfile""")


  target_dir = os.path.join(os.path.abspath("."),target_dir)
  ensure_dir(target_dir)
  outpath = os.path.join(target_dir,outfile)

  if(is_dp==1):
    exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
  else:
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


  job_cmd = job_template.substitute(cs=cs, tgs=tgs, fuse_wf=fuse_wf, tb=tb, perf_ctr=perf_ctr, nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel, ts=ts, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

  #sts = subprocess.call(job_cmd, shell=True)
  print job_cmd


def naive_nx512(perf_ctr); 
  is_dp = 1
  tb = 0
  ts = 0
  nx = 512
  kernel = 2
  cs = 11000
  tgs = 0
  fuse_wf = 0
  for N in range(32,528, 16):
    outfile=(exp_name + 'isDP%d_ts%d_kernel%d_cs%d_tgs%d_fuse_wf%d_N%d_TB%d' % (is_dp, ts, kernel, cs, tgs, fuse_wf,  N, tb))
    ny = N
    nz = N
    nt = 10 * 4e9/(nx*ny*nz*3)
    submit_emmy_experiment(kernel, tb=tb, ts=ts, perf_ctr=perf_ctr, nx=nx, ny=ny, nz=nz, nt=nt, is_dp=is_dp, 
                           outfile=outfile, target_dir=target_dir, tgs, fuse_wf=fuse_wf, cs=cs)


def main():

  from utils import create_project_tarball 

  import socket
  hostname = socket.gethostname()

  exp_name = "emmy_strongscaling_nx512_at_%s" % (hostname)  
  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)

  target_dir='results/' + exp_name 

  perf_ctr = 'MEM'

  naive_nx512(perf_ctr)


if __name__ == "__main__":
  main()
