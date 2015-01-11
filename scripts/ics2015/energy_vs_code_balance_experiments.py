#!/usr/bin/env python

def submit_experiment(kernel, Dw, nwf, ts, nx, ny, nz, nt, is_dp, outfile, target_dir, tgs, cs, exe_cmd):
  import os
  import subprocess
  from string import Template
  from utils import ensure_dir    

  job_template=Template(
"""$exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --nt $nt --target-kernel $kernel --cache-size $cs --thread-group-size $tgs --t-dim $Dw --num-wavefronts $nwf --mwd-type 2 | tee $outfile""")

  target_dir = os.path.join(os.path.abspath("."),target_dir)
  ensure_dir(target_dir)
  outpath = os.path.join(target_dir,outfile)

  if(is_dp==1):
    exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
  else:
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")

  job_cmd = job_template.substitute(Dw=Dw, nwf=nwf, cs=cs, tgs=tgs, nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel, ts=ts, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

  job_cmd = exe_cmd + job_cmd
 
  print job_cmd
  sts = subprocess.call(job_cmd, shell=True)



def ts_test(target_dir, exp_name, cs, ts, kernel, tgs, N, group, dw_l):
  is_dp = 1

  if kernel==5: nwf = 20   
  if kernel==1: nwf = 35 
  for Dw in dw_l:
    exe_cmd = "export OMP_NUM_THREADS=10; likwid-perfctr -m -C S0:0-9 -s 0x03 -g %s " % (group)
    outfile=('ts%d_kernel%d_tgs%d_N%d_%s_TB%d.txt' % (ts, kernel, tgs,  N, group,Dw))
    nx = N
    ny = N
    nz = N
    nt =40 # int(max(10 * 4e9/(nx*ny*nz*3), 50))
    submit_experiment(kernel, Dw=Dw, nwf=nwf, ts=ts, nx=nx, ny=ny, nz=nz, nt=nt, is_dp=is_dp, 
         outfile=outfile, target_dir=target_dir, tgs=tgs, cs=cs, exe_cmd=exe_cmd)


def main():

  from utils import create_project_tarball 

  import socket
  #hostname = socket.gethostname()
  exp_name = "emmy_"
  exp_name = exp_name + "energy_vs_code_balance"#_at_%s" % (hostname)  

  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)

  target_dir='results/' + exp_name 

  cs = 4096
  dw_l = [1,3,5,7]
  ts_test(target_dir, exp_name, cs, ts=2, kernel=1, tgs=5, N=960, group='ENERGY', dw_l=dw_l) 
  ts_test(target_dir, exp_name, cs, ts=2, kernel=1, tgs=5, N=960, group='MEM', dw_l=dw_l) 

  ts_test(target_dir, exp_name, cs, ts=2, kernel=5, tgs=5, N=480, group='ENERGY', dw_l=dw_l) 
  ts_test(target_dir, exp_name, cs, ts=2, kernel=5, tgs=5, N=480, group='MEM', dw_l=dw_l) 

if __name__ == "__main__":
  main()
