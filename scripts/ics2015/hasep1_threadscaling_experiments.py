#!/usr/bin/env python

def submit_experiment(kernel, ts, nx, ny, nz, nt, is_dp, outfile, target_dir, tgs, cs, exe_cmd):
  import os
  import subprocess
  from string import Template
  from utils import ensure_dir    

  job_template=Template(
"""$exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --nt $nt --target-kernel $kernel --cache-size $cs --thread-group-size $tgs --mwd-type 2 | tee $outfile""")

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



def ts_test(target_dir, exp_name, cs, ts, kernel, tgs, N, th_l): 
  is_dp = 1
  for th in th_l:
    if th%tgs == 0:
      exe_cmd = "export OMP_NUM_THREADS=%d; likwid-perfctr -m -C S0:0-%d -s 0x01 -g DATA " % (th, th-1)
      outfile=('ts%d_kernel%d_tgs%d_N%d_th%d.txt' % (ts, kernel, tgs,  N, th))
      nx = N
      ny = N
      nz = N

      min_nt = 50 if ts==2 else 30
      nt = int(max((14.0 * 4e9 / th) /(nx*ny*nz*3), min_nt))
      submit_experiment(kernel, ts=ts, nx=nx, ny=ny, nz=nz, nt=nt, is_dp=is_dp, 
         outfile=outfile, target_dir=target_dir, tgs=tgs, cs=cs, exe_cmd=exe_cmd)


def main():

  from utils import create_project_tarball 

  import socket
#  hostname = socket.gethostname()
  exp_name = "hasep1_"
  exp_name = exp_name + "threadscaling"#_at_%s" % (hostname)  

  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)

  target_dir='results/' + exp_name 

  th_l = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
  cs = 4096
  for tgs in [1,2,7,14]:
#    ts_test(target_dir, exp_name, cs, ts=2, kernel=1, tgs=tgs, N=960, th_l=th_l) 
#    ts_test(target_dir, exp_name, cs, ts=2, kernel=4, tgs=tgs, N=480, th_l=th_l) 
    ts_test(target_dir, exp_name, cs, ts=2, kernel=5, tgs=tgs, N=680, th_l=th_l) 


  cs = 16384
#  ts_test(target_dir, exp_name[:-8], cs, ts=0, kernel=1, tgs=0, N=960, th_l=[8]) 
#  ts_test(target_dir, exp_name[:-8], cs, ts=0, kernel=4, tgs=0, N=480, th_l=[3, 7, 10]) 
  ts_test(target_dir, exp_name, cs, ts=0, kernel=5, tgs=1, N=680, th_l=th_l) 



if __name__ == "__main__":
  main()
