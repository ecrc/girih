#!/usr/bin/env python

def submit_experiment(outfile, target_dir, dry_run=0, ntest=3, kernel=1, ts=0, nx=480, ny=480, nz=480, npx=1, npy=1, npz=1, nt=50, is_dp=1, tgs=1, cs=8192, exe_cmd='', tb=-1, nwf=-1, mwd_type=0, bsx=100000):
  import os
  import subprocess
  from string import Template
  from utils import ensure_dir    

  job_template=Template(
"""$exec_path --n-tests $ntest --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --nt $nt --target-kernel $kernel --cache-size $cs --thread-group-size $tgs --t-dim $tb --num-wavefronts $nwf --bsx $bsx --mwd-type $mwd_type | tee $outfile""")

  target_dir = os.path.join(os.path.abspath("."),target_dir)
  ensure_dir(target_dir)
  outpath = os.path.join(target_dir,outfile)

  if(is_dp==1):
    exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
  else:
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")

  job_cmd = job_template.substitute(cs=cs, tgs=tgs, nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel, ts=ts, outfile=outpath, exec_path=exec_path, target_dir=target_dir, tb=tb, nwf=nwf, mwd_type=mwd_type, ntest=ntest, npx=npx, npy=npy, npz=npz, bsx=bsx)

  job_cmd = exe_cmd + job_cmd
 
  print job_cmd
  if dry_run==0: sts = subprocess.call(job_cmd, shell=True)


def get_knum(k):
  # add the stencil operator
  if  k['Stencil Kernel coefficients'] in 'constant':
    if  int(k['Stencil Kernel semi-bandwidth'])==4:
      stencil = 0
    else:   
      stencil = 1
  elif  'no-symmetry' in k['Stencil Kernel coefficients']:
    stencil = 5
  elif  'sym' in k['Stencil Kernel coefficients']:
    if int(k['Stencil Kernel semi-bandwidth'])==1:
      stencil = 3
    else:   
      stencil = 4
  else:
    stencil = 2
  return stencil
 

def mwd(target_dir, exp_name, cs, N, ts, kernel, tgs, th, mwd_type, groups, tb=-1, nwf=-1, bsx=100000):
  is_dp = 1
  for group in groups:
    exe_cmd = "export OMP_NUM_THREADS=%d; likwid-perfctr -m -C S0:0-%d -s 0x01 -g %s -- " % (th, (th-1), group)
    outfile=(exp_name + 'ts%d_mwdtype%d_kernel%d_tgs%d_N%d_th%d_group%s.txt' % (ts, mwd_type, kernel, tgs,  N, th, group))
    nt = max(int(1.0*th*4e9 / (N**3*3)), 30)
    if kernel==5: nt=nt/3

    if tb==-1 or group!='DATA' or th!=14:
      submit_experiment(dry_run=0, kernel=kernel, ts=ts, nx=N, ny=N, nz=N, nt=nt, is_dp=is_dp, 
                      outfile=outfile, target_dir=target_dir, tgs=tgs, cs=cs,
                      exe_cmd=exe_cmd, tb=tb, nwf=nwf, mwd_type=mwd_type, bsx=bsx)



def main():
  from utils import create_project_tarball 
  import socket
  from csv import DictReader
  import os
   

  hostname = socket.gethostname()
  exp_name = "hasep1_"
  exp_name = exp_name + "prof_mwd_cache_"
  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)
  target_dir='results/' + exp_name 

  grid_set = {1:[480, 960], 5:[256, 640]}


# # use auto tuning to find the best params
#  th=14
#  cs = 6000 
#  ts = 2
#  for kernel in [1, 5]:
#    for N in grid_set[kernel]:
#      for tgs, mwt_l in [(1,[2]), (2, [2,3]), (7, [2,3]), (14, [2,3])]:
#        for mwdt in [2,3]:
#          mwd(target_dir, exp_name, cs, N=N, ts=ts, kernel=kernel, tgs=tgs, th=th, mwd_type=mwdt, groups=['DATA']) 
#

 
  # parse the results to obtain the selected parameters by the auto tuner
  data = []
  data_file = os.path.join(target_dir, 'summary.csv')
  with open(data_file, 'rb') as output_file:
    raw_data = DictReader(output_file)
    for k in raw_data:
      k['stencil'] = get_knum(k)
      k['method'] = 2 if 'Diamond' in k['Time stepper orig name'] else 0
      if k['method'] == 2:
        if k['Wavefront parallel strategy'] == 'Relaxed synchronization wavefront with fixed execution':
          k['mwdt'] = 3
        elif k['Wavefront parallel strategy'] == 'Relaxed synchronization wavefront':
          k['mwdt'] = 2
        if int(k['Thread group size']) == 1:
          k['mwdt'] = 0

      data.append(k)

  params = dict()
  for k in data:
    try:
      if k['method']!=0:
        params[(k['mwdt'], k['method'], k['stencil'], int(k['Thread group size']), int(k['Global NX']))]= (int(k['Time unroll']), int(k['Multi-wavefront updates']), int(k['Block size in X']))
    except:
      print k
      raise


  cs = 16000
  th_l = [1,2,4,6,8,10,12,14]
  for kernel in [1, 5]:
    for N in grid_set[kernel]:
      for tgs, ths, mwdt_l in [(0,th_l, [0]), (1, th_l, [0]), (2,th_l[1:], [2,3]), (7, [7,14], [2,3]), (14, [14], [2,3])]:
        for mwdt in mwdt_l:
          for th in ths:

            if tgs==0: # spatial blocking
              ts = 0
              tb, nwf, bsx = (-1,-1,100000)
            else:
              ts = 2
              tb, nwf, bsx = params[(mwdt, ts, kernel, tgs, N)]

#            mwd(target_dir, exp_name, cs, N=N, ts=ts, kernel=kernel, tgs=tgs, th=th, mwd_type=mwdt, groups=['DATA', 'L2', 'L3', 'MEM'], tb=tb, nwf=nwf, bsx=bsx) 
            mwd(target_dir, exp_name, cs, N=N, ts=ts, kernel=kernel, tgs=tgs, th=th, mwd_type=mwdt, groups=['DATA', 'L2', 'L3'], tb=tb, nwf=nwf, bsx=bsx) 


if __name__ == "__main__":
  main()
