#!/usr/bin/env python

def submit_experiment(outfile, target_dir, ntest=4, kernel=1, ts=0, nx=480, ny=480, nz=480, npx=1, npy=1, npz=1, nt=50, is_dp=1, tgs=1, cs=8192, exe_cmd='', tb=-1, nwf=-1, mwd_type=0, bsx=100000):
  import os
  import subprocess
  from string import Template

  job_template=Template(
"""$exec_path --n-tests $ntest --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --nt $nt --target-kernel $kernel --cache-size $cs --thread-group-size $tgs --t-dim $tb --num-wavefronts $nwf --bsx $bsx --mwd-type $mwd_type | tee $outfile""")

  outpath = os.path.join(target_dir,outfile)

  if(is_dp==1):
    exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
  else:
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")

  job_cmd = job_template.substitute(cs=cs, tgs=tgs, nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel, ts=ts, outfile=outpath, exec_path=exec_path, target_dir=target_dir, tb=tb, nwf=nwf, mwd_type=mwd_type, ntest=ntest, npx=npx, npy=npy, npz=npz, bsx=bsx)

  job_cmd = exe_cmd + job_cmd
  return job_cmd


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
 

def mwd(target_dir, exp_name, cs, N, ts, kernel, tgs, th, mwd_type, groups, aff, tb=-1, nwf=-1, bsx=100000):
  is_dp = 1
  cmds =[]
  for group in groups:
    exe_cmd = "export OMP_NUM_THREADS=%d; likwid-perfctr -m -C %s -g %s -- " % (th, aff, group)
    outfile=(exp_name + '_ts%d_mwdtype%d_kernel%d_tgs%d_N%d_th%d_aff%s_group%s.txt' % (ts, mwd_type, kernel, tgs,  N, th, str.replace(aff,':','_'), group))

    nt =min(100,max(th,10))
    if tb==-1 or group!='CPI':
      cmd = submit_experiment(kernel=kernel, ts=ts, nx=N, ny=N, nz=N, nt=nt, is_dp=is_dp, 
                      outfile=outfile, target_dir=target_dir, tgs=tgs, cs=cs,
                      exe_cmd=exe_cmd, tb=tb, nwf=nwf, mwd_type=mwd_type, bsx=bsx)
      cmds.append([cmd])
  return cmds


def main():
  from utils import create_project_tarball 
  from csv import DictReader
  from utils import ensure_dir    
  import os
   

  exp_name = "mic_mwd_thread_scaling"
  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)
  target_dir='results/' + exp_name 
  target_dir = os.path.join(os.path.abspath("."),target_dir)
  ensure_dir(target_dir)


  # use auto tuning to find the best params at one group
  cs = 64
  ts = 2
  kernel=1; N=768

  # cache blocks within L2 caches
  cmds = []
  mwdt=2
  for tgs in [1,2,4]:
    th_l = [tgs*t  for t in range(6,61,6)]
    for th in th_l:
      aff="E:N:%d:%d:%d" % (th, tgs, 4)
      cmds = cmds + mwd(target_dir, exp_name, cs, N=N, ts=ts, kernel=kernel, tgs=tgs, th=th, mwd_type=mwdt, groups=['CPI'], aff=aff)

  cmds = [c for sublist in cmds for c in sublist]
  fname = target_dir + "/ts%d_kernel%d_mwdt%d_CPI_autotune.sh"%(ts, kernel,mwdt)
  with open(fname, 'w') as fn:
    fn.write("#!/bin/sh\n\n")
    for c in cmds:
      cp = 'echo ' + str.replace(c,';',"';'")
      cp = str.replace(cp, '|', '#|')
      fn.write('echo "'+c+'"\n')
      fn.write(c+'\n\n')




  # parse the results to obtain the selected parameters by the auto tuner
#  data = []
#  data_file = os.path.join(target_dir, 'summary.csv')
#  with open(data_file, 'rb') as output_file:
#    raw_data = DictReader(output_file)
#    for k in raw_data:
#      k['stencil'] = get_knum(k)
#      k['method'] = 2 if 'Diamond' in k['Time stepper orig name'] else 0
#      if k['method'] == 2:
#        if k['Wavefront parallel strategy'] == 'Relaxed synchronization wavefront with fixed execution':
#          k['mwdt'] = 3
#        elif k['Wavefront parallel strategy'] == 'Relaxed synchronization wavefront':
#          k['mwdt'] = 2
#        if int(k['Thread group size']) == 1:
#          k['mwdt'] = 0
#
#      data.append(k)
#
#  params = dict()
#  for k in data:
#    try:
#      if k['method']!=0:
#        params[(k['mwdt'], k['method'], k['stencil'], int(k['Thread group size']), int(k['Global NX']))]= (int(k['Time unroll']), int(k['Multi-wavefront updates']), int(k['Block size in X']))
#    except:
#      print k
#      raise
#
#
#  cs = 16000
#  th_l = [1,2,3,4,5,6,7,8,9,10]
#    for tgs, ths, mwdt_l in [(0,th_l, [0]), (1, th_l, [0]), (2,[2,4,6,8,10], [2,3]), (5, [5,10], [2,3]), (10, [10], [2,3])]:
#      for mwdt in mwdt_l:
#        for th in ths:
#
#          if tgs==0: # spatial blocking
#            ts = 0
#            tb, nwf, bsx = (-1,-1,100000)
#          else:
#            ts = 2
#            tb, nwf, bsx = params[(mwdt, ts, kernel, tgs, N)]
#
#          mwd(target_dir, exp_name, cs, N=N, ts=ts, kernel=kernel, tgs=tgs, th=th, mwd_type=mwdt, groups=['DATA', 'L2', 'L3', 'MEM'], tb=tb, nwf=nwf, bsx=bsx) 


if __name__ == "__main__":
  main()
