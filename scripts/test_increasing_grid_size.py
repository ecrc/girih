#!/usr/bin/env python
def igs_test(target_dir, exp_name, tgs_l, th, params={}): 
  from scripts.conf.conf import machine_conf, machine_info
  from scripts.utils import run_test
  import itertools

  dry_run = 0
  is_dp=0

  cs = machine_info['cache_size']/4
  kr = [4,1,1,1,4,1,1]
  k_time_scale = [1,1,1,1,2,2,20]

  points = list(range(32, 5000, 128))
  points = sorted(list(set(points)))
  if is_dp ==1:
    kernels_limits = [1057, 1057, 0, 0, 545, 680, 289]
  else:
    kernels_limits = [1350, 0, 0, 0, 801, 0, 0]

  if(machine_info['hostname']=='Haswell_18core'):
    if is_dp == 1:
      kernels_limits = [1600, 1600, 0, 0, 960, 1000, 500]
    else:
      kernels_limits = [2100, 0, 0, 0, 1200, 0, 0]

  count=0
  for kernel in [0, 1, 4, 5, 6]:
    for tgs in tgs_l:
      ts = 2
      if tgs == 0:
        ts = 0
        mwdt_list=[0]
      elif tgs == 1:
        mwdt_list=[0]
      else:
        mwdt_list=[0,1,2,3]
      for mwdt in mwdt_list:
        for N in points:
          if (N < kernels_limits[kernel]) and (ts==0 or N >= (th/tgs)*4*kr[kernel]) and not (kernel==6 and mwdt==3):
              tb, nwf, bsx = (-1,-1,100000)
              key = (mwdt, kernel, tgs, N)
              if key in params.keys():
                tb, nwf, bsx = params[key]
#              if tgs==0 or tb !=-1: continue

              outfile=('kernel%d_isdp%d_ts%d_mwdt%d_tgs%d_N%d.txt' % (kernel, is_dp, ts, mwdt, tgs,  N))
              nt = max(int(th * 4e9/(N**3*3)/k_time_scale[kernel]), 30)
              run_test(dry_run=dry_run, is_dp=is_dp, th=th,  kernel=kernel, ts=ts, nx=N, ny=N, nz=N, nt=nt,
                                 outfile=outfile, target_dir=target_dir, tgs=tgs, cs=cs, mwdt=mwdt, tb=tb, nwf=nwf)
              count = count+1
  print "experiments count =" + str(count)


def main():
  from scripts.utils import create_project_tarball, get_stencil_num
  from scripts.conf.conf import machine_conf, machine_info
  import os
  from csv import DictReader

  sockets=1 # number of processors to use in the experiments

  #update the pinning information to use all cores
  th = machine_info['n_cores']*sockets

  if sockets == 1:
    pin_str = "S0:0-%d "%(th-1)
  if sockets == 2:
    pin_str = "S0:0-%d@S1:0-%d -i "%(th/2-1, th/2-1)



  if(machine_info['hostname']=='Haswell_18core'):
    machine_conf['pinning_args'] = "-m -g MEM -C " + pin_str + machine_conf['pinning_args']
    tgs_l = [0, 1, 3, 6, 9, 18]
  elif(machine_info['hostname']=='IVB_10core'):
    machine_conf['pinning_args'] = "-m -g MEM -C " + pin_str + machine_conf['pinning_args']
#    machine_conf['pinning_args'] = "-c " + pin_str + machine_conf['pinning_args']
    tgs_l = [0, 1, 2, 5, 10]

  exp_name = "increasing_grid_size_SP_sockets_%d_at_%s" % (sockets,machine_info['hostname'])  

  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)
  target_dir='results/' + exp_name 


  # parse the results to obtain the selected parameters by the auto tuner
  data = []
  data_file = os.path.join(target_dir, 'summary.csv')
  try:
    with open(data_file, 'rb') as output_file:
      raw_data = DictReader(output_file)
      for k in raw_data:
        k['stencil'] = get_stencil_num(k)
        k['method'] = 2 if 'Diamond' in k['Time stepper orig name'] else 0
        if k['method'] == 2:

          if k['Wavefront parallel strategy'] == 'Relaxed synchronization wavefront with fixed execution':
            k['mwdt'] = 3
          elif k['Wavefront parallel strategy'] == 'Relaxed synchronization wavefront':
            k['mwdt'] = 2
          elif k['Wavefront parallel strategy'] == 'Wavefront':
            k['mwdt'] = 0
          elif k['Wavefront parallel strategy'] == 'Fixed execution wavefronts':
            k['mwdt'] = 1
          if int(k['Thread group size']) == 1:
            k['mwdt'] = 0
        data.append(k)
  except:
     pass
  params = dict()
  for k in data:
    try:
      if k['method']==2:
        params[(k['mwdt'], k['stencil'], int(k['Thread group size']), int(k['Global NX']))]= (int(k['Time unroll']), int(k['Multi-wavefront updates']), int(k['Block size in X']))
    except:
      print k
      raise


  igs_test(target_dir, exp_name, tgs_l=tgs_l, th=th, params=params) 


if __name__ == "__main__":
  main()
