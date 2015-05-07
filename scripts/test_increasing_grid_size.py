#!/usr/bin/env python
def igs_test(target_dir, exp_name, tgs_l): 
  from scripts.conf import conf
  from scripts.utils import run_test
  import itertools

  cs = conf.machine_info['cache_size']/4
  th = conf.machine_info['n_cores']
  kr = [4,1,1,1,4,1,1]
  k_time_scale = [1,1,1,1,2,2,20]
  #points = list(range(40, 961, 40)) + list(range(32, 1025, 32))
  points = list(range(32, 1025, 64)) + [1024]
  points = sorted(list(set(points)))
  kernels_limits = [1057, 1057, 0, 0, 545, 680, 289]
  print points
  count=0

  for kernel in [0, 1, 4, 5, 6]:
    for tgs in tgs_l:
      if tgs == 0:
        ts = 0
        mwdt_list=[0]
      else:
        ts = 2
        mwdt_list=[1,2,3]
      for mwdt in mwdt_list:
        for N in points:
          if (N < kernels_limits[kernel]) and (ts==0 or N > (th/tgs)*4*kr[kernel]) and not (kernel==6 and mwdt==3):
              outfile=('kernel%d_ts%d_mwdt%d_tgs%d_N%d.txt' % (kernel, ts, mwdt, tgs,  N))
              nt = max(int(th * 4e9/(N**3*3)/k_time_scale[kernel]), 30)
              run_test(dry_run=0, th=th,  kernel=kernel, ts=ts, nx=N, ny=N, nz=N, nt=nt,
                                 outfile=outfile, target_dir=target_dir, tgs=tgs, cs=cs, mwdt=mwdt)
              count = count+1
  print "experiments count =" + str(count)


def main():
  from scripts.utils import create_project_tarball 
  from scripts.conf.conf import machine_conf, machine_info

  #update the pinning information to use all cores
  th = machine_info['n_cores']

  if(machine_info['hostname']=='Haswell_18core'):
    machine_conf['pinning_args'] = "-c S0:0-%d "%(th-1) + machine_conf['pinning_args']
    tgs_l = [0, 1, 3, 6, 9, 18]
  elif(machine_info['hostname']=='IVB_10core'):
    machine_conf['pinning_args'] = "-C S0:0-%d -g MEM "%(th-1) + machine_conf['pinning_args']
    tgs_l = [0, 1, 2, 5, 10]

  exp_name = "increasing_grid_size_at_%s" % (machine_info['hostname'])  

  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)
  target_dir='results/' + exp_name 

  igs_test(target_dir, exp_name, tgs_l=tgs_l) 


if __name__ == "__main__":
  main()
