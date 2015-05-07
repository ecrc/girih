#!/usr/bin/env python

def igs_test(target_dir, exp_name, tgs_l): 
  from scripts import conf
  from scripts.utils import run_test
  import itertools

  cs = conf.machine_info['cache_size']/4
  th = conf.machine_info['n_cores']
  kr = [4,1,1,1,4,1,1]
  #points = list(range(40, 961, 40)) + list(range(32, 1025, 32))
  points = list(range(32, 1025, 64)) + [1024]
  points = sorted(list(set(points)))
  kernels_limits = [1025, 1025, 0, 0, 521, 681, 240]
  print points
  count=0

  for kernel in [0, 1, 4, 5, 6]:
    for (ts, tgs_list, mwdt_list) in [(0, [1], [0]), (2, tgs_l, [1,2,3])]:
      for mwdt in mwdt_list:
        for tgs in tgs_list:
          for N in points:
            if (N < kernels_limits[kernel]) and (ts==0 or N > (th/tgs)*4*kr[kernel]) and not (kernel==6 and mwdt==3):
                outfile=('kernel%d_ts%d_mwdt%d_tgs%d_N%d.txt' % (kernel, ts, mwdt, tgs,  N))
                nt = max(int(th * 4e9/(N**3*3)), 30)
                run_test(dry_run=0, th=th,  kernel=kernel, ts=ts, nx=N, ny=N, nz=N, nt=nt,
                                   outfile=outfile, target_dir=target_dir, tgs=tgs, cs=cs, mwdt=mwdt)
                count = count+1
  print "experiments count =" + str(count)

def main():
  from scripts.utils import create_project_tarball 
  from scripts.conf import machine_conf, machine_info

  #update the pinning information to use all cores
  th = machine_info['n_cores']
  machine_conf['pinning_args'] = "-c S0:0-%d "%(th-1) + machine_conf['pinning_args']
 
  exp_name = "increasing_grid_size_at_%s" % (machine_info['hostname'])  

  tarball_dir='results/'+exp_name
  create_project_tarball(tarball_dir, "project_"+exp_name)
  target_dir='results/' + exp_name 

  igs_test(target_dir, exp_name, tgs_l=[3]) 


if __name__ == "__main__":
  main()
