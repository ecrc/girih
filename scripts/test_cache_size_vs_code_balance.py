#!/usr/bin/env python
def igs_test(target_dir, exp_name, th, group='', dry_run=0): 
  from scripts.conf.conf import machine_conf, machine_info
  from scripts.utils import run_test
  import itertools

  cs = 8192

  th = th

  # Test using rasonable time
  # T = scale * size / perf
  # scale = T*perf/size
  desired_time = 20
  if(machine_info['hostname']=='Haswell_18core'):
    k_perf_order = {0:150, 1:500, 4:40, 5:200 ,6:20}
  elif(machine_info['hostname']=='IVB_10core'):
    k_perf_order = {0:120, 1:300, 4:35, 5:150 ,6:20}


  k_time_scale = {n: desired_time*k_perf_order[n] for n in k_perf_order.keys()}

  #exp = is_dp, ts, k, N, bs_z,  tb_l
  exp_l = []
  # spatial blocking
  exp_l = exp_l + \
          [(0, 0, 0, 960, 0, [-1])
          ,(1, 0, 0, 960, 0, [-1])
          ,(1, 0, 1, 960, 0, [-1])
          ,(1, 0, 4, 480, 0, [-1])
          ,(1, 0, 5, 680, 0, [-1])
          ]
  # 1WD
  exp_l = exp_l + \
          [(0, 2, 0, 960, 1, [1, 3, 5])
          ,(1, 2, 0, 960, 1, [1, 3, 5])
          ,(1, 2, 1, 960, 1, [1, 3, 5, 7, 9, 11, 15, 19, 23, 29])
          ,(1, 2, 4, 480, 1, [1, 3, 5])
          ,(1, 2, 5, 680, 1, [1, 3, 9, 19])
          ]

# Solar kernel
  exp_l = exp_l + \
          [(1, 2, 6, 480, 1, [1, 3, 5, 7])
          ,(1, 2, 6, 480, 2, [1, 3, 5, 7])
          ,(1, 2, 6, 480, 3, [1, 3, 5, 7])
          ,(1, 2, 6, 480, 6, [1, 3, 5, 7])
          ,(1, 2, 6, 480, 9, [1, 3, 5, 7])]

  mwdt=1
  tgs, thx, thy, thz = (1,1,1,1)
  count=0
  for is_dp, ts, kernel, N, bs_z, tb_l in exp_l:
    for tb in tb_l:
      outfile=('kernel%d_isdp%d_ts%d_bsz$d_tb%d_N%d_%s_%s.txt' % (kernel, is_dp, ts, bs_z, tb, N, group, exp_name[-13:]))
      nt = max(int(k_time_scale[kernel]/(N**3/1e6)), 30)
#      print outfile, ts, kernel, tb, N  
      run_test(ntests=1,dry_run=dry_run, is_dp=is_dp, th=th, tgs=tgs, thx=thx, thy=thy, thz=thz, kernel=kernel, ts=ts, nx=N, ny=N, nz=N, nt=nt, outfile=outfile, target_dir=target_dir, cs=cs, mwdt=mwdt, tb=tb, nwf=bs_z)
      count = count+1

  return count


def main():
  from scripts.utils import create_project_tarball, get_stencil_num, parse_results
  from scripts.conf.conf import machine_conf, machine_info
  import os, sys
  import time,datetime

  # user params
  dry_run = 1   if len(sys.argv)<2 else int(sys.argv[1]) # dry run

  time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H_%M')
  exp_name = "cache_size_vs_code_balance_at_%s_%s" % (machine_info['hostname'], time_stamp)  

  tarball_dir='results/'+exp_name
  if(dry_run==0): create_project_tarball(tarball_dir, "project_"+exp_name)
  target_dir='results/' + exp_name 

  th = 1
  pin_str = "S0:0-%d "%(th-1)

  count=0
  group = 'MEM'
  if( (machine_info['hostname']=='IVB_10core') and (group=='TLB_DATA') ): group='TLB'
  machine_conf['pinning_args'] = "-m -g " + group + " -C " + pin_str + ' -s 0x03 --'

  count= count + igs_test(target_dir, exp_name, th=th, group=group, dry_run=dry_run) 

  print "experiments count =" + str(count)


if __name__ == "__main__":
  main()
