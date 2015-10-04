#!/usr/bin/env python
def thread_scaling_test(target_dir, exp_name, group='', params={}, dry_run=0, is_tgs_only=0): 
  from scripts.conf.conf import machine_conf, machine_info
  from scripts.utils import run_test
  import itertools

  is_dp=1

  cs = 4096


  # Test using rasonable time
  # T = scale * size / perf
  # scale = T*perf/size
  desired_time = 20
  if(machine_info['hostname']=='Haswell_18core'):
    k_perf_order = {0:1500, 1:5000, 4:400, 5:2000 ,6:100}
  elif(machine_info['hostname']=='IVB_10core'):
    k_perf_order = {0:1200, 1:3000, 4:350, 5:1500 ,6:80}
  k_time_scale = {n: desired_time*k_perf_order[n] for n in k_perf_order.keys()}


  points = dict()
  points[0] = [896]
  points[1] = [896]
  points[4] = [768]
  points[5] = [768]
  points[6] = [384]

  count=0
  th_l = [machine_info['n_cores']]
  th_max = machine_info['n_cores']

  for ts, tgs_r in [(0,0), (2,1), (2,-1)]:
    for kernel, mwdt in [(0,1), (1,2), (4,1), (5,2)]: #, (6,1)]:
      for N in points[kernel]:
        tb, nwf, tgs, thx, thy, thz = (-1,-1,tgs_r,-1,-1,-1)
#        key = (mwdt, kernel, N, tgs_r, group)
#        if key in params.keys():
#          continue #already computed
        if tgs==1: mwdt=-1
        key = (mwdt, kernel, N, tgs_r, 'MEM')
        if key in params.keys(): #reuse the conf of existing test
          tb, nwf, tgs, thx, thy, thz = params[key]
        elif ts==0:
          pass
        else:
          print "full socket test required for ts:%d tgs:%d kernel:%d N:%d "%(ts, tgs_r, kernel, N)
          continue

        for th in list(range(1, 1+th_max)):
          if (tgs>1  and th%tgs!=0):
            continue
          machine_conf['pinning_args'] = "-m -g %s -C S1:0-%d -s 0x03 --"%(group, (th-1))
          outfile=('kernel%d_isdp%d_ts%d_mwdt%d_tgs%d_N%d_th%d_%s_%s.txt' % (kernel, is_dp, ts, mwdt, tgs_r, N, th, group, exp_name[-13:]))
          nt = max(int( float(th)/float(th_max)* k_time_scale[kernel]/(N**3/1e6)), 30)

          c_mwdt = mwdt if mwdt!=-1 else 1

          print outfile, th, tb, nwf, tgs, thx, thy, thz
          run_test(dry_run=dry_run, is_dp=is_dp, th=th, tgs=tgs, thx=thx, thy=thy, thz=thz, kernel=kernel, ts=ts, nx=N, ny=N, nz=N, nt=nt, outfile=outfile, target_dir=target_dir, cs=cs, mwdt=c_mwdt, tb=tb, nwf=nwf)
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
  exp_name = "thread_scaling_at_%s_%s" % (machine_info['hostname'], time_stamp)  

  tarball_dir='results/'+exp_name
  if(dry_run==0): create_project_tarball(tarball_dir, "project_"+exp_name)
  target_dir='results/' + exp_name 

  # parse the results to obtain the selected parameters by the auto tuner
  params = parse_results()

  count=0
  for group in ['MEM']:
    if( (machine_info['hostname']=='IVB_10core') and (group=='TLB_DATA') ): group='TLB'

#    for k,v in params.iteritems():
#      if k[2]==896: print k,v

    count= count + thread_scaling_test(target_dir, exp_name, params=params, group=group, dry_run=dry_run) 

  print "experiments count =" + str(count)


if __name__ == "__main__":
  main()
