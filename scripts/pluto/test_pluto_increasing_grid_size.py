#!/usr/bin/env python
def tee(fp, string):
  print(string)
  fp.write(string)

def run_pluto_test(dry_run, kernel, nx, ny, nz, nt, params, outfile='', pinning_cmd=-1, pinning_args=-1, auto_tuning=0):
  import os, subprocess
  from os.path import join as jpath
  from string import Template
  from scripts.conf.conf import machine_conf, machine_info
  import time

  if(pinning_cmd==-1): pinning_cmd = machine_conf['pinning_cmd']
  if(pinning_args==-1): pinning_args = machine_conf['pinning_args']

  job_template=Template(
"""$pinning_cmd $pinning_args $exec_path $nx $ny $nz $nt $outfile""")

  if(outfile!=''):
    outfile = ' | tee -a ' + outfile

  # set the executable
  exec_name = 'lbpar_' + kernel
  #add the tile size parameters to the executable name
  exec_dir = exec_name + "%d_%d_%d_%d"%(params[0], params[0], params[1], params[2])

  exec_path = jpath(os.path.abspath("."),'pluto_examples', 'gen_kernels', exec_dir,  exec_name)
 
   
  job_cmd = job_template.substitute(nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel,
                       outfile=outfile, exec_path=exec_path, pinning_cmd=pinning_cmd, 
                       pinning_args=pinning_args)
 
  tstart = time.time()
  test_str=''
  if(auto_tuning):
    proc = subprocess.Popen(job_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    test_str = test_str + proc.stdout.read()
  else:
    print job_cmd
    test_str = job_cmd + '\n'
    if(dry_run==0): 
      sts = subprocess.call(job_cmd, shell=True)
      proc = subprocess.Popen(job_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      test_str = test_str + proc.stdout.read()
  tend = time.time()
  return test_str, (tend-tstart)


def get_performance(res_str):
  for line in res_str.split('\n'):
    if "RANK0 MStencil/s  MAX" in line:
      return float(line.split(':')[1])
  return -1
def set_runtime(kernel, nx, ny, nz, pinning_cmd, pinning_args, fp, params):
  tee(fp, "[AUTO TUNE] ====================================\n[AUTO TUNE] Finding good number of time steps\n")
  nt=32
  prev_perf = 0.
  while(1):
    res, telapsed = run_pluto_test(dry_run=0, kernel=kernel, nx=nx, ny=ny, nz=nz, nt=nt, params=params, pinning_cmd=pinning_cmd, pinning_args=pinning_args, auto_tuning=1)
    cur_perf = get_performance(res)
    perf_err = abs((cur_perf-prev_perf)/cur_perf)
    tee(fp, "[AUTO TUNE] Testing Nt: %d   elapsed_time:%ds   performance variation: %05.2f\n"%(nt, telapsed, 100*perf_err))
    if( (perf_err < 0.01) or (telapsed>100) ):
      return nt
    else:
      nt = nt*2
      prev_perf = cur_perf

def pluto_tuner(kernel, nx, ny, nz, fp):
  from scripts.conf.conf import machine_conf, machine_info
  from scripts.pluto.gen_kernels import param_space 
  from operator import itemgetter
  import time

  # Use pinning without HW counters measurements
  th = machine_info['n_cores']
  pinning_cmd = 'likwid-pin'
  pinning_args = " -q -c S0:0-%d "%(th-1)
  # get representative run time
  nt = set_runtime(kernel, nx, ny, nz, pinning_cmd, pinning_args, fp, params=[32,32,1024])

  param_l = param_space[kernel]
  max_perf = -1

  # prune test cases
  next_nx=0
  param_l = sorted(param_l, key=itemgetter(2))
  for i in param_l:
    if(i[2]>nx): 
      next_nx = i[2]
      break
  if(next_nx!=0): param_l = [ i for i in param_l if i[2]<=next_nx]


  best_param = []
  n_params = len(param_l)
  tstart = time.time()
  for num, param in enumerate(param_l):
    res, telapsed = run_pluto_test(dry_run=0, kernel=kernel, nx=nx, ny=ny, nz=nz, nt=nt, params=param, pinning_cmd=pinning_cmd, pinning_args=pinning_args, auto_tuning=1)
    perf = get_performance(res)
    tee(fp, "[AUTO TUNE] Test %d/%d  elapsed_time:%ds   kernel: %s Nx:%d Ny:%d Nz: %d  Peformance: %08.3f  params:%s\n"% (num+1, n_params, telapsed,  kernel, nx, ny, nz, perf, str(param).strip('[]')) )
    if(perf>max_perf):
      max_perf = perf
      best_param = param
  if(max_perf == -1):
    tee(fp, "Tuner failed\n")
    raise
  tend = time.time()

  tee(fp, "[AUTO TUNE] Best performance - kernel: %s Nx:%d Ny:%d Nz: %d  Peformance: %08.3f  params:%s\n"% (kernel, nx, ny, nz, max_perf, str(best_param).strip('[]')) )
  tee(fp, "[AUTO TUNE] elapsed time: %d", (tend-tstart) )
  return best_param, nt


def igs_test(dry_run, target_dir, exp_name, group='', setup=[]): 
  from scripts.conf.conf import machine_info
  import itertools, os
  from os.path import join as jpath

  target_dir = jpath(os.path.abspath("."),target_dir)

  # Test using rasonable time
  # T = scale * size / perf
  # scale = T*perf/size
#  desired_time = 10
#  if(machine_info['hostname']=='Haswell_18core'):
#    k_perf_order = {'3d25pt':1000, '3d7pt':2500, '3d25pt_var':200, '3d7pt_var':1000}
#  else:
#    k_perf_order = {'3d25pt':700, '3d7pt':1500, '3d25pt_var':200, '3d7pt_var':1000}
#  k_time_scale={}
#  for k, v in k_perf_order.items():
#    k_time_scale[k] = desired_time*v

  points = list(range(32, 5000, 128))
  points = sorted(list(set(points)))

  kernels_limits = {'3d25pt':1057, '3d7pt':1200, '3d25pt_var':545, '3d7pt_var':680}

  if(machine_info['hostname']=='Haswell_18core'):
    kernels_limits = {'3d25pt':1600, '3d7pt':1600, '3d25pt_var':960, '3d7pt_var':1000}

  count=0
  for kernel in ['3d7pt_var']:# '3d25pt', '3d25pt_var', '3d7pt_var']:
    for N in points:
      if (N < kernels_limits[kernel]):
        outfile=('pluto_kernel_%s_N%d_%s_%s.txt' % (kernel, N, group, exp_name[-13:]))
        outfile = jpath(target_dir, outfile)
        if(dry_run==0): fp = open(outfile, 'w')
#        nt = max(int(k_time_scale[kernel]/(N**3/1e6)), 30)
        param =[32,32,1024]
        nt =30
        key = (kernel, N)
        if key not in setup:
          if(dry_run==0): param, nt = pluto_tuner(kernel=kernel, nx=N, ny=N, nz=N, fp=fp)
#         continue

        if(dry_run==0): tee(fp, outfile)
        print outfile
        test_str, telapsed = run_pluto_test(dry_run=dry_run, kernel=kernel, nx=N, ny=N, nz=N, nt=nt, params=param, outfile=outfile)
        if(dry_run==0): tee(fp, test_str)
        if(dry_run==0): fp.close()
        count = count+1
  return count

def main():
  from scripts.utils import create_project_tarball, get_stencil_num
  from scripts.conf.conf import machine_conf, machine_info
  import os
  from csv import DictReader
  import time,datetime

  dry_run = 0

  time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H_%M')
  exp_name = "pluto_increasing_grid_size_at_%s_%s" % (machine_info['hostname'], time_stamp)  

  tarball_dir='results/'+exp_name
  if(dry_run==0): create_project_tarball(tarball_dir, "test_"+exp_name)
  target_dir='results/' + exp_name 

  # parse the results to find out which of the already exist
  data = []
  data_file = os.path.join('results', 'summary.csv')
  try:
    with open(data_file, 'rb') as output_file:
      raw_data = DictReader(output_file)
      for k in raw_data:
        kernel = get_stencil_num(k)
        if(kernel==0):
          k['stencil'] ='3d25pt'
        elif(kernel==1):
          k['stencil'] ='3d7pt'
        elif(kernel==4):
          k['stencil'] ='3d25pt_var'
        elif(kernel==5):
          k['stencil'] ='3d7pt_var'
        else:
          raise
        data.append(k)
  except:
     pass
  setup = set()
  for k in data:
    try:
      setup.add( (k['stencil'], int(k['Global NX'])) )
    except:
      print k
      raise

  #update the pinning information to use all cores
  th = machine_info['n_cores']

  pin_str = "0-%d "%(th-1)

  count = 0
  for group in ['MEM']: #, 'DATA', 'TLB_DATA', 'L2', 'L3', 'ENERGY']:
    if(machine_info['hostname']=='Haswell_18core'):
      machine_conf['pinning_args'] = " -m -g " + group + " -C S1:" + pin_str
    elif(machine_info['hostname']=='IVB_10core'):
      if group=='TLB_DATA': group='TLB' 
      machine_conf['pinning_args'] = " -g " + group + " -C S0:" + pin_str
#    for k in setup: print k
    count = count + igs_test(dry_run, target_dir, exp_name, setup=setup, group=group) 

  print "experiments count =" + str(count)


if __name__ == "__main__":
  main()
