#!/usr/bin/env python
def set_likwid(th, group):
  from scripts.conf.conf import machine_conf, machine_info
  if(machine_info['hostname']=='Haswell_18core'):
    machine_conf['pinning_args'] = " -m -g " + group + " -c " + "%d-%d "%(th, 2*th-1) + '-- numactl --physcpubind=%d-%d'%(th,2*th-1)
  elif(machine_info['hostname']=='IVB_10core'):
    if group=='TLB_DATA': group='TLB' 
    machine_conf['pinning_args'] = " -m -g " + group + " -c " + "%d-%d "%(0, th-1) + '-- numactl --physcpubind=%d-%d'%(0,th-1)


def run_pochoir_test(dry_run, th, kernel, nx, ny, nz, nt, target_dir, outfile, pinning_cmd, pinning_args):
  import os
  import subprocess
  from string import Template
  from scripts.utils import ensure_dir    

  job_template=Template(
"""echo 'OpenMP Threads: $th' | tee $outpath; $pinning_cmd $pinning_args $exec_path $nx $ny $nz $nt | tee $outpath""")

  # set the output path
  target_dir = os.path.join(os.path.abspath("."),target_dir)
  ensure_dir(target_dir)
  outpath = os.path.join(target_dir, outfile)

  # set the executable
  if(kernel==0):
    exec_name ='3dfd'
  elif(kernel==1):
    exec_name ='3d7pt'
  elif(kernel==4):
    exec_name ='3d25pt_var'
  elif(kernel==5):
    exec_name ='3d7pt_var'
  else:
    raise
  exec_path = os.path.join(os.path.abspath("."),exec_name)
 
   
  job_cmd = job_template.substitute(th=th, nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel, outpath=outpath, 
                                    exec_path=exec_path, pinning_cmd=pinning_cmd, pinning_args=pinning_args)
 
  print job_cmd
  if(dry_run==0): sts = subprocess.call(job_cmd, shell=True)

  return job_cmd


def thread_scaling_test(dry_run, target_dir, exp_name, group='', params=[]): 
  from scripts.conf.conf import machine_conf, machine_info
  import itertools


  # Test using rasonable time
  # T = scale * size / perf
  # scale = T*perf/size
  desired_time = 5
  if(machine_info['hostname']=='Haswell_18core'):
    k_perf_order = {0:1000, 1:2500, 4:200, 5:1000}
  else:
    k_perf_order = {0:500, 1:1000, 4:100, 5:600}
  k_time_scale={}
  for k, v in k_perf_order.items():
    k_time_scale[k] = desired_time*v


  radius = {0:4, 1:1, 4:4, 5:1}

  points = {}
  points[0] = [896]
  points[1] = [896]
  points[4] = [768]
  points[5] = [768]

  #update the pinning information to use all cores
  th_max = machine_info['n_cores']

  count=0
  for kernel in [0, 1, 4, 5]:
    for N in points[kernel]:
      for th in list(range(1,1+th_max)):
        set_likwid(th, group)
        outfile=('pochoir_kernel%d_N%d_th%d_%s_%s.txt' % (kernel, N, th, group, exp_name[-13:]))
        nt = max(int(k_time_scale[kernel]/(N**3/1e6)*float(th)/float(th_max)), 30)
        N = N + 2 * radius[kernel] # Pochoir takes the whole size including the halo region
#        print outfile
        run_pochoir_test(dry_run=dry_run, th=th, kernel=kernel, nx=N, ny=N, nz=N, nt=nt, outfile=outfile, target_dir=target_dir, pinning_cmd=machine_conf['pinning_cmd'], pinning_args=machine_conf['pinning_args'])
        count = count+1
  return count

def main():
  from scripts.utils import create_project_tarball, get_stencil_num
  from scripts.conf.conf import machine_conf, machine_info
  import os
  from csv import DictReader
  import time,datetime

  dry_run = 1 if len(sys.argv)<2 else int(sys.argv[1])

  time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H_%M')
  exp_name = "pochoir_thread_scaling_at_%s_%s" % (machine_info['hostname'], time_stamp)  

  tarball_dir='results/'+exp_name
  if(dry_run==0): create_project_tarball(tarball_dir, "project_"+exp_name)
  target_dir='results/' + exp_name 

  # parse the results to find out which of the already exist
  data = []
  data_file = os.path.join('results', 'summary.csv')
  try:
    with open(data_file, 'rb') as output_file:
      raw_data = DictReader(output_file)
      for k in raw_data:
        k['stencil'] = get_stencil_num(k)
        data.append(k)
  except:
     pass
  params = set()
  for k in data:
    try:
      params.add( (k['stencil'], int(k['Global NX'])) )
    except:
      print k
      raise

  count = 0
  for group in ['MEM', 'L2']:# 'TLB_DATA', 'L3', 'DATA', 'ENERGY']:
#    for k in params: print k
    count = count + thread_scaling_test(dry_run, target_dir, exp_name, params=params, group=group) 

  print "experiments count =" + str(count)


if __name__ == "__main__":
  main()
