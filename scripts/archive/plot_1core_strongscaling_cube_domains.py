#!/usr/bin/env python
def main():
  import sys

  raw_data = load_csv(sys.argv[1])

 #   for ts in ['Naive', 'Dynamic-Intra-Diamond']  
 #   for k in [0, 1, 2]:
 #       for is_dp in [0,1]:
  plot_lines(raw_data, 'perf')
  plot_lines(raw_data, 'BW')

        
def plot_lines(raw_data, measurment):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import matplotlib
  import pylab

  req_fields = [('Global NX', 1), ('MStencil/s  MAX', 2), ('LIKWID performance counter',0), 
    ('Sustained Memory BW', 2), ('Total Memory Transfer', 2), 
    ('Sustained L3 BW', 2), ('Total L3 Transfer',2),
    ('Sustained L2 BW', 2), ('Total L2 Transfer',2),
    ('Number of time steps',1), ('Number of tests',1)]

  data = []
  for k in raw_data:
    tup = {}
    # add the general fileds
    for f in req_fields:
      v = k[f[0]]
      if f[1]==1: v = int(k[f[0]]) 
      if f[1]==2: 
        v = 0 if k[f[0]] == '' else float(k[f[0]]) 
      tup[f[0]] = v
    data.append(tup)

  data = sorted(data, key=itemgetter('Global NX', 'LIKWID performance counter'))
  #for i in data: print i


  fig, ax1 = plt.subplots()
  lns = []
  #for perfctr in ['MEM', 'L3', 'L2']:
  x = []
  y = []
  x_mem_bw = []
  y_mem_bw = []
  x_l3_bw = []
  y_l3_bw = []
  for k in data:
    if k['LIKWID performance counter']=='MEM' or k['LIKWID performance counter']=='L3':
      x.append(k['Global NX'])
      y.append(k['MStencil/s  MAX'])
    if k['LIKWID performance counter']=='L3':
      N = k['Global NX']
      tot_LUPs = (2*N**2 + N**2*(N+2)) * k['Number of time steps'] * k['Number of tests']
      bw = k['Total L3 Transfer'] * 1e9 / tot_LUPs
      x_l3_bw.append(N)
      y_l3_bw.append(bw)
    if k['LIKWID performance counter']=='MEM':
      N = k['Global NX']
      tot_LUPs = (2*N**2 + N**2*(N+2)) * k['Number of time steps'] * k['Number of tests']
      bw = k['Total Memory Transfer'] * 1e9 / tot_LUPs
      x_mem_bw.append(N)
      y_mem_bw.append(bw)



  lns = lns + ax1.plot(x, y, color='k', linestyle='-', label='Performance')
  plt.ylim([0, 1100])

  if measurment == 'perf':
    # add limits
    proc_speed = 2.2*1e9

   # Overlap cache with mem access
    ECM_MEM =  proc_speed / (34/16) / ( 1024 * 1024)
    lns = lns + ax1.plot([128, max(x)], [ECM_MEM, ECM_MEM], color='g', linestyle='-', label='Cache/Mem access overlap')

    # N > 1024
    ECM_MEM_2l =  proc_speed / (48.1/16) / ( 1024 * 1024)
    lns = lns + ax1.plot([1024, max(x)], [ECM_MEM_2l, ECM_MEM_2l], color='r', linestyle='-', label='2 layers L3-MM')
    # 128 < N < 1024
    ECM_L3 =  proc_speed / (44.6/16) / ( 1024 * 1024)
    lns = lns + ax1.plot([128, 1024], [ECM_L3, ECM_L3], color='b', linestyle='-', label='1 layer L3-MM')
    # 45 < N < 128
    ECM_L2 =  proc_speed / (40.6/16) / ( 1024 * 1024)
    lns = lns + ax1.plot([45, 128], [ECM_L2, ECM_L2], color='m', linestyle='-', label='1 layer L2-L3')
  elif measurment == 'BW':
    # Plot the achieved Bytes/LUP at different cache levels
    ax2 = ax1.twinx()
    lns = lns + ax2.plot(x_mem_bw, y_mem_bw, color='r',  marker='*', linestyle='-', label='obtained L3-MEM')
    lns = lns + ax2.plot(x_l3_bw, y_l3_bw, color='r', marker='x', linestyle='-', label='obtained L2-L3')

    # Plot predicted Bytes/LUP
    # N > 1024
    lns = lns + ax2.plot([1024, max(x)], [16, 16], color='r', linestyle='-', label='L3-MEM')
    # N < 1024
    ax2.plot([0, 1024], [12, 12], color='r', linestyle='-', label='L3-MEM')

    # N < 128
    lns = lns + ax2.plot([0, 128], [12, 12], color='b', linestyle='--', label='L2-L3')
    # N > 128
    ax2.plot([128, 256], [20, 20], color='b', linestyle='--', label='L2-L3')

    ax2.set_ylabel('Bytes/LUP', color='r')
    for tl in ax2.get_yticklabels():
      tl.set_color('r')

    plt.sca(ax2)
    plt.yticks([0,4,8,12,16,20,24])
 

  title = 'Cube domain Ivy Bridge 1-core strongscaling, 7_pt const. stencil SP ' + measurment
  f_name = title
  ax1.set_xlabel('Size in each dimension')
  ax1.set_ylabel('MLUP/s')
  plt.title(title)

  labs = [l.get_label() for l in lns]
  ax1.legend(lns, labs, loc='best')

  ax1.grid()
  pylab.savefig(f_name+'.png')
  pylab.savefig(f_name+'.pdf', format='pdf')

  #plt.show()     
  plt.clf()
      
def load_csv(data_file):
  from csv import DictReader
  with open(data_file, 'rb') as output_file:
    data = DictReader(output_file)
    data = [k for k in data]
  return data
    
    
if __name__ == "__main__":
  main()
