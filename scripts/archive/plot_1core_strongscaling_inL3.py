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
  import pylab

  req_fields = [('Global NY', 1), ('Global NX', 1), ('Global NZ', 1), ('MStencil/s  MAX', 2), ('LIKWID performance counter',0), 
    ('Sustained Memory BW', 2), ('Total Memory Transfer', 2), 
#    ('Sustained L3 BW', 2), ('Total L3 Transfer',2),
#    ('Sustained L2 BW', 2), ('Total L2 Transfer',2),
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

  data = sorted(data, key=itemgetter('Global NY', 'LIKWID performance counter'))
  #for i in data: print i


  fig, ax1 = plt.subplots()
  lns = []
  x = []
  y = []
  y_mem_bw = []
  for k in data:
    x.append(k['Global NY'])
    y.append(k['MStencil/s  MAX'])
    y_mem_bw.append(k['Sustained Memory BW']/1024.0)


  lns = lns + ax1.plot(x, y, color='k', linestyle='-', label='Performance')

  proc_speed = 2.2*1e9
  ECM_inL3 =  proc_speed / ((17.7+10+10)/16) / ( 1024 * 1024)
  lns = lns + ax1.plot([0, max(x)], [ECM_inL3, ECM_inL3], color='k', linestyle='--', label='ECM in L3')
  plt.ylim([0, 1000])

  ax2 = ax1.twinx()
  lns = lns + ax2.plot(x, y_mem_bw, color='r', marker='x', linestyle='-', label='obtained L3-MEM BW')

  ax1.set_xlabel('Size in Y')
  ax1.set_ylabel('MLUP/s')

  ax1.grid()
 
  title = 'In L3 problem_nx=512_nz=32_Ivy Bridge 1-core strongscaling, 7_pt const. stencil SP '
  plt.title(title)

  ax2.set_ylabel('GB/s', color='r')
  for tl in ax2.get_yticklabels():
    tl.set_color('r')
  plt.sca(ax2)

  labs = [l.get_label() for l in lns]
  ax1.legend(lns, labs, loc='best',prop={'size':9})
  pylab.savefig(title+'.png')

#  plt.show()     
  plt.clf()
      
def load_csv(data_file):
  from csv import DictReader
  with open(data_file, 'rb') as output_file:
    data = DictReader(output_file)
    data = [k for k in data]
  return data
    
    
if __name__ == "__main__":
  main()
