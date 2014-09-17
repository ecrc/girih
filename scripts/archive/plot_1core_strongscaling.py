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
  import matplotlib

  req_fields = [('Global NY', 1), ('Global NX', 1), ('Global NZ', 1), ('MStencil/s  MAX', 2), ('LIKWID performance counter',0), 
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

  data = sorted(data, key=itemgetter('Global NY', 'LIKWID performance counter'))
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
      x.append(k['Global NY'])
      y.append(k['MStencil/s  MAX'])
    if k['LIKWID performance counter']=='L3':
      ny = k['Global NY']
      nx = k['Global NX']
      nz = k['Global NZ']
      tot_LUPs = (2*ny*nz + ny*nz*(nx+2)) * k['Number of time steps'] * k['Number of tests']
      bw = k['Total L3 Transfer'] * 1e9 / tot_LUPs
      x_l3_bw.append(ny)
      y_l3_bw.append(bw)
    if k['LIKWID performance counter']=='MEM':
      ny = k['Global NY']
      nx = k['Global NX']
      nz = k['Global NZ']
      tot_LUPs = (2*ny*nz + ny*nz*(nx+2)) * k['Number of time steps'] * k['Number of tests']
      bw = k['Total Memory Transfer'] * 1e9 / tot_LUPs
      x_mem_bw.append(ny)
      y_mem_bw.append(bw)



  #lns = lns + ax1.semilogx(x, y, color='k', linestyle='-', label='Performance', basex=2)
  lns = lns + ax1.plot(x, y, color='k', linestyle='-', label='Performance')
  plt.ylim([0,1100])

  lns1 =[]
  lns2 =[]
  if measurment == 'perf':
    # add limits
    proc_speed = 2.2*1e9

    # ny > 4096
    ECM_MEM_3l =  proc_speed / (51.6/16) / ( 1024 * 1024)
    lns2 = lns2 + ax1.plot([4096, max(x)], [ECM_MEM_3l, ECM_MEM_3l], color='g', linestyle='-', label='3 layers L3-MM')

    #  2048 < ny < 4096
    ECM_MEM_2l =  proc_speed / (48.1/16) / ( 1024 * 1024)
    lns2 = lns2 + ax1.plot([2048, 4096], [ECM_MEM_2l, ECM_MEM_2l], color='r', linestyle='-', label='2 layers L3-MM')

    # 64 < ny < 2048
    # pessimistic bounds
    ECM_MEM_1l =  proc_speed / (44.6/16) / ( 1024 * 1024)
    lns2 = lns2 + ax1.plot([256, 2048], [ECM_MEM_1l, ECM_MEM_1l], color='b', linestyle='-', label='1 layer L3-MM')
    # Overlap cache with mem access
    ECM_MEM_1l =  proc_speed / (34/16) / ( 1024 * 1024)
    lns2 = lns2 + ax1.plot([256, max(x)], [ECM_MEM_1l, ECM_MEM_1l], color='m', linestyle='-', label='Cache/Mem access overlap')
 
    # 32 < ny < 64
    ECM_L3_1l =  proc_speed / (42.6/16) / ( 1024 * 1024)
    lns1 = lns1 + ax1.plot([32, 64], [ECM_L3_1l, ECM_L3_1l], color='m', linestyle='-', label='2 layers L2-L3')
    # ny < 32
    ECM_L2 =  proc_speed / (40.6/16) / ( 1024 * 1024)
    lns1 = lns1 + ax1.plot([0, 32], [ECM_L2, ECM_L2], color='c', linestyle='-', label='1 layer L2-L3')


  elif measurment == 'BW':
    # Plot the achieved Bytes/LUP at different cache levels
    ax2 = ax1.twinx()
    lns2 = lns2 + ax2.plot(x_mem_bw, y_mem_bw, color='r',  marker='*', linestyle='-', label='obtained L3-MEM')
    lns1 = lns1 + ax2.plot(x_l3_bw, y_l3_bw, color='r', marker='x', linestyle='-', label='obtained L2-L3')

    # Plot predicted Bytes/LUP
    # ny > 4096
    lns2 = lns2 + ax2.plot([4096, max(x)], [20, 20], color='b', linestyle='-', label='3 layers L3-MEM')
    # 2048 < ny < 4096
    lns2 = lns2 + ax2.plot([2048, 4096], [16, 16], color='g', linestyle='-', label='2 layers L3-MEM')
    # 64 < ny < 2048
    lns2 = lns2 + ax2.plot([256, 2048], [12, 12], color='m', linestyle='-', label='1 layer L3-MEM')
    lns1 = lns1 + ax2.plot([64, 256], [20, 20], color='b', linestyle='--', label='3 layers L2-L3')
    # 32 < ny < 64
    lns1 = lns1 + ax2.plot([32, 64], [16, 16], color='g', linestyle='--', label='2 layers L2-L3')
    # ny < 32
    lns1 = lns1 + ax2.plot([0, 32], [12, 12], color='m', linestyle='--', label='1 layer L2-L3')

    ax2.set_ylabel('Bytes/LUP', color='r')
    for tl in ax2.get_yticklabels():
      tl.set_color('r')
    plt.sca(ax2)
    plt.yticks([0,4,8,12,16,20,24])
 


  ax1.set_xlabel('Size in Y')
  ax1.set_ylabel('MLUP/s')

  ax1.grid()
 
  title = 'nx=512 fixed domain Ivy Bridge 1-core strongscaling 7_pt const. stencil SP ' + measurment
  plt.title(title)

  plt.xlim([257,max(x)])
  labs = [l.get_label() for l in lns + lns2]
  ax1.legend(lns+lns2, labs, loc='best',prop={'size':9})
  pylab.savefig(title+'.png')
  pylab.savefig(title+'.pdf', format='pdf')


  plt.xlim([0,256])
  labs = [l.get_label() for l in lns + lns1]
  ax1.legend(lns+lns1, labs, loc='best',prop={'size':9})
  pylab.savefig(title+'_with_limits.png')
  pylab.savefig(title+'_with_limits.pdf', format='pdf')

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
