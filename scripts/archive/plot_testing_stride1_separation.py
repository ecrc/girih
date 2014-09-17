#!/usr/bin/env python
def main():
  import sys

  raw_data = load_csv(sys.argv[1])

 #   for ts in ['Naive', 'Dynamic-Intra-Diamond']  
 #   for k in [0, 1, 2]:
 #       for is_dp in [0,1]:
  plot_lines(raw_data, 'Naive')
  plot_lines(raw_data, 'Diamond')

        
def plot_lines(raw_data, ts):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import pylab


  req_fields = [('Time stepper orig name', 0), ('Using separate call to 1-stride loop',1), ('Global NY', 1), ('Global NX', 1), ('Global NZ', 1), ('MStencil/s  MAX', 2), 
                ('Time unroll',1)]

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

    # add the stencil operator
    if  int(k['Stencil Kernel semi-bandwidth'])==4:
      stencil = 0
    elif  k['Stencil Kernel coefficients'] in 'constant':
      stencil = 1
    else:
      stencil = 2
    tup['kernel'] = stencil

    data.append(tup)

  data = sorted(data, key=itemgetter('Time stepper orig name', 'kernel', 'Using separate call to 1-stride loop', 'Global NY'))
  #for i in data: print i

  k_n = ['25-pt-const', '7-pt-const', '7-pt-var']
  sp_n = ['w/o split', 'w/ split']
  k_col = ['k', 'b', 'r']
  sp_l = ['-', '--']
  fig, ax1 = plt.subplots()
  for kernel in [0,1,2]:
    for split in [0, 1]:
      lns = []
      x = []
      y = []
      for k in data:
        if ts in k['Time stepper orig name'] and k['kernel']==kernel and k['Using separate call to 1-stride loop'] == split:
          x.append(k['Global NY'])
          y.append(k['MStencil/s  MAX'])
      if x:
        lns = lns + ax1.plot(x, y, color=k_col[kernel], linestyle=sp_l[split], label='%s_%s'%(k_n[kernel], sp_n[split]))


    
  ax1.set_xlabel('Size in each dimension')
  ax1.set_ylabel('MLUP/s')

  ax1.grid()
 
  title = 'Splitting vs. no splitting of stride-1 loop in %s TS' % ts
  plt.title(title)

  ax1.legend(loc='best',prop={'size':9})
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
