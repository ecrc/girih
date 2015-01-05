#!/usr/bin/env python
def main():
  import sys

  raw_data = load_csv(sys.argv[1])
  create_table(raw_data)

def get_stencil_num(k):
  # add the stencil operator
  if  k['Stencil Kernel coefficients'] in 'constant':
    if  int(k['Stencil Kernel semi-bandwidth'])==4:
      stencil = 0
    else:
      stencil = 1
  elif  'no-symmetry' in k['Stencil Kernel coefficients']:
    stencil = 5
  elif  'sym' in k['Stencil Kernel coefficients']:
    if int(k['Stencil Kernel semi-bandwidth'])==1:
      stencil = 3
    else:
      stencil = 4
  else:
    stencil = 2
  return stencil


def sort_cols(data,items):
  items2 = []
  for i in items:
    if i in data:
      data.remove(i)
      items2.append(i)
  data = items2 + list(data)
  return data


def create_table(raw_data):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import pylab
  from csv import DictWriter


  req_fields = [('Thread group size', int), ('WD main-loop RANK0 MStencil/s  MAX', float),
                ('Time stepper orig name', str), ('OpenMP Threads', int),
                ('MStencil/s  MAX', float), ('Time unroll',int),
                ('Intra-diamond prologue/epilogue MStencils', int), ('Global NX', int),
                ('Number of time steps', int), ('Number of tests', int),
                ('Total Memory Transfer', float), ('L2 data volume sum', float),
                ('L3 data volume sum', float), ('LIKWID performance counter',str)]
  data = []
  for k in raw_data: 
    tup = {}
    # add the general fileds
    for f in req_fields:
      if k[f[0]] != '':
        try:
          tup[f[0]] = map(f[1], [k[f[0]]] )[0]
        except:
          print k[f[0]], f[0], f[1] 
          raise
    # add the stencil operator
    tup['stencil'] = get_stencil_num(k)
    
    # add the precision information
    if k['Precision'] in 'DP':
      p = 1
    else:
      p = 0
    tup['Precision'] = p
    data.append(tup)
#    for i in data: print i


  # aggrigate the data
  tgs_l = set()
  kernel_l = set()
  th_l = set()
  group_l = set()
  size_l = set()
  for k in data:
    tgs_l.add(k['Thread group size'])
    kernel_l.add(k['stencil'])
    th_l.add(k['OpenMP Threads'])
    group_l.add(k['LIKWID performance counter'])
    size_l.add(k['Global NX'])

  data2 = []
  for tgs in tgs_l:
    for kernel in kernel_l:
      for size in size_l:
        for th in th_l:
          k = {}
          k['Thread group size'] = tgs
          k['stencil'] = kernel
          k['N'] = size
          k['OpenMP Threads'] = th
          for tup in data:
            try:
              if (tup['Thread group size']==tgs ) and   (tup['stencil']==kernel ) and (tup['Global NX']==size ) and (tup['OpenMP Threads']==th ):
                if tup['LIKWID performance counter'] == 'MEM':
                  k['Total Memory Transfer'] = tup['Total Memory Transfer']
                elif tup['LIKWID performance counter'] == 'L2':
                  k['L2 data volume sum'] = tup['L2 data volume sum']
                elif tup['LIKWID performance counter'] == 'L3':
                  k['L3 data volume sum'] = tup['L3 data volume sum']
                k[tup['LIKWID performance counter']+'_perf'] = tup['MStencil/s  MAX']
                k['Number of time steps'] = tup['Number of time steps']
                k['Intra-diamond prologue/epilogue MStencils'] = tup['Intra-diamond prologue/epilogue MStencils']
                k['Number of tests'] = tup['Number of tests']
                k['D width'] = (tup['Time unroll']+1)*2
            except:
              print tup
              raise
            
          if 'Total Memory Transfer' in k.keys():
            data2.append(k)


  for k in data2:
    glups = (k['Number of time steps'] * k['N']**3 - k['Intra-diamond prologue/epilogue MStencils']*10**6 ) * k['Number of tests'] / 10**9

    k['MEM Bytes/LUP'] = k['Total Memory Transfer']/glups
    k['L2 Bytes/LUP']  = k['L2 data volume sum']/glups
    k['L3 Bytes/LUP']  = k['L3 data volume sum']/glups

  from operator import itemgetter
  data2 = sorted(data2, key=itemgetter('stencil', 'Thread group size', 'N', 'OpenMP Threads'))

  fields = sort_cols(data2[0].keys(), ['stencil', 'Thread group size', 'N', 'OpenMP Threads', 'D width', 'MEM Bytes/LUP', 'L3 Bytes/LUP', 'L2 Bytes/LUP'])
  with open('cache_transfer.csv', 'w') as output_file:
    r = DictWriter(output_file, fieldnames=fields)
    r.writeheader()
    for k in data2:
      k2 = dict()
      for f in k.keys():
        for f2 in fields:
          if f == f2:
            k2[f] = k[f]
      r.writerow(k2)


def load_csv(data_file):
  from csv import DictReader
  with open(data_file, 'rb') as output_file:
    data = DictReader(output_file)
    data = [k for k in data]
  return data
  
    
if __name__ == "__main__":
  main()

