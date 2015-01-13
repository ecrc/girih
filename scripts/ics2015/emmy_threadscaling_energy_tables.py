#!/usr/bin/env python
def main():
  import sys
  from utils import select_fields, load_csv
  raw_data = load_csv(sys.argv[1])

  stencil='7_pt_const'
  rows = [
           {'Thread group size':'0' , 'Stencil Kernel coefficients':'constant', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'6' },
           {'Thread group size':'1' , 'Stencil Kernel coefficients':'constant', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'10'},
           {'Thread group size':'2' , 'Stencil Kernel coefficients':'constant', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'10'},
           {'Thread group size':'5' , 'Stencil Kernel coefficients':'constant', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'10'},
           {'Thread group size':'10', 'Stencil Kernel coefficients':'constant', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'10'}
         ]
  create_table(raw_data, rows, stencil)


  stencil='7_pt_var'
  rows = [
           {'Thread group size':'0' , 'Stencil Kernel coefficients':'variable no-symmetry', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'6'},
           {'Thread group size':'1' , 'Stencil Kernel coefficients':'variable no-symmetry', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'8'},
           {'Thread group size':'2' , 'Stencil Kernel coefficients':'variable no-symmetry', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'10'},
           {'Thread group size':'5' , 'Stencil Kernel coefficients':'variable no-symmetry', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'10'},
           {'Thread group size':'10', 'Stencil Kernel coefficients':'variable no-symmetry', 'Stencil Kernel semi-bandwidth':'1', 'OpenMP Threads':'10'}
         ]
  create_table(raw_data, rows, stencil)


  stencil='25_pt_var'
  rows = [
           {'Thread group size':'0' , 'Stencil Kernel coefficients':'variable axis-symmetric', 'Stencil Kernel semi-bandwidth':'4', 'OpenMP Threads':'8' },
           {'Thread group size':'1' , 'Stencil Kernel coefficients':'variable axis-symmetric', 'Stencil Kernel semi-bandwidth':'4', 'OpenMP Threads':'9'},
           {'Thread group size':'2' , 'Stencil Kernel coefficients':'variable axis-symmetric', 'Stencil Kernel semi-bandwidth':'4', 'OpenMP Threads':'8'},
           {'Thread group size':'5' , 'Stencil Kernel coefficients':'variable axis-symmetric', 'Stencil Kernel semi-bandwidth':'4', 'OpenMP Threads':'10'},
           {'Thread group size':'10', 'Stencil Kernel coefficients':'variable axis-symmetric', 'Stencil Kernel semi-bandwidth':'4', 'OpenMP Threads':'10'}
         ]
  create_table(raw_data, rows, stencil)



def create_table(raw_data, rows, stencil):

  from csv import DictWriter
  from utils import select_fields, load_csv
  from ics_utils import models, get_stencil_num
  from operator import itemgetter


  cols_format = [('Time stepper orig name', str), ('Stencil Kernel coefficients', str), ('Thread group size', int), ('Stencil Kernel semi-bandwidth', int), ('OpenMP Threads', int), ('Energy', float), ('Energy DRAM', float), ('Power', float), ('Power DRAM', float), ('MStencil/s  MAX', float), ('Global NX', int), ('Local NY', int), ('Global NY', int), ('Global NZ', int), ('Number of time steps', int), ('Number of tests', int), ('Intra-diamond prologue/epilogue MStencils', int), ('Total cache block size (kiB)', int), ('Block size in X', int), ('Precision', str), ('Time unroll',int), ('Intra-diamond width', int), ('Multi-wavefront updates', int)]

  cols = [f[0] for f in cols_format]

  data  =  select_fields(raw_data, rows, cols)

  for k in data:
    if k['Block size in X'] == '':
      k['Block size in X'] = '100000'
    for val, fmt in cols_format:
      try:
        k[val] = map(fmt, [k[val]])[0]
      except:
        print val, k[val]
  # compute derived values
  for k in data:
    lups = k['Number of tests'] * (k['Global NX']*k['Global NY']*k['Global NZ']*k['Number of time steps'] - k['Intra-diamond prologue/epilogue MStencils'])

    k['pJ/LUP CPU'] = k['Energy']/lups*1e9
    k['pJ/LUP DRAM'] = k['Energy DRAM']/lups*1e9
    k['pJ/LUP Total'] = k['pJ/LUP CPU'] + k['pJ/LUP DRAM']
    k['Power CPU'] = k['Power']
    k['Power Total'] = k['Power CPU'] + k['Power DRAM']

    k['Threads'] = k['OpenMP Threads']
    k['MLUP/s'] = k['MStencil/s  MAX']

    k['Thread group size'] = k['Thread group size'] 
    tgs = k['Thread group size'] 
    if tgs == 0:
      k['Method'] = 'Spt. Blk.'
    else:
      k['Method'] = str(k['Thread group size'])+'WD'

    k['Kernel'] = get_stencil_num(k)
    k['Model Bytes/LUP'] = models(k)

    k['Dw|Nf'] = str(k['Intra-diamond width'])+'|'+str(k['Multi-wavefront updates'])

    k['Cache blk. [MiB]'] = k['Total cache block size (kiB)']/1024
    if k['Block size in X'] < k['Global NX']:
      k['Cache blk. [MiB]'] = k['Cache blk. [MiB]'] * k['Global NX'] / k['Block size in X']

  data = sorted(data, key=itemgetter('Thread group size'))

  fields = ['Method', 'Threads','MLUP/s', 'Cache blk. [MiB]', 'Model Bytes/LUP', 'Dw|Nf', 'Power CPU', 'Power DRAM', 'Power Total', 'pJ/LUP CPU', 'pJ/LUP DRAM', 'pJ/LUP Total']
  with open(stencil+'_threadscaling_table.csv', 'w') as output_file:
    r = DictWriter(output_file, fieldnames=fields)
    r.writeheader()
    for k in data:
      k2 = dict()
      for f in k.keys():
        for f2 in fields:
          if f == f2:
            k2[f] = k[f]
      r.writerow(k2)

if __name__ == "__main__":
  main()

