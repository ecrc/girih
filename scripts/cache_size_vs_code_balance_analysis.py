#!/usr/bin/env python

def load_csv(data_file):
  from csv import DictReader
  with open(data_file, 'rb') as output_file:
    data = DictReader(output_file)
    data = [k for k in data]
  return data

def main():
  import sys
  from ics_utils import get_stencil_num

  raw_data = load_csv(sys.argv[1])

  k_l = set()
  for k in raw_data:
    k_l.add(get_stencil_num(k))
  k_l = list(k_l)

  n_l = set()
  for k in raw_data:
    n_l.add(k['Global NX'])
  n_l = list(n_l)


  for k in k_l:
    for N in n_l:
      gen_res(raw_data, int(k), int(N))

#    gen_res(raw_data)


def gen_res(raw_data, stencil_kernel, N):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import pylab
  from csv import DictWriter
  from operator import itemgetter
  from ics_utils import models, get_bs, get_stencil_num, get_nd

  #fig_width = 8.588*0.393701 # inches
  fig_width = 5.5*0.393701 # inches
  fig_height = 0.68*fig_width #* 210.0/280.0#433.62/578.16

  fig_size =  [fig_width,fig_height]
  params = {
         'axes.labelsize': 5,
         'axes.linewidth': 0.5,
         'lines.linewidth': 1,
         'text.fontsize': 5,
         'legend.fontsize': 5,
         'xtick.labelsize': 5,
         'ytick.labelsize': 5,
         'lines.markersize': 5,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)


  req_fields = [('Total cache block size (kiB)', int), ('MStencil/s  MAX', float), ('Time stepper orig name', str), ('Stencil Kernel semi-bandwidth', int), ('Stencil Kernel coefficients', str), ('Precision', str), ('Time unroll',int), ('Number of time steps',int), ('Number of tests',int), ('Local NX',int), ('Local NY',int), ('Local NZ',int), ('Total Memory Transfer', float), ('Thread group size' ,int), ('Intra-diamond prologue/epilogue MStencils',int), ('Multi-wavefront updates', int), ('Intra-diamond width', int)]
  data = []
  for k in raw_data:

    # Use single field to represent the performance
    if 'Total RANK0 MStencil/s MAX' in k.keys():
      if(k['Total RANK0 MStencil/s MAX']!=''):
        k['MStencil/s  MAX'] = k['MWD main-loop RANK0 MStencil/s MAX']
    # temporary for deprecated format
    if 'RANK0 MStencil/s  MAX' in k.keys():
      if k['RANK0 MStencil/s  MAX']!='':
        k['MStencil/s  MAX'] = k['RANK0 MStencil/s  MAX']

    tup = dict()
    # add the general fileds
    for f in req_fields:
      try:
        tup[f[0]] = map(f[1], [k[f[0]]] )[0]
      except:
        print f[0]
    # add the stencil operator
    tup['Kernel'] = get_stencil_num(k)
    data.append(tup)

    # add precision information
    tup['word size'] = 8 if k['Precision'] in 'DP' else 4

  #for i in data: print i

  data2 = []
  for tup in data:
    tup['Actual Bytes/LUP'] = actual_BpU(tup)
    tup['Model'] = models(tup)
    # model error
    tup['Err %'] = 100 * (tup['Model'] - tup['Actual Bytes/LUP'])/tup['Actual Bytes/LUP']
    tup['D_width'] = tup['Intra-diamond width']
    tup['Performance'] = tup['MStencil/s  MAX']
    tup['Cache block'] = get_bs(Dw=tup['D_width'], Nd=get_nd(tup['Kernel']), Nf=tup['Multi-wavefront updates'], Nx=tup['Local NX'], WS=tup['word size'], R=tup['Stencil Kernel semi-bandwidth'])
    print "cache block C:", tup['Total cache block size (kiB)']/1024.0, " Python:", tup['Cache block']
    data2.append(tup)
#    try: print "%6.3f  %6.3f  %6.3f" % (tup['Cache block'], tup['Total cache block size (kiB)']/1024.0,tup['Cache block']- tup['Total cache block size (kiB)']/1024.0)
#    except: pass

  #for i in data2: print i
  data2 = sorted(data2, key=itemgetter('Kernel', 'Local NX', 'D_width'))


  cs=[]
  cb=[]
  cb_meas=[]
  Dw=[]
  for k in data2:
    if k['Kernel']==stencil_kernel and k['Local NX']==N:
      cs.append(k['Cache block'])
      cb.append(k['Model'])
      cb_meas.append(k['Actual Bytes/LUP'])
      Dw.append(k['D_width'])

  #for i in range(len(cs)):
  #  print Dw[i], cs[i], cb_meas[i], cb[i]

  if Dw==[]: return

  fig, ax = plt.subplots()
  ax.plot(cs, cb     , marker='+', linestyle='-', color='k', label="Model")
  ax.plot(cs, cb_meas, marker='x', linestyle='--', color='b', label="Measured")

  # show the usable cache size limits
  ax.plot([12.5, 12.5], [0, 0.7*cb[0]], linestyle='-', color='r', label="Usable cache size")

  ax.set_ylabel('Code balance (Bytes/LUP)')
  ax.set_xlabel('Cache block size (MiB) PER THREAD')
  ax.set_ylim([0, max(cb_meas+cb)+1])
  ax.set_xlim([0, max(cs)+0.5])
  ax2 = ax.twiny()
  ax2.set_xticks(cs)
  ax2.set_xlabel('Diamond width')
  ax2.set_xlim(ax.get_xlim())

  if stencil_kernel==1:
    Dw = map(str,Dw)
    Dw[1]=''
    Dw[3]=''
    Dw[5]=''
  ax2.set_xticklabels(Dw)

#  for i, d in enumerate(Dw):
    #if ((d+4)%8 == 0):
#    ax.annotate(d, (cs[i], cb[i]))  

  title = '_code_balance_vs_cache_size_N'+str(N)
  if stencil_kernel == 0:
      title = '25_pt_const' + title
  elif stencil_kernel == 1:
      title = '7_pt_const' + title
  elif stencil_kernel == 4:
      title = '25_pt_var' + title
  elif stencil_kernel == 5:
      title = '7_pt_var' + title

  ax.legend(loc='best')
  ax.grid()
  pylab.savefig(title+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
  plt.clf()



def actual_BpU(tup):
  total_mem = tup['Total Memory Transfer']
  R = tup['Stencil Kernel semi-bandwidth']
  nt = tup['Number of time steps']

  nx = tup['Local NX']
  ny = tup['Local NY']
  nz = tup['Local NZ']

  oh = 0#tup['Intra-diamond prologue/epilogue MStencils']

  stencil_size = nx*ny*nz#2*ny*nz + ny*nz*(nx+2*R) 
  BpU = (total_mem * 10**9) / ( (stencil_size * nt - oh*10**6)*tup['Number of tests'])

  #print BpU, total_mem, stencil_size, nt, oh, tup['Number of tests']
  return BpU


if __name__ == "__main__":
  main()

