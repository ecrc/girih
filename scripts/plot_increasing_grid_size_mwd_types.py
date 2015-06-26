#!/usr/bin/env python
def main():
  import sys
  from scripts.utils import get_stencil_num, load_csv
  raw_data = load_csv(sys.argv[1])

  k_l = set()
#  mwdt_l = set()
  prec_l = set()
  prof_mem = True 
  prof_energy = True 
  for k in raw_data:

    # Use single field to represent the performance
    if 'Total RANK0 MStencil/s MAX' in k.keys():
      k['MStencil/s  MAX'] = k['MWD main-loop RANK0 MStencil/s MAX'] 
    # temporary for deprecated format
    if 'RANK0 MStencil/s  MAX' in k.keys():
      if k['RANK0 MStencil/s  MAX']!='':
        k['MStencil/s  MAX'] = k['RANK0 MStencil/s  MAX'] 


    # add stencil operator
    k['stencil'] = get_stencil_num(k)
    k_l.add(k['stencil'])

    # add mwd type
    mwd = k['Wavefront parallel strategy'].lower()
    k['mwdt']='none'
    if('fixed' in mwd) and ('relaxed' in mwd):
      k['mwdt'] = 'fers'
    elif('fixed' in mwd):
      k['mwdt'] = 'fe'
    elif('relaxed' in mwd):
      k['mwdt'] = 'rs'
    elif('wavefront' in mwd):
      k['mwdt'] = 'block'
#d    mwdt_l.add(k['mwdt'])

    # disable bandwidth and energy plots if missing from any entry
    if 'Sustained Memory BW' in k.keys():
      if k['Sustained Memory BW'] =='':
        prof_mem=False
    else: 
      prof_mem=False
    if 'Energy' not in k.keys(): prof_energy=False

    # add the approach
    if int(k['Thread group size']) == 0:
      k['method'] = 'Spt.blk.'
    elif ( (int(k['Thread group size']) == 1) and (int(k['Threads along z-axis']) == 0) ):
      k['method'] = 'CATS2'
    else:
      k['method'] = 'MWD'

    # add precision information
    p = 1 if k['Precision'] in 'DP' else 0
    k['Precision'] = p
    prec_l.add(p)
  k_l = list(k_l)
#  mwdt_l = list(mwdt_l)
  prec_l = list(prec_l)

  for k in k_l:
#    for mwdt in mwdt_l:
    for is_dp in prec_l:
      plot_lines(raw_data, k, is_dp, prof_mem=prof_mem, prof_energy=prof_energy)

        
def plot_lines(raw_data, stencil_kernel, is_dp, prof_mem, prof_energy):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import matplotlib
  import pylab
  from pylab import arange,pi,sin,cos,sqrt
  from scripts.utils import get_stencil_num

  m = 3.0
  fig_width = 4.0*0.393701*m # inches
  fig_height = 1.0*fig_width #* 210.0/280.0#433.62/578.16

  fig_size =  [fig_width,fig_height]
  params = {
         'axes.labelsize': 6*m,
         'axes.linewidth': 0.25*m,
         'lines.linewidth': 0.75*m,
         'font.size': 7*m,
         'legend.fontsize': 4*m,
         'xtick.labelsize': 6*m,
         'ytick.labelsize': 6*m,
         'lines.markersize': 1,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)

  

  req_fields = [('method', str), ('MStencil/s  MAX', float), ('mwdt', str), ('Precision', int), ('stencil', int), ('Intra-diamond prologue/epilogue MStencils', int), ('Global NX', int), ('Number of time steps', int), ('Number of tests', int), ('Thread group size', int)]
  if(prof_mem): req_fields = req_fields + [ ('Sustained Memory BW', float)]
  if(prof_energy):
    req_fields = req_fields + [('Energy', float), ('Energy DRAM', float), ('Power',float), ('Power DRAM', float)]
  
  data = []
  for k in raw_data:
    tup = {}
    # add the general fileds
    for f in req_fields:
      try:
        tup[f[0]] = map(f[1], [k[f[0]]] )[0]
      except:
        print("ERROR: results entry missing essential data at file:%s"%(k['file_name']))
        print f[0]
        print k
        return
    data.append(tup)


  method_l = set()
  mwdt_l = set()
  for k in data: 
    method_l.add(k['method'])
    if (k['mwdt']!=''): mwdt_l.add(k['mwdt'])
  method_l = list(method_l)
  mwdt_l = list(mwdt_l)
  


  data = sorted(data, key=itemgetter('method', 'mwdt', 'Global NX'))
#  for i in data: print i


  # compute pJ/LUP
  if(prof_energy):
    for k in data:
      ext_lups = k['Intra-diamond prologue/epilogue MStencils']*1e6
      NT = k['Global NX']**3 * k['Number of time steps']
      Neff = NT - ext_lups
      N = k['Number of tests'] * Neff / 1e9
      k['cpu energy pj/lup'] = k['Energy']/N
      k['dram energy pj/lup'] = k['Energy DRAM']/N
      k['total energy pj/lup'] = k['cpu energy pj/lup'] + k['dram energy pj/lup']


  max_x = 0
  max_y = 0

  marker_s = 7
  line_w = 1
  line_s = '-' 
  cols = 'kbgmrcy'
  markers = '+xo^vx*'
  idx=0
  for method in method_l:

    mwdt_str_l = mwdt_l if(method=='MWD') else ['']

    for mwdt_str in mwdt_str_l:
      marker = markers[idx]
      col = cols[idx]
      x = []
      y = []
      y_m = []
      y_e = []
      for k in data:
        if ( (k['method'] == method) and ( k['method']!='MWD' or (k['method']=='MWD' and k['mwdt'] == mwdt_str) ) and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp) ):
          if(prof_mem): y_m.append(k['Sustained Memory BW']/10**3)
          if(prof_energy): y_e.append(k['total energy pj/lup'])
          x.append(k['Global NX'])
          y.append(k['MStencil/s  MAX']/10**3)
      ts2 = method + mwdt_str
      if(x):
        plt.figure(0)
        plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=ts2)
        max_y = max(max(y), max_y)
        max_x = max(max(x), max_x)
      if( (prof_mem) and (x) ):
        plt.figure(1)
        plt.plot(x, y_m, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=ts2)

      if( (prof_energy) and (x) ):
        plt.figure(2)
        plt.plot(x, y_e, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=ts2)
      idx=idx+1

#  # add limits
#  sus_mem_bw = 40.0 #IB
#  # divide by streams number
#  if stencil_kernel == 0:
#    mem_limit = sus_mem_bw/4.0
#  elif stencil_kernel == 1:
#    mem_limit = sus_mem_bw/3.0
#  elif stencil_kernel == 2:
#    mem_limit = sus_mem_bw/5.0
#  elif stencil_kernel == 3:
#    mem_limit = sus_mem_bw/6.0
#  elif stencil_kernel == 4:
#    mem_limit = sus_mem_bw/16.0
#  elif stencil_kernel == 5:
#    mem_limit = sus_mem_bw/10.0
#  # divide by word size
#  if is_dp == 1: 
#    mem_limit = mem_limit/8.0
#  else:
#    mem_limit = mem_limit/4.0
#  if(prof_mem):
#    plt.plot([1, max_x], [mem_limit, mem_limit], color='.2', linestyle='--', label='Spt.lim.')

#  if mwdt==0 or mwdt==-1:
#    mwd_str = 'block'
#  elif mwdt==1:
#    mwd_str = 'fe'
#  elif mwdt==2:
#    mwd_str = 'rs'
#  elif mwdt==3:
#    mwd_str = 'fers'

  title = '_inc_grid_size'
  if stencil_kernel == 0:
    title = '25_pt_const' + title
  elif stencil_kernel == 1:
    title = '7_pt_const' + title
  elif stencil_kernel == 4:
    title = '25_pt_var' + title
  elif stencil_kernel == 5:
    title = '7_pt_var' + title
  elif stencil_kernel == 6:
    title = 'solar' + title


  plt.figure(0)
  plt.ylabel('GLUP/s')
  plt.grid()
  plt.ylim([0,max_y*1.1])
  f_name = title
  plt.xlabel('Size in each dimension')
  plt.legend(loc='best')
  pylab.savefig('perf_'+f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
  plt.clf()

  if prof_mem:
    plt.figure(1)
    plt.ylabel('GB/s')
    plt.grid()
    f_name = title
    plt.xlabel('Size in each dimension')
    plt.legend(loc='best')
    pylab.savefig('bw_'+f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    plt.clf()

  if prof_energy:
    plt.figure(2)
    plt.ylim([0, 150])
    plt.ylabel('pJ/LUP')
    plt.grid()
    f_name = title
    plt.xlabel('Size in each dimension')
    plt.legend(loc='best')
    pylab.savefig('energy_'+f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    plt.clf()


 
if __name__ == "__main__":
  main()
