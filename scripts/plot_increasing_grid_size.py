#!/usr/bin/env python
def main():
  import sys
  from scripts.utils import get_stencil_num, load_csv
  raw_data = load_csv(sys.argv[1])

  k_l = set()
  mwdt_l = set()
  use_prof = True 
  for k in raw_data:
    k_l.add(get_stencil_num(k))
    # add mwd type
    mwd = k['Wavefront parallel strategy'].lower()
    k['mwdt']=-1
    if('fixed' in mwd) and ('relaxed' in mwd):
      k['mwdt'] = 3
    elif('fixed' in mwd):
      k['mwdt'] = 1
    elif('relaxed' in mwd):
      k['mwdt'] = 2
    if(k['mwdt']>=0): mwdt_l.add(k['mwdt'])
    if 'Sustained Memory BW' not in k.keys(): use_prof=False
  k_l = list(k_l)
  mwdt_l = list(mwdt_l)


  for k in k_l:
    for mwdt in mwdt_l:
      plot_lines(raw_data, k, mwdt, use_prof=use_prof)

        
def plot_lines(raw_data, stencil_kernel, mwdt, use_prof=True, is_dp=1):
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
         'text.fontsize': 7*m,
         'legend.fontsize': 5*m,
         'xtick.labelsize': 6*m,
         'ytick.labelsize': 6*m,
         'lines.markersize': 1,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)

  

  req_fields = [('Local NX', int), ('Thread group size', int), ('Time stepper orig name', str), ('MStencil/s  MAX', float), ('mwdt', int)]
  if(use_prof): req_fields = req_fields + [ ('Sustained Memory BW', float)]
  data = []
  for k in raw_data:
    tup = {}
    # add the general fileds
    for f in req_fields:
      try:
        tup[f[0]] = map(f[1], [k[f[0]]] )[0]
      except:
        print("ERROR: results entry missing essential data")
        print f[0]
        print k
        return

    # add the stencil operator
    tup['stencil'] = get_stencil_num(k)

    # add the precision information
    if k['Precision'] in 'DP':
      p = 1
    else:
      p = 0
    tup['Precision'] = p
    data.append(tup)

  data = sorted(data, key=itemgetter('mwdt', 'Thread group size', 'Local NX'))
  #for i in data: print i


  tgs_l = set()
  for k in data:
    tgs_l.add(k['Thread group size'])
  tgs_l = list(tgs_l)

  max_x = 0
  max_y = 0

  marker_s = 7
  line_w = 1
  line_s = '-' 
  cols = 'kbgmrcy'
  markers = '+xo^vx*'
  for idx, tgs in enumerate(tgs_l):
    marker = markers[idx]
    col = cols[idx]
    x = []
    y = []
    y_m = []
    for k in data:
      if ( (k['mwdt']==-1 or  k['mwdt'] == mwdt) and (k['Thread group size'] == tgs)  and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp) ):
        if(use_prof): y_m.append(k['Sustained Memory BW']/10**3)
        x.append(k['Local NX'])
        y.append(k['MStencil/s  MAX']/10**3)
    ts2 = str(tgs) + 'WD' if tgs!=0 else 'Spt.blk.'

    if(x):
      plt.figure(0)
      plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=ts2)
      max_y = max(max(y), max_y)
      max_x = max(max(x), max_x)
    if( (use_prof) and (x) ):
      plt.figure(1)
      plt.plot(x, y_m, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=ts2)

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
#  if(use_prof):
#    plt.plot([1, max_x], [mem_limit, mem_limit], color='.2', linestyle='--', label='Spt.lim.')


  if mwdt==1:
    mwd_str = 'fe'
  elif mwdt==2:
    mwd_str = 'rs'
  elif mwdt==3:
    mwd_str = 'fers'

  title = '_inc_grid_size_mwdt_'+mwd_str
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

  if use_prof:
    plt.figure(1)
    plt.ylabel('GB/s')
    plt.grid()
    f_name = title
    plt.xlabel('Size in each dimension')
    plt.legend(loc='best')
    pylab.savefig('bw_'+f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    plt.clf()

 
if __name__ == "__main__":
  main()
