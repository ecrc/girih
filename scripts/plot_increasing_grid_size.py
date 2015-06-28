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


  req_fields = [('method', str), ('MStencil/s  MAX', float), ('mwdt', str), ('Precision', int), ('stencil', int), ('Global NX', int), ('Number of time steps', int), ('Number of tests', int), ('Thread group size', int), ('LIKWID performance counter', str), ('stencil_name', str)]
  hw_ctr_fields = {
                    '':[],
                    'TLB':[('L1 DTLB miss rate sum', float)],
                    'DATA':[('Load to Store ratio avg', float)],
                    'L2':[('L2 data volume sum', float)],
                    'L3':[('L3 data volume sum', float)],
                    'MEM':[('Total Memory Transfer', float),('Sustained Memory BW', float)],
                    'ENERGY':[('Energy', float), ('Energy DRAM', float), ('Power',float), ('Power DRAM', float)]}
  hw_ctr_labels = {
                    '':(),
                    'TLB':('L1 DTLB miss rate sum', 'tlb_'),
                    'DATA':('Load to Store ratio avg', 'cpu_'),
                    'L2':('Bytes/LUP', 'L2_'),
                    'L3':('Bytes/LUP', 'L3_'),
                    'MEM':('GB/s', 'bw_'),
                    'ENERGY':('pJ/LUP', 'energy_') }
  

  plots = dict()
  for k in raw_data:

    # Use single field to represent the performance
    if 'Total RANK0 MStencil/s MAX' in k.keys():
      if(k['Total RANK0 MStencil/s MAX']!=''):
        k['MStencil/s  MAX'] = k['MWD main-loop RANK0 MStencil/s MAX'] 
    # temporary for deprecated format
    if 'RANK0 MStencil/s  MAX' in k.keys():
      if k['RANK0 MStencil/s  MAX']!='':
        k['MStencil/s  MAX'] = k['RANK0 MStencil/s  MAX'] 


    # add stencil operator
    k['stencil'] = get_stencil_num(k)
    k_l.add(k['stencil'])
    if   k['stencil'] == 0:
      k['stencil_name'] = '25_pt_const'
    elif k['stencil'] == 1:
      k['stencil_name'] = '7_pt_const'
    elif k['stencil'] == 4:
      k['stencil_name']  = '25_pt_var'
    elif k['stencil'] == 5:
      k['stencil_name']  = '7_pt_var'
    elif k['stencil'] == 6:
      k['stencil_name']  = 'solar'


    # disable bandwidth and energy plots if missing from any entry
    if 'Sustained Memory BW' in k.keys():
      if k['Sustained Memory BW'] =='':
        prof_mem=False
    else: 
      prof_mem=False
    if 'Energy' not in k.keys(): prof_energy=False

    # add the approach
    if(k['Time stepper orig name'] == 'Spatial Blocking'):
      k['method'] = 'Spt.blk.'
    elif(k['Time stepper orig name'] in ['PLUTO', 'Pochoir']):
      k['method'] = k['Time stepper orig name']
    elif(k['Time stepper orig name'] == 'Diamond'):
      if('_tgs1_' in k['file_name']):
        k['method'] = 'CATS2'
      else:
        k['method'] = 'MWD'
    else:
      print("ERROR: Unknow time stepper")
      raise

    # add mwd type
    k['mwdt']='none'
    if(k['method'] == 'MWD'):
      mwd = k['Wavefront parallel strategy'].lower()
      if('fixed' in mwd) and ('relaxed' in mwd):
        k['mwdt'] = 'fers'
      elif('fixed' in mwd):
        k['mwdt'] = 'fe'
      elif('relaxed' in mwd):
        k['mwdt'] = 'rs'
      elif('wavefront' in mwd):
        k['mwdt'] = 'block'


    # add precision information
    p = 1 if k['Precision'] in 'DP' else 0
    k['Precision'] = p
    prec_l.add(p)


    entry = {}
    # add the general fileds
    for f in req_fields + hw_ctr_fields[k['LIKWID performance counter']]:
      try:
        entry[f[0]] = map(f[1], [k[f[0]]] )[0]
      except:
        print("ERROR: results entry missing essential data at file:%s"%(k['file_name']))
        print f[0]
        print k
        return

#    for m,n in entry.iteritems(): print m,n
    plot_key = (entry['Precision'], entry['stencil_name'], entry['LIKWID performance counter'])
    line_key = (k['mwdt'], k['method'])
    if plot_key not in plots.keys():
      plots[plot_key] = {}
    if line_key not in plots[plot_key].keys():
      plots[plot_key][line_key] = [[], [], []]


    # append the data
    plots[plot_key][line_key][0].append(entry['Global NX'])
    plots[plot_key][line_key][1].append(entry['MStencil/s  MAX'])
    N = entry['Global NX']**3 * entry['Number of time steps'] * entry['Number of tests']/1e9
    # Memory
    if k['LIKWID performance counter'] == 'MEM':
      plots[plot_key][line_key][2].append(entry['Sustained Memory BW']/1e3)
    # Energy
    if k['LIKWID performance counter'] == 'ENERGY':
      entry['cpu energy pj/lup'] = entry['Energy']/N
      entry['dram energy pj/lup'] = entry['Energy DRAM']/N
      entry['total energy pj/lup'] = entry['cpu energy pj/lup'] + entry['dram energy pj/lup']
      if (entry['total energy pj/lup'] > 1e5): entry['total energy pj/lup'] = 0
      plots[plot_key][line_key][2].append(entry['total energy pj/lup'])
    # TLB
    if k['LIKWID performance counter'] == 'TLB':
      plots[plot_key][line_key][2].append(entry['L1 DTLB miss rate sum'])
    # L2
    if k['LIKWID performance counter'] == 'L2':
      plots[plot_key][line_key][2].append(entry['L2 data volume sum']/N)
    #L3
    if k['LIKWID performance counter'] == 'L3':
      plots[plot_key][line_key][2].append(entry['L3 data volume sum']/N)
    #CPU
    if k['LIKWID performance counter'] == 'DATA':
      plots[plot_key][line_key][2].append(entry['Load to Store ratio avg'])
 
  del raw_data
  #sort the plot lines
  for p in plots:
    for l in plots[p]:
      x = plots[p][l][0]
      y = plots[p][l][1]
      z = plots[p][l][2]
      xyz = sorted(zip(x,y,z))
      plots[p][l][0] = [x for (x,y,z) in xyz]
      plots[p][l][1] = [y for (x,y,z) in xyz]
      plots[p][l][2] = [z for (x,y,z) in xyz]

#  for m,n in plots.iteritems(): 
#    print "##############",m
#    for i,j in n.iteritems():
#      print i,j

  for p in plots:
    y_label, file_prefix = hw_ctr_labels[p[2]]
    stencil = p[1]
    plot_line(plots[p], stencil, y_label, file_prefix)



def plot_line(p, stencil, y_label, file_prefix):
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

  marker_s = 7
  line_w = 1
  line_s = '-' 
  cols = 'kbgmrcy'
  markers = '+xo^vx*'
  idx=0
  for idx, l in enumerate(p):
    label = l[1]
    marker = markers[idx]
    col = cols[idx]
    x = p[l][0]
    y_p = p[l][1]
    y = p[l][2]

    plt.figure(0)
    plt.plot(x, y_p, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=label)

    plt.figure(1)
    plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=label)


  f_name = stencil+'_inc_grid_size'

  plt.figure(0)
  plt.ylabel('GLUP/s')
  plt.grid()
  plt.xlabel('Size in each dimension')
  plt.legend(loc='best')
  pylab.savefig(file_prefix + 'perf_' + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
  plt.clf()


  plt.figure(1)
  plt.ylabel(y_label)
  plt.grid()
  plt.xlabel('Size in each dimension')
  plt.legend(loc='best')
  pylab.savefig(file_prefix + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
  plt.clf()


if __name__ == "__main__":
  main()
