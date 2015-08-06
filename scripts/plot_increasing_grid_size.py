#!/usr/bin/env python

import pylab

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
       'xtick.labelsize': 5*m,
       'ytick.labelsize': 6*m,
       'lines.markersize': 1,
       'text.usetex': True,
       'figure.figsize': fig_size}
pylab.rcParams.update(params)


marker_s = 7
line_w = 1
line_s = '-' 
method_style = {'MWD':('k','x'), '1WD':('r','+'),
                'Spt.blk.':('g','o'), 'PLUTO':('m','*'), 'Pochoir':('b','^'),
                     '2WD':('g','o'),   '3WD':('m','*'),     '6WD':('b','^'), '9WD':('c','v'), '18WD':('y','>'), '5WD':('m','*')}


hw_ctr_labels = {
                    '':(),
                    'TLB':[('L1 DTLB miss rate sum', 'tlb_', 'tlb')],
                    'DATA':[('Load to Store ratio avg', 'cpu_', 'data')],
                    'L2':[('L2 Bytes/LUP', 'L2_', 'l2 vol')],
                    'L3':[('L3 Bytes/LUP', 'L3_', 'l3 vol')],
                    'MEM':[('MEM GB/s', 'mem_bw_', 'mem bw'), ('MEM Bytes/LUP', 'mem_vol_', 'mem vol')],
                    'ENERGY':[('Total pJ/LUP', 'energy_', 'total energy')] }
 
def main():
  import sys
  from scripts.utils import get_stencil_num, load_csv
  from collections import OrderedDict

  raw_data = load_csv(sys.argv[1])
  is_tgs_only = 0 if len(sys.argv)<3 else int(sys.argv[2])

  req_fields = [('MStencil/s  MAX', float), ('Precision', int), ('Global NX', int), ('Number of time steps', int), ('Number of tests', int), ('Thread group size', int)]

  hw_ctr_fields = {
                    '':[],
                    'TLB':[('L1 DTLB miss rate sum', float)],
                    'DATA':[('Load to Store ratio avg', float)],
                    'L2':[('L2 data volume sum', float)],
                    'L3':[('L3 data volume sum', float)],
                    'MEM':[('Total Memory Transfer', float),('Sustained Memory BW', float)],
                    'ENERGY':[('Energy', float), ('Energy DRAM', float), ('Power',float), ('Power DRAM', float)]}

 
  duplicates = set()
  plots = dict()
  perf_fig = dict()
  machine_name = ''
  for k in raw_data:

    # get processor name from the file names
    if(k['OpenMP Threads']!=''):
      if(int(k['OpenMP Threads']) == 10):
        machine_name = 'ivb10'
      elif(int(k['OpenMP Threads']) == 18):
        machine_name = 'hw18'


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


    # add the approach
    if(k['Time stepper orig name'] == 'Spatial Blocking'):
      k['method'] = 'Spt.blk.'
    elif(k['Time stepper orig name'] in ['PLUTO', 'Pochoir']):
      k['method'] = k['Time stepper orig name']
    elif(k['Time stepper orig name'] == 'Diamond'):
      if('_tgs1_' in k['file_name']):
        k['method'] = '1WD'
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


    # TLB measurement for LIKWID 4
    if 'L1 DTLB load miss rate avg' in k.keys():
      if k['L1 DTLB load miss rate avg']!='':
        hw_ctr_fields['TLB'] =  [('L1 DTLB load miss rate avg', float)]
        hw_ctr_labels['TLB'] =  [('L1 DTLB load miss rate avg', 'tlb_', 'tlb')]

    entry = {}
    # parse the general fileds' format
    for f in req_fields + hw_ctr_fields[k['LIKWID performance counter']]:
      try:
        entry[f[0]] = map(f[1], [k[f[0]]] )[0]
      except:
        print("ERROR: results entry missing essential data at file:%s"%(k['file_name']))
        print f[0]
        print k
        return


    k['tgsl'] = entry['Thread group size']
    if(is_tgs_only==0): # regular mode for all MWD
      if(k['method']=='MWD'):
        k['tgsl'] = 100


    #find repeated data
    key = (entry['Precision'], k['stencil_name'], k['LIKWID performance counter'], k['mwdt'], k['method'], k['tgsl'], entry['Global NX'])
    if key not in duplicates:
      duplicates.add(key)
    else:
      print("Repeated result at: %s"%(k['file_name']))
      continue


    # Initialize plot entry if does not exist for current data entry
#    for m,n in entry.iteritems(): print m,n
    measure_list = ['n', 'perf', 'total energy', 'tlb', 'mem bw', 'l2 bw', 'l3 bw', 'mem vol', 'l2 vol', 'l3 vol', 'data', 'tgs', 'thx', 'thy', 'thz', 'blk size', 'diam width', 'pluto z tile', 'pluto y tile', 'pluto x tile']
    plot_key = (entry['Precision'], k['stencil_name'], k['LIKWID performance counter'])
    line_key = (k['mwdt'], k['method'], k['tgsl'])
    if plot_key not in plots.keys():
      plots[plot_key] = {}
    if line_key not in plots[plot_key].keys():
      plots[plot_key][line_key] = {meas:[] for meas in measure_list}

    # append the measurement data
    plots[plot_key][line_key]['n'].append(entry['Global NX'])
    plots[plot_key][line_key]['perf'].append(entry['MStencil/s  MAX']/1e3)
    N = entry['Global NX']**3 * entry['Number of time steps'] * entry['Number of tests']/1e9
    # Memory
    if k['LIKWID performance counter'] == 'MEM':
      plots[plot_key][line_key]['mem bw'].append(entry['Sustained Memory BW']/1e3)
      plots[plot_key][line_key]['mem vol'].append(entry['Total Memory Transfer']/N)
    # Energy
    elif k['LIKWID performance counter'] == 'ENERGY':
      entry['cpu energy pj/lup'] = entry['Energy']/N
      entry['dram energy pj/lup'] = entry['Energy DRAM']/N
      entry['total energy pj/lup'] = entry['cpu energy pj/lup'] + entry['dram energy pj/lup']
      if (entry['total energy pj/lup'] > 1e5): entry['total energy pj/lup'] = 0
      plots[plot_key][line_key]['total energy'].append(entry['total energy pj/lup'])
    # TLB
    elif k['LIKWID performance counter'] == 'TLB':
      plots[plot_key][line_key]['tlb'].append(entry[ hw_ctr_fields['TLB'][0][0] ])
    # L2
    elif k['LIKWID performance counter'] == 'L2':
      plots[plot_key][line_key]['l2 vol'].append(entry['L2 data volume sum']/N)
    #L3
    elif k['LIKWID performance counter'] == 'L3':
      plots[plot_key][line_key]['l3 vol'].append(entry['L3 data volume sum']/N)
    #CPU
    elif k['LIKWID performance counter'] == 'DATA':
      plots[plot_key][line_key]['data'].append(entry['Load to Store ratio avg'])
    #Diamond tiling data
    if(k['method'] == '1WD' or k['method'] == 'MWD'):
      plots[plot_key][line_key]['diam width'].append(int(k['Intra-diamond width']))
      plots[plot_key][line_key]['tgs'].append(int(k['Thread group size']))
      plots[plot_key][line_key]['thx'].append(int(k['Threads along x-axis']))
      plots[plot_key][line_key]['thy'].append(int(k['Threads along y-axis']))
      plots[plot_key][line_key]['thz'].append(int(k['Threads along z-axis']))
      plots[plot_key][line_key]['blk size'].append(int(k['Total cache block size (kiB)'])/1024.0)

    #PLUTO data
    if(k['method'] == 'PLUTO'):
      plots[plot_key][line_key]['pluto z tile'].append(int(k['PLUTO tile size of loop 1']))
      plots[plot_key][line_key]['pluto y tile'].append(int(k['PLUTO tile size of loop 3']))
      plots[plot_key][line_key]['pluto x tile'].append(int(k['PLUTO tile size of loop 4']))

    # append the performance data
    plot_key = (entry['Precision'], k['stencil_name'])
    line_key = (k['mwdt'], k['method'], k['tgsl'])
    if plot_key not in perf_fig.keys(): # figure
      perf_fig[plot_key] = dict()

    perf_line = perf_fig[plot_key]
    if line_key not in perf_line.keys(): # line
      perf_line[line_key] = dict()

    perf_point = perf_line[line_key]
    nx = entry['Global NX'] 
    if nx not in perf_point.keys(): # points
      perf_point[nx] = [entry['MStencil/s  MAX']/1e3]
    else:
      perf_point[nx].append(entry['MStencil/s  MAX']/1e3)


  del raw_data

  #sort performance results
  for k,v in perf_fig.iteritems():
    for k2,v2 in perf_fig[k].iteritems():
      perf_line = perf_fig[k][k2]
      perf_fig[k][k2] = OrderedDict(sorted(perf_fig[k][k2].iteritems(), key=lambda x:x[0]))
#  for k,v in perf_fig.iteritems():
#    print(k, "##########")
#    for k2,v2 in perf_fig[k].iteritems():
#      print(k2,v2)


  #sort the plot lines
  for p in plots:
    for l in plots[p]:
      pl = plots[p][l]
      #remove unused fields
      empty = []
      for key, val in pl.iteritems():
        if(val==[]):
          empty.append(key)
      for key in empty:
          del pl[key]
      lines = []
      [lines.append(pl[val]) for val in measure_list if val in pl.keys()]

      lines = sorted(zip(*lines))
      idx=0
      for val in measure_list:
        if(val in pl.keys()):
          if(pl[val]):
            pl[val] = [x[idx] for x in lines]
            idx = idx+1

#  for m,n in plots.iteritems(): 
#    print "##############",m
#    for i,j in n.iteritems():
#      print i,j

  # Plot performance
  for p in perf_fig:
    plot_perf_fig(perf_fig[p], stencil=p[1], machine_name=machine_name, is_tgs_only=is_tgs_only)
 

  # Plot other measurements
  for p in plots:
    plot_meas_fig(plots[p], stencil=p[1], plt_key=p[2], machine_name=machine_name, is_tgs_only=is_tgs_only)


def plot_perf_fig(p, stencil, machine_name, is_tgs_only):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import matplotlib
  import pylab, itertools
  from pylab import arange,pi,sin,cos,sqrt
  from scripts.utils import get_stencil_num

  f_name = stencil+'_inc_grid_size'


  # performance
  plt.figure(0)
  for l in p:
    label = l[1]
    tgsl = l[2]
    if((label=='MWD') and (is_tgs_only==1)):
      label= str(tgsl) + 'WD'
    col, marker = method_style[label]
    x = []
    y = []
    for xp, y_l in p[l].iteritems():
      for x1, y1 in (itertools.product([xp], y_l)):
        x.append(x1)
        y.append(y1)
    plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=label)

  plt.ylabel('GLUP/s')
  plt.grid()
  plt.xlabel('Size in each dimension')
  plt.legend(loc='best')
  plt.gca().set_ylim(bottom=0)
  pylab.savefig(machine_name + '_perf_' + f_name + '.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
  plt.clf()


def plot_meas_fig(p, stencil, plt_key, machine_name, is_tgs_only):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import matplotlib
  import pylab
  from pylab import arange,pi,sin,cos,sqrt
  from scripts.utils import get_stencil_num

  f_name = stencil+'_inc_grid_size'


  # HW measurements
  plt_idx=1
  for y_label, file_prefix, measure in hw_ctr_labels[plt_key]:
    plt.figure(plt_idx)
    plt_idx = plt_idx + 1
    for l in p:
      label = l[1]
      tgsl = l[2]
      if((label=='MWD') and (is_tgs_only==1)):
        label= str(tgsl) + 'WD'
      col, marker = method_style[label]
      x = p[l]['n']
      y = p[l][measure]
      plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=label)

    plt.ylabel(y_label)
    plt.grid()
    plt.xlabel('Size in each dimension')
    plt.legend(loc='best')
    plt.gca().set_ylim(bottom=0)
    pylab.savefig(machine_name + '_' + file_prefix + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    plt.clf()


  # PLUTO tiling information
  if (any(method[1] == 'PLUTO' for method in p) and plt_key=='MEM'):
    plt.figure(plt_idx)
    plt_idx = plt_idx + 1
    for l in p:
      method=l[1]
      if(method == 'PLUTO'):
        tgs_labels =(
               ('pluto z tile', 'b', '^', 'Z'),
               ('pluto y tile', 'r', '+', 'Y'),
               ('pluto x tile', 'm', '*', 'X') )
        for measure, col, marker, label in tgs_labels:
          x = p[l]['n']
          y = p[l][measure]
          plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=label)

    plt.ylabel('PLUTO tile size')
    plt.grid()
    plt.xlabel('Size in each dimension')
    plt.legend(loc='upper left')
    plt.gca().set_ylim(bottom=0)
    pylab.savefig(machine_name + '_pluto_tile_size_' + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    plt.clf()


  # Diamond tiling information
  if (any(method[1] in ['MWD', '1WD'] for method in p) and plt_key=='MEM'):

    #Cache block size and diamond width
    for measure, y_label, f_prefix in [('blk size', 'Cache block size (MiB)', 'cache_block_size_'),
                                        ('diam width', 'Diamond width', 'diamond_width_')]:
      plt.figure(plt_idx)
      plt_idx = plt_idx + 1
      for l in p:
        method=l[1]
        tgsl = l[2]
        if(method in ['MWD', '1WD']):
          if(is_tgs_only==1):
            method= str(tgsl) + 'WD'
          col, marker = method_style[method]
          x = p[l]['n']
          y = p[l][measure]
          plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=method)

      plt.ylabel(y_label)
      plt.grid()
      plt.xlabel('Size in each dimension')
      plt.legend(loc='best')
      plt.gca().set_ylim(bottom=0)
      pylab.savefig(machine_name + '_' + f_prefix + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
      plt.clf()

    if(is_tgs_only==1): return
    # Thread group size information
    plt.figure(plt_idx)
    plt_idx = plt_idx + 1
    for l in p:
      method=l[1]
      if(method == 'MWD'):
        tgs_labels =(
               ('tgs', 'b', '^', 'MWD Group'),
               ('thx', 'r', '+', 'Along x'),
               ('thy', 'g', 'o', 'Along y'),
               ('thz', 'm', '*', 'Along z') )
        for measure, col, marker, label in tgs_labels:
          x = p[l]['n']
          y = p[l][measure]
          plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=label)

    plt.ylabel('Intra-tile threads')
    plt.grid()
    plt.xlabel('Size in each dimension')
    plt.legend(loc='upper left')
    plt.gca().set_ylim(bottom=0)
    pylab.savefig(machine_name + '_thread_group_size_' + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    plt.clf()


if __name__ == "__main__":
  main()
