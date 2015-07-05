#!/usr/bin/env python

method_style = {'Spt.blk.':('g','o'), 'MWD':('k','x'), 'CATS2':('r','+'),
                'PLUTO':('m','*'), 'Pochoir':('b','^')}


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


  req_fields = [('MStencil/s  MAX', float), ('Precision', int), ('Global NX', int), ('Number of time steps', int), ('Number of tests', int)]

  hw_ctr_fields = {
                    '':[],
                    'TLB':[('L1 DTLB miss rate sum', float)],
                    'DATA':[('Load to Store ratio avg', float)],
                    'L2':[('L2 data volume sum', float)],
                    'L3':[('L3 data volume sum', float)],
                    'MEM':[('Total Memory Transfer', float),('Sustained Memory BW', float)],
                    'ENERGY':[('Energy', float), ('Energy DRAM', float), ('Power',float), ('Power DRAM', float)]}


  duplicates = set()
  meas_figs = dict()
  perf_fig = dict()
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

    #find repeated data
    key = (entry['Precision'], k['stencil_name'], k['LIKWID performance counter'], k['mwdt'], k['method'], entry['Global NX'])
    if key not in duplicates:
      duplicates.add(key)
    else:
      print("Repeated result at: %s"%(k['file_name']))
      continue


    # Initialize plot entry if does not exist for current data entry
#    for m,n in entry.iteritems(): print m,n
    measure_list = ['n', 'perf', 'total energy', 'tlb', 'mem bw', 'l2 bw', 'l3 bw', 'mem vol', 'l2 vol', 'l3 vol', 'data', 'tgs', 'thx', 'thy', 'thz', 'blk size', 'diam width', 'bs_z']
    plot_key = (entry['Precision'], k['stencil_name'], k['LIKWID performance counter'])
    line_key = (k['mwdt'], k['method'])
    if plot_key not in meas_figs.keys():
      meas_figs[plot_key] = {}
    if line_key not in meas_figs[plot_key].keys():
      meas_figs[plot_key][line_key] = {meas:[] for meas in measure_list}

    # append the measurement data
    meas_figs[plot_key][line_key]['n'].append(entry['Global NX'])
#    meas_figs[plot_key][line_key]['perf'].append(entry['MStencil/s  MAX']/1e3)
    N = entry['Global NX']**3 * entry['Number of time steps'] * entry['Number of tests']/1e9
    # Memory
    if k['LIKWID performance counter'] == 'MEM':
      meas_figs[plot_key][line_key]['mem bw'].append(entry['Sustained Memory BW']/1e3)
      meas_figs[plot_key][line_key]['mem vol'].append(entry['Total Memory Transfer']/N)
    # Energy
    elif k['LIKWID performance counter'] == 'ENERGY':
      entry['cpu energy pj/lup'] = entry['Energy']/N
      entry['dram energy pj/lup'] = entry['Energy DRAM']/N
      entry['total energy pj/lup'] = entry['cpu energy pj/lup'] + entry['dram energy pj/lup']
      if (entry['total energy pj/lup'] < 3e3):
#        entry['total energy pj/lup'] = 0
        meas_figs[plot_key][line_key]['total energy'].append(entry['total energy pj/lup'])
      else:
        del meas_figs[plot_key][line_key]['n'][-1]
    # TLB
    elif k['LIKWID performance counter'] == 'TLB':
      meas_figs[plot_key][line_key]['tlb'].append(entry[ hw_ctr_fields['TLB'][0][0] ])
    # L2
    elif k['LIKWID performance counter'] == 'L2':
      meas_figs[plot_key][line_key]['l2 vol'].append(entry['L2 data volume sum']/N)
    #L3
    elif k['LIKWID performance counter'] == 'L3':
      meas_figs[plot_key][line_key]['l3 vol'].append(entry['L3 data volume sum']/N)
    #CPU
    elif k['LIKWID performance counter'] == 'DATA':
      meas_figs[plot_key][line_key]['data'].append(entry['Load to Store ratio avg'])
    #Diamond tiling data
    if(k['method'] == 'CATS2' or k['method'] == 'MWD'):
      meas_figs[plot_key][line_key]['diam width'].append(int(k['Intra-diamond width']))
      meas_figs[plot_key][line_key]['tgs'].append(int(k['Thread group size']))
      meas_figs[plot_key][line_key]['thx'].append(int(k['Threads along x-axis']))
      meas_figs[plot_key][line_key]['thy'].append(int(k['Threads along y-axis']))
      meas_figs[plot_key][line_key]['thz'].append(int(k['Threads along z-axis']))
      meas_figs[plot_key][line_key]['blk size'].append(int(k['Total cache block size (kiB)'])/1024.0)
      meas_figs[plot_key][line_key]['bs_z'].append(int(k['Multi-wavefront updates']))

    # append the performance data
    plot_key = (entry['Precision'], k['stencil_name'])
    line_key = (k['mwdt'], k['method'])
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
  for p in meas_figs:
    for l in meas_figs[p]:
      pl = meas_figs[p][l]
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

#  for m,n in meas_figs.iteritems():
#    print "##############",m
#    for i,j in n.iteritems():
#      print i,j

  plot_all(perf_fig, meas_figs)


def plot_setup():
  import matplotlib.pyplot as plt
  import pylab

  # figure size and params
  fig_width = 20.0*0.393701 # inches
  fig_height = 1.0*fig_width #* 210.0/280.0#433.62/578.16
  fig_size =  [fig_width,fig_height]
  params = {
         'lines.linewidth': 1,
         'lines.markersize': 3,
         'lines.linestyle': '-',
         'axes.labelsize': 7,
         'axes.linewidth': 0.25,
         'font.size': 7,
         'legend.fontsize': 7,
         'xtick.labelsize': 6,
         'ytick.labelsize': 6,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)


  #create figures instance
  n_plt = 12
  plt_rows = 4
  plt_cols = 3
  fig, axarr = plt.subplots(plt_rows, plt_cols)

  # set figures information
  plt_info = {
    'perf':      (0,0, 'GLUP/s'),
    'mem bw':    (0,1, hw_ctr_labels['MEM'][0][0]),
    'tlb':       (0,2, hw_ctr_labels['TLB'][0][0]),

    'mem vol':   (1,0, hw_ctr_labels['MEM'][1][0]),
    'l3 vol':    (1,1, hw_ctr_labels['L3'][0][0]),
    'l2 vol':    (1,2, hw_ctr_labels['L2'][0][0]),

    'diam width':(2,0, 'Diamond width'),
    'bs_z':      (2,1, 'Block size along z-axis'),
    'blk size':  (2,2, 'Cache block size (MiB)'),

    'tgs':       (3,0, 'MWD Intra-tile threads'),
    'data':      (3,1, hw_ctr_labels['DATA'][0][0]),
  'total energy':(3,2, hw_ctr_labels['ENERGY'][0][0])}

  axl = dict()
  for (meas, (r,c, y_label)) in plt_info.iteritems():
    axl[meas] = axarr[r,c]

    # set the default parameters
    ax = axl[meas]
    ax.grid()
    ax.set_ylabel(y_label, labelpad=1)
    if(r==plt_rows-1): ax.set_xlabel('Size in each dimension')
    if(meas in ['tlb', 'total energy']):
      ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

  return axl, fig


def plot_finalize(stencil, axl, fig):
  import matplotlib.pyplot as plt
  import pylab

  # set figures limits
  for meas in axl:
    if meas in ['blk size', 'diam width', 'bs_z', 'tgs']:
      axl[meas].set_ylim(bottom=0)

  title = stencil.replace('_','-')
  fig.get_axes()[0].annotate(title, (0.5, 0.94), xycoords='figure fraction', ha='center', fontsize=12)
  pylab.savefig('all_'+ stencil + '_perf_inc_grid_size' + '.pdf', format='pdf', bbox_inches="tight", pad_inches=0.1)
  plt.clf()


def plot_all(perf_fig, meas_figs):
  stencils = [p[1] for p in perf_fig]

  for stencil in stencils:
    print(stencil)
    axl, fig = plot_setup()
    # Plot performance
    for p in perf_fig:
      if p[1] == stencil:
        plot_perf_fig(perf_fig[p], stencil, axl)

    # Plot hardware counters measurements
    for p in meas_figs:
      if p[1] == stencil:
        plot_meas_fig(meas_figs[p], stencil=p[1], plt_key=p[2], axl=axl)

    # Plot used parameters
    for p in meas_figs:
      if (p[1]==stencil and p[2]!=''):
        plot_params_fig(meas_figs[p], stencil=p[1], plt_key=p[2], axl=axl)
        break # because it will be the same for all HW counter data

    plot_finalize(stencil, axl, fig)


def plot_perf_fig(p, stencil, axl):
  import matplotlib.pyplot as plt
  import itertools
  from scripts.utils import get_stencil_num


  # performance
  for l in p:
    label = l[1]
    col, marker = method_style[label]
    x = []
    y = []
    for xp, y_l in p[l].iteritems():
      for x1, y1 in (itertools.product([xp], y_l)):
        x.append(x1)
        y.append(y1)
    axl['perf'].plot(x, y, color=col, marker=marker, label=label)


def plot_meas_fig(p, stencil, plt_key, axl):
  import matplotlib.pyplot as plt
  from scripts.utils import get_stencil_num

  # HW measurements
  for y_label, file_prefix, measure in hw_ctr_labels[plt_key]:
    for l in p:
      label = l[1]
      col, marker = method_style[label]
      x = p[l]['n']
      y = p[l][measure]
      axl[measure].plot(x, y, color=col, marker=marker, label=label)

    if(measure == 'tlb'):
      handles, labels = axl['perf'].get_legend_handles_labels()
      axl[measure].legend(handles, labels, loc='best')


def plot_params_fig(p, stencil, plt_key, axl):
  import matplotlib.pyplot as plt
  from scripts.utils import get_stencil_num


  # Diamond tiling information
  # Thread group size information
  tgs_labels =(
         ('tgs', 'b', '^', 'TGS'),
         ('thx', 'r', '+', 'x'),
         ('thy', 'g', 'o', 'y'),
         ('thz', 'm', '*', 'z') )
  for l in p:
    method=l[1]
    if(method == 'MWD'):
      for measure, col, marker, label in tgs_labels:
        x = p[l]['n']
        y = p[l][measure]
        axl['tgs'].plot(x, y, color=col, marker=marker, label=label)

      axl['tgs'].legend(loc='best')

  #Cache block size and diamond width
  blk_width_labels = [('blk size', 'Cache block size (MiB)', 'cache_block_size_'),
                      ('diam width', 'Diamond width', 'diamond_width_'),
                      ('bs_z', 'Block size along z-axis', 'bs_z_')]
  for measure, y_label, f_prefix in blk_width_labels:
    for l in p:
      method=l[1]
      if(method in ['MWD', 'CATS2']):
        col, marker = method_style[method]
        x = p[l]['n']
        y = p[l][measure]
        axl[measure].plot(x, y, color=col, marker=marker, label=method)


if __name__ == "__main__":
  main()
