#!/usr/bin/env python
marker_s = 7
line_w = 1
line_s = '-' 
method_style = {'MWD':('k','x'), '1WD':('r','+'),
                'Spt.blk.':('g','o'), 'PLUTO':('m','*'), 'Pochoir':('b','^'),
                     '2WD':('g','o'),   '3WD':('m','*'),     '6WD':('b','^'), '9WD':('c','v'), '18WD':('y','>'), '5WD':('m','*'), '10WD':('b','^')}


hw_ctr_labels = {
                    '':(),
                    'TLB':[('L1 DTLB miss rate sum', 'tlb_', 'tlb')],
                    'DATA':[('Load to Store ratio avg', 'cpu_', 'data')],
                    'L2':[('L2 Bytes/LUP', 'L2_', 'l2 vol')],
                    'L3':[('L3 Bytes/LUP', 'L3_', 'l3 vol')],
                    'MEM':[('MEM GB/s', 'mem_bw_', 'mem bw'), ('MEM Bytes/LUP', 'mem_vol_', 'mem vol')],
                    'ENERGY':[('CPU pJ/LUP', 'energy_cpu_', 'cpu energy'),
                              ('DRAM pJ/LUP', 'energy_dram_', 'dram energy'),
                              ('Total pJ/LUP', 'energy_total_', 'total energy')] }

req_fields = [('MStencil/s  MAX', float), ('Precision', int), ('Global NX', int), ('Number of time steps', int), ('Number of tests', int), ('Thread group size', int)]

hw_ctr_fields = {
                    '':[],
                    'TLB':[('L1 DTLB miss rate sum', float)],
                    'DATA':[('Load to Store ratio avg', float)],
                    'L2':[('L2 data volume sum', float)],
                    'L3':[('L3 data volume sum', float)],
                    'MEM':[('Total Memory Transfer', float),('Sustained Memory BW', float)],
                    'ENERGY':[('Energy', float), ('Energy DRAM', float), ('Power',float), ('Power DRAM', float)]}


measure_list = ['n', 'perf', 'cpu energy', 'dram energy', 'total energy', 'tlb', 'mem bw', 'l2 bw', 'l3 bw', 'mem vol', 'l2 vol', 'l3 vol', 'data', 'tgs', 'thc', 'thx', 'thy', 'thz', 'blk size', 'diam width', 'wavefront width', 'pluto bs_x']

def init_figs():
  import pylab

  fig_scale = 3.0
  fig_width = 4.0*0.393701*fig_scale # inches
  fig_height = 1.0*fig_width #* 210.0/280.0#433.62/578.16

  fig_size =  [fig_width,fig_height]
  params = {
         'axes.labelsize': 6*fig_scale,
         'axes.linewidth': 0.25*fig_scale,
         'lines.linewidth': 0.75*fig_scale,
         'font.size': 7*fig_scale,
         'legend.fontsize': 4*fig_scale,
         'xtick.labelsize': 5*fig_scale,
         'ytick.labelsize': 6*fig_scale,
         'lines.markersize': 1,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)



def parse_entry_info(k, is_tgs_only):
  from scripts.utils import get_stencil_num, load_csv

#    # get processor name from the file names
#    if(k['OpenMP Threads']!=''):
#      if(int(k['OpenMP Threads']) == 10):
#        machine_name = 'ivb10'
#      elif(int(k['OpenMP Threads']) == 18):
#        machine_name = 'hw18'
#

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

  # parse the general fileds' format
  for f in req_fields + hw_ctr_fields[k['LIKWID performance counter']]:
    try:
      k[f[0]] = map(f[1], [k[f[0]]] )[0]
    except:
      print("ERROR: results entry missing essential data at file:%s"%(k['file_name']))
      print f[0]
      print k
      return


  k['tgsl'] = k['Thread group size']
  if(is_tgs_only==0): # regular mode for all MWD
    if(k['method']=='MWD'):
      k['tgsl'] = 100


def init_plot_entry(plots, k):
  # Initialize plot entry if does not exist for current data entry
#    for m,n in k.iteritems(): print m,n
  plot_key = (k['Precision'], k['stencil_name'], k['LIKWID performance counter'])
  line_key = (k['mwdt'], k['method'], k['tgsl'])
  if plot_key not in plots.keys():
    plots[plot_key] = {}
  if line_key not in plots[plot_key].keys():
    plots[plot_key][line_key] = {meas:[] for meas in measure_list}


def append_meas_data(perf_fig, plots, k):
  # append the measurement data
  plot_key = (k['Precision'], k['stencil_name'], k['LIKWID performance counter'])
  line_key = (k['mwdt'], k['method'], k['tgsl'])

  plots[plot_key][line_key]['n'].append(k['Global NX'])
  plots[plot_key][line_key]['perf'].append(k['MStencil/s  MAX']/1e3)
  N = k['Global NX']**3 * k['Number of time steps'] * k['Number of tests']/1e9
  # Memory
  if k['LIKWID performance counter'] == 'MEM':
    plots[plot_key][line_key]['mem bw'].append(k['Sustained Memory BW']/1e3)
    plots[plot_key][line_key]['mem vol'].append(k['Total Memory Transfer']/N)
  # Energy
  elif k['LIKWID performance counter'] == 'ENERGY':
    k['cpu energy pj/lup'] = k['Energy']/N
    k['dram energy pj/lup'] = k['Energy DRAM']/N
    k['total energy pj/lup'] = k['cpu energy pj/lup'] + k['dram energy pj/lup']
    if (k['cpu energy pj/lup'] < 1e5): 
      plots[plot_key][line_key]['cpu energy'].append(k['cpu energy pj/lup'])
    if (k['dram energy pj/lup'] < 1e5): 
      plots[plot_key][line_key]['dram energy'].append(k['dram energy pj/lup'])
    if (k['total energy pj/lup'] < 1e5): 
      plots[plot_key][line_key]['total energy'].append(k['total energy pj/lup'])
  # TLB
  elif k['LIKWID performance counter'] == 'TLB':
    plots[plot_key][line_key]['tlb'].append(k[ hw_ctr_fields['TLB'][0][0] ])
  # L2
  elif k['LIKWID performance counter'] == 'L2':
    plots[plot_key][line_key]['l2 vol'].append(k['L2 data volume sum']/N)
  #L3
  elif k['LIKWID performance counter'] == 'L3':
    plots[plot_key][line_key]['l3 vol'].append(k['L3 data volume sum']/N)
  #CPU
  elif k['LIKWID performance counter'] == 'DATA':
    plots[plot_key][line_key]['data'].append(k['Load to Store ratio avg'])
  #Diamond tiling data
  if(k['method'] == '1WD' or k['method'] == 'MWD'):
    plots[plot_key][line_key]['diam width'].append(int(k['Intra-diamond width']))
    plots[plot_key][line_key]['wavefront width'].append(int(k['Multi-wavefront updates']))
    plots[plot_key][line_key]['tgs'].append(int(k['Thread group size']))
    plots[plot_key][line_key]['thc'].append(int(k['Threads per cell    ']))
    plots[plot_key][line_key]['thx'].append(int(k['Threads along x-axis']))
    plots[plot_key][line_key]['thy'].append(int(k['Threads along y-axis']))
    plots[plot_key][line_key]['thz'].append(int(k['Threads along z-axis']))
    plots[plot_key][line_key]['blk size'].append(int(k['Total cache block size (kiB)'])/1024.0)

  #PLUTO data
  if(k['method'] == 'PLUTO'):
    plots[plot_key][line_key]['diam width'].append(int(k['PLUTO tile size of loop 1']))
    plots[plot_key][line_key]['wavefront width'].append(int(k['PLUTO tile size of loop 3']))
    plots[plot_key][line_key]['pluto bs_x'].append(int(k['PLUTO tile size of loop 4']))

  # append the performance data
  plot_key = (k['Precision'], k['stencil_name'])
  line_key = (k['mwdt'], k['method'], k['tgsl'])
  if plot_key not in perf_fig.keys(): # figure
    perf_fig[plot_key] = dict()

  perf_line = perf_fig[plot_key]
  if line_key not in perf_line.keys(): # line
    perf_line[line_key] = dict()

  perf_point = perf_line[line_key]
  nx = k['Global NX'] 
  if nx not in perf_point.keys(): # points
    perf_point[nx] = [k['MStencil/s  MAX']/1e3]
  else:
    perf_point[nx].append(k['MStencil/s  MAX']/1e3)




def gen_plot_data(data_file, is_tgs_only):
  from scripts.utils import load_csv

  raw_data = load_csv(data_file)

  duplicates = set()
  plots = dict()
  perf_fig = dict()
  for k in raw_data:

    parse_entry_info(k, is_tgs_only)

    # Handle repeated data
    key = (k['Precision'], k['stencil_name'], k['LIKWID performance counter'], k['mwdt'], k['method'], k['tgsl'], k['Global NX'])
    if key not in duplicates:
      duplicates.add(key)
    else:
      print("Repeated result at: %s"%(k['file_name']))
      continue

    init_plot_entry(plots, k)
    append_meas_data(perf_fig, plots,k)

  del raw_data
  return plots, perf_fig


def sort_perf_fig(perf_fig):
  from collections import OrderedDict
  #sort performance results
  for k,v in perf_fig.iteritems():
    for k2,v2 in perf_fig[k].iteritems():
      perf_line = perf_fig[k][k2]
      perf_fig[k][k2] = OrderedDict(sorted(perf_fig[k][k2].iteritems(), key=lambda x:x[0]))
#  for k,v in perf_fig.iteritems():
#    print(k, "##########")
#    for k2,v2 in perf_fig[k].iteritems():
#      print(k2,v2)


def sort_meas_fig(plots):
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




