#!/usr/bin/env python
def main():
  import sys

  raw_data = load_csv(sys.argv[1])

  k_l = set()
  for k in raw_data:
    k_l.add(get_stencil_num(k))
  k_l = list(k_l)

  for k in k_l:
    plot_stacked_clustered_bars(raw_data, k)
  #plot_lines(raw_data, plot_title)


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


def plot_stacked_clustered_bars(raw_data, kernel):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import pylab
  import numpy as np
  from pylab import arange,pi,sin,cos,sqrt

  sec_fontsize = 12

  fig_width = 2*7.4*0.393701 # inches
  fig_height = 0.75*fig_width #* 210.0/280.0#433.62/578.16
  fig_size =  [fig_width,fig_height]
  params = {
         'axes.labelsize': 7,
#         'axes.linewidth': 0.5,
#         'lines.linewidth': 0.75,
         'text.fontsize': 7,
         'legend.fontsize': 7,
         'xtick.labelsize': 7,
         'ytick.labelsize': 7,
#         'lines.markersize': 3,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)


  cluster_space = 2 
  
  req_fields = [('OpenMP Threads', int), ('Thread group size', int), ('WD main-loop RANK0 MStencil/s  MAX', float), ('Time stepper orig name', str), ('MPI size', int), ('MStencil/s  MAX', float)]

  data = []
  for k in raw_data:
    tup = {}
    # add the general fileds
    for f in req_fields:
      tup[f[0]] = map(f[1], [k[f[0]]] )[0]

    # add the stencil operator
    tup['stencil'] = get_stencil_num(k)

    # add the precision information
    if k['Precision'] in 'DP':
      p = 1
    else:
      p = 0
    tup['Precision'] = p

# add MWD profiling
    measures = ['Wavefront barrier wait [s]:', 'Wavefront barrier wait [%]:',
               'Wavefront steady state [s]:', 'Wavefront steady state [%]:',
               'Wavefront startup/end [s]:', 'Wavefront startup/end [%]:',
               'Wavefront communication [s]:', 'Wavefront communication [%]:',
               'Wavefront others [s]:', 'Wavefront others [%]:',
               'Group spin-wait [s]:', 'Group spinn-wait [%]:',
               'Resolved diamonds:']
    if tup['Time stepper orig name'] == 'Dynamic-Intra-Diamond':
      for i in range(tup['OpenMP Threads']):
        for m in measures:
          key =m[:-1] + ' group %d'%i
          if key in k:
            if k[key] != '':
              tup[key] = float(k[key]) if '.' in k[key] else int(k[key])

    measures =['RANK0 Total',
                'RANK0 Computation',
                'RANK0 Communication',
                'RANK0 Waiting',
                'RANK0 Other']
    if 'Naive' in tup['Time stepper orig name']:
      for m in measures:
        tup[m] = float(k[m])

    if tup['MPI size'] not in [12, 24]: # exclusion list
      data.append(tup)


  d2 = data
  data =[]
  for k in d2: # clean up the data and create the bars information
    k['ntg'] = 0
    k['comp_l'] = [0]
    k['comm_l'] = [0]
    k['idle_l'] = [0]
    k['others_l'] = [0]
    if k['Time stepper orig name'] == 'Dynamic-Intra-Diamond':
      k['ts'] = str(k['Thread group size']) + 'WD'
      k['ntg'] = k['OpenMP Threads']/k['Thread group size']
      comp_l = []
      comm_l = []
      idle_l = []
      others_l = []
      for i in range(k['ntg']):
        comp_l.append(k['Wavefront steady state [%] group '+str(i)]+
                 k['Wavefront startup/end [%] group '+str(i)]) 
        comm_l.append(k['Wavefront communication [%] group '+str(i)])
        idle_l.append(k['Group spinn-wait [%] group '+str(i)])
        others_l.append(100 - comp_l[-1] - comm_l[-1] - idle_l[-1])
      k['comp_l']= comp_l
      k['comm_l']= comm_l
      k['idle_l']= idle_l
      k['others_l']= others_l
      k['Performance'] = k['WD main-loop RANK0 MStencil/s  MAX'] 

    else:
      k['ts'] = 'Spt.blk.'
      k['comp_l'] = [k['RANK0 Computation']/k['RANK0 Total'] * 100.0]
      k['comm_l'] = [k['RANK0 Communication']/k['RANK0 Total'] * 100.0]
      k['Performance'] = k['MStencil/s  MAX'] 

#    del k['MStencil/s  MAX'] 
#    del k['WD main-loop RANK0 MStencil/s  MAX']
#    del k['Time stepper orig name']

    k['comp_std'] = np.std(k['comp_l'])
    k['comm_std'] = np.std(k['comm_l'])
    k['idle_std'] = np.std(k['idle_l'])
    k['others_std'] = np.std(k['others_l'])
    
    # Remove errors under certain threashold
    threash = 3
    if k['comp_std']<threash: k['comp_std'] = 0
    if k['comm_std']<threash: k['comm_std'] = 0
    if k['idle_std']<threash: k['idle_std'] = 0
    if k['others_std']<threash: k['others_std'] = 0

    k['comp_mean'] = np.mean(k['comp_l'])
    k['comm_mean'] = np.mean(k['comm_l'])
    k['idle_mean'] = np.mean(k['idle_l'])
    k['others_mean'] = np.mean(k['others_l'])

    k['comm_btm']   = np.mean(k['comp_mean'])
    k['idle_btm']   = np.mean(k['comm_btm']+k['comm_mean'])
    k['others_btm'] = np.mean(k['idle_btm']+k['idle_mean'])

    data.append(k)

  data = [k for k in data if k['stencil']==kernel]
  data = sorted(data, key=itemgetter('MPI size', 'Thread group size'))
  #for i in data: print i
#  for i in data: 
#    print "#####", i['MPI size'], i['ntg']
#    print i['comp_l'], i['comm_l'], i['idle_l'], i['others_l']

  np_l = [k['MPI size'] for k in data]
  ts_l = [k['ts'] for k in data]
  tsn_l = [k['Thread group size'] for k in data]
  perf_l = [k['Performance'] for k in data]
  nnp = len(set(np_l))
  nts = len(set(tsn_l))
  np_set = sorted(list(set(np_l)))
  tsn_set = sorted(list(set(tsn_l)))
  ts_set = []
  for i in tsn_set: 
    if i == 0: ts_set.append('Spt.blk.')
    else: ts_set.append(str(i)+'WD')
  np_i = dict(zip(np_set, range(nnp)))
  ts_i = dict(zip(ts_set, range(nts)))
  cluster_size = nts+cluster_space
  nx = cluster_size*nnp - cluster_space  
  x = [0.0]*nx

  fig = plt.figure()
  fig.subplots_adjust(bottom=0.25)  # space on the bottom for second time axis
  host = fig.add_subplot(111)  # setup plot

 
  p1=[]; p2=[]; p3=[]; p4=[]
  width =0.8
  ecolor = 'k'
  xtk = np.arange(nx)
  for i, k in enumerate(data):
    idx = np_i[k['MPI size']]*cluster_size + ts_i[k['ts']]
    p1 = host.bar(idx, k['comp_mean'], width, color='0.95', align='center',yerr=[(k['comp_std'],),[0]], hatch="/", error_kw={'elinewidth':2, 'ecolor':ecolor,'capsize':3})
    p2 = host.bar(idx, k['comm_mean'], width, color='0.65', align='center',bottom=k['comm_btm'], yerr=[(k['comm_std'],),[0]], hatch="", error_kw={'elinewidth':2, 'ecolor':ecolor,'capsize':3})
    p3 = host.bar(idx, k['idle_mean'], width, color='0.45', align='center',bottom=k['idle_btm'], yerr=[(k['idle_std'],),[0]], hatch="\\", error_kw={'elinewidth':2, 'ecolor':ecolor,'capsize':3})
    #p4 = host.bar(idx, k['others_mean'], width, color='r', align='center',bottom=k['others_btm'], yerr=[(k['others_std'],),[0]], hatch="//")
    p4 = host.bar(idx, k['others_mean'], width, color='0.0', align='center',bottom=k['others_btm'], hatch="")

  for i in range(nnp):
    xs = cluster_size*i + cluster_size-1.5
    host.plot([xs, xs], [0, 105], color='k', linestyle='--')

  for i in range(nx):
    if i%cluster_size >= nts:
      xtk[i]=0
  host.set_xticks(xtk)

  np_tickl = ['']*nx
  for i in range(nnp):
    idx = int(cluster_size*(i+0.5) - cluster_space/2)
    np_tickl[idx] = np_set[i]
  host.set_xticklabels(np_tickl)


  host.set_ylabel('Time (\%)', fontsize=sec_fontsize)
  host.set_xlabel('Processors \#', fontsize=sec_fontsize)
  host.tick_params(axis='both', which='major', labelsize=sec_fontsize)
  host.tick_params(axis='both', which='minor', labelsize=sec_fontsize)

  #plt.title(plot_title)
  #plt.grid() 

  # Insert the time steppers names at the X-axis 
  newax = host.twiny()  # create new axis
  newax.xaxis.set_ticks_position('bottom')
  newax.spines['bottom'].set_position(('outward', 15))

#  newax.patch.set_visible(False)
#  newax.xaxis.set_label_position('bottom')
#  newax.set_frame_on(False)
#  newax.tick_params('x', width=0)

  ticks = np.arange(nx)
  for i in range(nx):
    if (i >= nts) and (i < nx-nts):
      ticks[i] = 0
  newax.set_xticks(ticks)

  ts_tickl = ['']*nx
  for i in range(nts):
    ts_tickl[i] = ts_l[i]
    ts_tickl[i+nx-nts] = ts_l[i]
  newax.set_xticklabels(ts_tickl, rotation=90, size='small')

  newax.axis((0.0, float(nx), 0.0, 105))

  newax.set_xlim((-width, nx))
  host.set_xlim((-width, nx))


  # set color according to stencil type
  #my_colors = 'rgbkmcy'
  #ts_col = dict(zip(ts_set,my_colors[:nts]))
  #for i, k in enumerate(data):
  #  idx = np_i[k['MPI size']]*cluster_size + ts_i[k['ts']]
  #  p1[idx].set_color(ts_col[k['ts']])
 
  #pylab.legend(p1[:nts], ts_set)
  host.legend( (p4, p3, p2, p1), ('Others', 'Idle','Communicate', 'Compute'), loc='center left', fontsize=sec_fontsize)

  title = ''
#  title = 'Distributed memory performance'
#  plt.title(title)
  if kernel == 0:
    title = title + ' 25_pt constant coeff. star stencil'
  elif kernel == 1:
    title = title + ' 7_pt constant coeff. star stencil'
  elif kernel == 2:
    title = title + ' 7_pt varialbe. coeff. star stencil'
  elif kernel == 3:
    title = title + ' 7_pt varialbe. coeff. axis symm star stencil'
  elif kernel == 4:
    title = '25_pt_var_dist_comp_ratio'
    #title = title + ' 25_pt varialbe. coeff. axis symm star stencil'
  elif kernel == 5:
    title = '7_pt_var_dist_comp_ratio'
    #title = title + ' 7_pt varialbe. coeff. no symm star stencil'

  f_name = title.replace(' ', '_')
  pylab.savefig(f_name+'.png', bbox_inches="tight", pad_inches=0.04)
  pylab.savefig(f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)

  return
      
def load_csv(data_file):
  from csv import DictReader
  with open(data_file, 'rb') as output_file:
    data = DictReader(output_file)
    data = [k for k in data]
  return data
  
  
if __name__ == "__main__":
  main()
