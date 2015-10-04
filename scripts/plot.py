#!/usr/bin/env python
def main():
  import sys
  from scripts.plot_utils import measure_list
  import scripts.plot_utils as pu
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument('--data-file', type=str, default='summary.csv', help='Input data file name')
  parser.add_argument('--tgs', default=False, action='store_true', help='indicate TGS experiment')
  parser.add_argument('--thread-scaling', default=False, action='store_true', help='indicate thread scaling experiment')
  parser.add_argument('--machine-name', type=str, default='', help='Machine name to prefix results files')
  args = parser.parse_args()
  data_file = args.data_file
  is_tgs_only = int(args.tgs)
  is_thread_scaling = int(args.thread_scaling)
  machine_name = args.machine_name

  pu.init_figs()

  if (is_thread_scaling==0):
    data_x_label = 'Global NX'
    x_value = 'n'
    x_label = 'Size in each dimension'
    exp_name = 'inc_grid_size'
  else:
    data_x_label = 'OpenMP Threads' # 'Global NX'
    x_value = 'threads' # 'n'
    x_label = 'Threads number' # 'Size in each dimension'
    exp_name = 'thread_scaling' # 'inc_grid_size'


  plots, perf_fig = pu.gen_plot_data(data_file, is_tgs_only, x_label=data_x_label)

  # Plot performance
  for p in perf_fig:
    print p[1], ' Perf'
    plot_perf_fig(perf_fig[p], stencil=p[1], machine_name=machine_name, is_tgs_only=is_tgs_only, x_value=x_value, x_label=x_label, exp_name=exp_name)
 

  # Plot other measurements
  for p in plots:
    print p[1], p[2]
    plot_meas_fig(plots[p], stencil=p[1], plt_key=p[2], machine_name=machine_name, is_tgs_only=is_tgs_only, x_value=x_value, x_label=x_label, exp_name=exp_name)


def plot_perf_fig(p, stencil, machine_name, is_tgs_only, x_value, x_label, exp_name):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import matplotlib
  import pylab, itertools
  from pylab import arange,pi,sin,cos,sqrt
  from scripts.utils import get_stencil_num
  from scripts.plot_utils import method_style, hw_ctr_labels, marker_s, line_w, line_s

  f_name = stencil+'_'+exp_name
  if(is_tgs_only==1):
    f_name = f_name+'_tgs'

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
  plt.xlabel(x_label)
#  plt.legend(loc='best')
  plt.gca().set_ylim(bottom=0)
  pylab.savefig(machine_name + '_perf_' + f_name + '.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
  plt.clf()


def plot_meas_fig(p, stencil, plt_key, machine_name, is_tgs_only, x_value, x_label, exp_name):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import matplotlib
  import pylab
  from pylab import arange,pi,sin,cos,sqrt
  from scripts.utils import get_stencil_num
  from scripts.plot_utils import method_style, hw_ctr_labels, marker_s, line_w, line_s

  f_name = stencil+'_'+exp_name
  if(is_tgs_only==1):
    f_name = f_name+'_tgs'

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
      x = p[l][x_value]
      y = p[l][measure]
      plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=label)

    plt.ylabel(y_label)
    plt.grid()
    plt.xlabel(x_label)
    if(file_prefix=='mem_bw_' or file_prefix=='tlb_'):
      # Sort legends and use them
      ax = plt.gca()
      handles, labels = ax.get_legend_handles_labels()
      lab_int = []
      for s in labels:
        lab_key = 100
        if(s=='Spt.blk.'):
          lab_key = 101
        elif(s=='MWD'):
          lab_key = -1
        elif(s=='PLUTO'):
          lab_key = 102
        elif(s=='Pochoir'):
          lab_key = 103
        else:
          lab_key = int(s.split('W')[0])
        lab_int.append(lab_key)

      hl = sorted(zip(handles, labels, lab_int), key=itemgetter(2))
      handles2, labels2, dummy = zip(*hl)
      ax.legend(handles2, labels2, loc='best')

    plt.gca().set_ylim(bottom=0)
    if (plt_key in ['TLB', 'ENERGY']):
      plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    pylab.savefig(machine_name + '_' + file_prefix + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    plt.clf()


  # PLUTO tiling information
  if (any(method[1] == 'PLUTO' for method in p) and plt_key=='MEM'):
    plt.figure(plt_idx)
    plt_idx = plt_idx + 1
    for l in p:
      method=l[1]
      if(method == 'PLUTO'):
        measure, col, marker, label = ('pluto bs_x', 'm', '*', 'X')
        x = p[l][x_value]
        y = p[l][measure]
        plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=label)

    plt.ylabel('PLUTO tile size')
    plt.grid()
    plt.xlabel(x_label)
    plt.legend(loc='upper left')
    plt.gca().set_ylim(bottom=0)
    pylab.savefig(machine_name + '_pluto_tile_size_' + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    plt.clf()


  # Tiling information
  if (any(method[1] in ['MWD', '1WD'] for method in p) and plt_key=='MEM'):

    #Cache block size and diamond width
    for measure, y_label, f_prefix in [('blk size', 'Cache block size (MiB)', 'cache_block_size_'),
                                        ('diam width', 'Diamond width', 'diamond_width_'),
                                        ('wavefront width', 'Wavefront width', 'wavefront_width_')]:
      plt.figure(plt_idx)
      plt_idx = plt_idx + 1
      for l in p:
        method=l[1]
        tgsl = l[2]
        if((method in ['MWD', '1WD']) or (method =='PLUTO' and measure!='blk size')):
          if(is_tgs_only==1):
            method= str(tgsl) + 'WD'
          col, marker = method_style[method]
          x = p[l][x_value]
          y = p[l][measure]
          plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=method)

      plt.ylabel(y_label)
      plt.grid()
      plt.xlabel(x_label)
#      plt.legend(loc='best')
      plt.gca().set_ylim(bottom=0)
      pylab.savefig(machine_name + '_' + f_prefix + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
      plt.clf()

    if(is_tgs_only==1): return
    # MWD Thread group size information
    plt.figure(plt_idx)
    plt_idx = plt_idx + 1
    for l in p:
      method=l[1]
      if(method == 'MWD'):
        tgs_labels =[
               ('tgs', 'b', '^', 'MWD Group'),
               ('thx', 'r', '+', 'Along x'),
               ('thz', 'm', '*', 'Along z') ]
        if(stencil == 'solar'):
          tgs_labels = tgs_labels + [('thc', 'g', 'o', 'In comp.')]
        else:
          tgs_labels = tgs_labels + [('thy', 'g', 'o', 'Along y')]
        for measure, col, marker, label in tgs_labels:
          x = p[l][x_value]
          y = p[l][measure]
          plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=label)

    plt.ylabel('Intra-tile threads')
    plt.grid()
    plt.xlabel(x_label)
    plt.legend(loc='upper left')
    plt.gca().set_ylim(bottom=0)
    pylab.savefig(machine_name + '_thread_group_size_' + f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    plt.clf()


if __name__ == "__main__":
  main()
