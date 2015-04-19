#!/usr/bin/env python

def plot_stacked_clustered_bars():
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import pylab
  import numpy as np
  from pylab import arange,pi,sin,cos,sqrt

  sec_fontsize = 13

  fig_width = 2*7.4*0.393701 # inches
  fig_height = 0.75*fig_width #* 210.0/280.0#433.62/578.16
  fig_size =  [fig_width,fig_height]
  params = {
         'axes.labelsize': 12,
#         'axes.linewidth': 0.5,
#         'lines.linewidth': 0.75,
         'text.fontsize': 12,
         'legend.fontsize': 12,
         'xtick.labelsize': 12,
         'ytick.labelsize': 12,
#         'lines.markersize': 3,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)


  stencils = ['7pt const.', '7pt var.', '25pt var.','','','','','','','','']
  processors = ['','Westmere','','', '','Ivy Bridge','', '','','Haswell','','']
  pochoir_perf = [806.0,  260.0,  93.0, 0,   2814.0, 753.0,  267.0, 0,   3716.0, 0966.0, 378.0]
  girih_perf   = [1475.0, 445.0, 148.0, 0,   3975.0, 1266.0, 363.0, 0,   6125.0, 1940.0, 559.0]
  pochoir_perf = [i/1e3 for i in pochoir_perf]
  girih_perf = [i/1e3 for i in girih_perf]

  speedup = [0]*len(stencils)
  for i in range(len(stencils)):
    speedup[i] = girih_perf[i]/pochoir_perf[i] if pochoir_perf[i] else 0.0
  nx = len(pochoir_perf)
  x = range(nx)
  cluster_size = 4

  fig = plt.figure()
  fig.subplots_adjust(bottom=0.25)  # space on the bottom for second time axis
  host = fig.add_subplot(111)  # setup plot  
  width =0.8

  p1 = host.bar(x, girih_perf,   width, color='0.65', align='center', hatch="")
  p2 = host.bar(x, pochoir_perf, width, color='0.95', align='center', hatch="/")

  for i, r in enumerate(p1):
    if speedup[i] > 0:
      height = r.get_height()
      host.text(r.get_x()+width/2., height+0.01, '%2.1fx'%speedup[i], ha='center', va='bottom')

  host.set_ylabel('GLUP/s', fontsize=sec_fontsize)
#  host.set_xlabel('Processor', fontsize=sec_fontsize)
  host.tick_params(axis='both', which='major', labelsize=sec_fontsize)
  host.tick_params(axis='both', which='minor', labelsize=sec_fontsize)

  xtk = np.arange(nx)
  for i in range(nx):
    if i%cluster_size >= 3:
      xtk[i]=0
  host.set_xticks(xtk)
  host.set_xticklabels(processors)

  # Insert the time steppers names at the X-axis 
  newax = host.twiny()  # create new axis
  newax.xaxis.set_ticks_position('bottom')
  newax.spines['bottom'].set_position(('outward', 20))

  #  newax.patch.set_visible(False)
  #  newax.xaxis.set_label_position('bottom')
  #  newax.set_frame_on(False)
  #  newax.tick_params('x', width=0)

  ticks = np.arange(nx)
  for i in range(nx):
    if (i >= 3):
      ticks[i] = 0
  newax.set_xticks(ticks)
  newax.set_xticklabels(stencils, rotation=90, size='medium')

  newax.axis((0.0, float(nx), 0.0, max(girih_perf)*1.1))

  newax.set_xlim((-width, nx))
  host.set_xlim((-width, nx))
 
  host.yaxis.grid() 

  #pylab.legend(p1[:nts], ts_set)
  host.legend( (p1, p2), ('MWD', 'Pochoir'), loc='center left', fontsize=sec_fontsize)

  f_name = "pochoir_comparison" 
  #pylab.savefig(f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
  pylab.savefig(f_name+'.eps', format='eps', bbox_inches="tight", pad_inches=0.02)

  return
      
if __name__ == "__main__":
  plot_stacked_clustered_bars()
