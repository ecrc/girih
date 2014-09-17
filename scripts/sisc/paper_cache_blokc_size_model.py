#!/usr/bin/env python

def get_blk_size(nx, to, r, diam_width, nwf, coeff_type, word_size):

  t_dim = diam_width/(2*r) - 1
  diam_h = t_dim*2*r + nwf;

  wf_updates = (t_dim+1)*(t_dim+1)*2 * r; # Y-T projection
  wf_elements = (wf_updates - diam_width) * r + diam_width + diam_width*(nwf-1);

  if coeff_type == 'const':
    total_points = ( (to+1)           *wf_elements + (diam_width + diam_h )*2*r) * nx * word_size

  elif coeff_type == 'var_all_sym':
    total_points = ( (to+1 + (1+  r) )*wf_elements + (diam_width + diam_h )*2*r) * nx * word_size

  elif coeff_type == 'var_ax_sym':
    total_points = ( (to+1 + (1+3*r) )*wf_elements + (diam_width + diam_h )*2*r) * nx * word_size

  elif coeff_type == 'var_no_sym':
    total_points = ( (to+1 + (1+6*r) )*wf_elements + (diam_width + diam_h )*2*r) * nx * word_size

  return total_points


def plot_series(exp, f_name, show_y_label):
  import matplotlib.pyplot as plt
  import matplotlib
  import pylab
  from pylab import arange,pi,sin,cos,sqrt

  fig_width = 4.3*0.393701 # inches
  fig_height = 0.6*fig_width #* 210.0/280.0#433.62/578.16

  fig_size =  [fig_width,fig_height]
  params = {
         'axes.labelsize': 6,
         'axes.linewidth': 0.5,
         'lines.linewidth': 0.75,
         'text.fontsize': 7,
         'legend.fontsize': 6,
         'xtick.labelsize': 7,
         'ytick.labelsize': 7,
         'lines.markersize': 2,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)

  max_width = 0
  for nx, r, coeff_type, to, name, m in exp:
    width_l = [i*4*r for i in range(1,21)]
    blk_l = [get_blk_size(nx=nx, to=to, r=r, diam_width=i, nwf=1, coeff_type=coeff_type, word_size=8) for i in width_l]
    blk_l = map(lambda x: x/(1024.0*1024.0), blk_l)
    try:
      max_idx = next(x[0] for x in enumerate(blk_l) if x[1] > 35)
    except:
      max_idx = len(blk_l)-1
    plt.plot(width_l[:max_idx], blk_l[:max_idx], marker=m, label=name)#+'@Nx:'+str(nx))
    if max_width < width_l[max_idx]:  max_width = width_l[max_idx] 

  plt.grid()
  plt.xticks(range(0,max_width, 8))
  plt.xlabel("Diamond tile width")
  if show_y_label ==1: plt.ylabel("Cache block size (MiB)") 
  plt.legend(loc='best')
  pylab.savefig(f_name+'.png', bbox_inches="tight", pad_inches=0.04)
  pylab.savefig(f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)

  plt.clf()


def main():
  import matplotlib.pyplot as plt

  exp = [(960, 1, 'const', 1, '$960$', '^'),
         (480, 1, 'const', 1, '$480$', 'o'),
         (240, 1, 'const', 1, '$240$', 'v')]
  plot_series(exp, '7_pt_const_cache_req', 1)

  exp = [(480, 4, 'var_ax_sym', 1, '$480$', '^'),
         (240, 4, 'var_ax_sym', 1, '$240$', 'o'),
         (120, 4, 'var_ax_sym', 1, '$120$', 'v')]
  plot_series(exp, '25_pt_var_cache_req', 0)

#  exp = [(960, 1, 'const', 1, '7pt_const', '^'),
##         (960, 4, 'const', 2, '25pt_const', 'v'),
#         (680, 1, 'var_no_sym', 1, '7pt_var', 'o'),
#         (480, 4, 'var_ax_sym', 1, '25pt_var', '*')]
#  plot_series(exp, 'all_cache_req_large')
#
#  exp = [(480, 1, 'const', 1, '7pt_const', '^'),
##         (480, 4, 'const', 2, '25pt_const', 'v'),
#         (340, 1, 'var_no_sym', 1, '7pt_var', 'o'),
#         (240, 4, 'var_ax_sym', 1, '25pt_var', '*')]
#  plot_series(exp, 'all_cache_req_small')



  return


if __name__ == "__main__":
  main()
