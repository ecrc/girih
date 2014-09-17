#!/usr/bin/env python
def main():
  import matplotlib.pyplot as plt
  import pylab
  import matplotlib
  from pylab import arange,pi,sin,cos,sqrt

  fig_width = 4.3*0.393701 # inches
  fig_height = 0.6*fig_width #* 210.0/280.0#433.62/578.16

  fig_size =  [fig_width,fig_height]
  params = {
         'axes.labelsize': 7,
         'axes.linewidth': 0.5,
         'lines.linewidth': 0.75,
         'text.fontsize': 7,
         'legend.fontsize': 7,
         'xtick.labelsize': 7,
         'ytick.labelsize': 7,
         'lines.markersize': 3,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)

  x_7pt_const = [0,
  4, 
  8,
  12,
  16,
  20,
  24,
  32]

  meas_7pt_const = [24.5,
15.1,
8.1,
5.4,
4.1,
3.3,
2.7,
2.0]

  model_7pt_const = [ 24.0,
16.0,
8.0,
5.3,
4.0,
3.2,
2.7,
2.0]

  plt.plot(x_7pt_const, meas_7pt_const,  linestyle='-', color='k', marker='o', label='Measured')
  plt.plot(x_7pt_const, model_7pt_const, linestyle='--', color='b', marker='^', label='Model')
  
  plt.legend(loc='best')
  plt.xlabel('Diamond width')
  plt.ylabel('Bytes/LUP')

  plt.grid()

  pylab.savefig('7pt_const_byte_req.pdf', format='pdf', bbox_inches="tight", pad_inches=0)

  plt.close()
###############################################################
  x_7pt_var = [0,
4,
8,
12,
16,
20,
24]
 
  meas_7pt_var = [81.3,
45.5,
23.0,
15.2,
11.3,
9.1,
7.6]

  model_7pt_var = [80,
44.0,
22.0,
14.7,
11.0,
8.8,
7.3]

  plt.plot(x_7pt_var, meas_7pt_var,  linestyle='-', color='k', marker='o', label='Measured')
  plt.plot(x_7pt_var, model_7pt_var, linestyle='--', color='b', marker='^', label='Model')
  
  #plt.legend(loc='best')
  plt.xlabel('Diamond width')
  plt.ylabel('Bytes/LUP')

  plt.grid()

  pylab.savefig('7pt_var_byte_req.pdf', format='pdf', bbox_inches="tight", pad_inches=0)



if __name__ == "__main__":
  main()

