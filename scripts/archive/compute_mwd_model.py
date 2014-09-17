#!/usr/bin/env python


def plot_series(exp, f_name):
  import matplotlib.pyplot as plt
  import matplotlib
  import pylab

  matplotlib.rcParams.update({'font.size': 16})

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
  plt.ylabel("Cache block size (MiB)") 
  plt.legend(loc='best')
  pylab.savefig(f_name+'.png')
  pylab.savefig(f_name+'.pdf', format='pdf')

  plt.clf()



def get_blk_size(N, r, diam_width, nwf, nbuf, word_size):
  t_dim = diam_width/(2*r) - 1
  diam_h = t_dim*2*r + nwf;

  wf_updates = (t_dim+1)*(t_dim+1)*2 * r; # Y-T projection
  wf_elements = (wf_updates - diam_width) * r + diam_width + diam_width*(nwf-1);

  return (nbuf*wf_elements + (diam_width + diam_h )*2*r) * N * word_size


def get_max_diam(N, CS, r, nbuf, word_size):
  nwf = 1
  ntg = 1

  max_diam_width = 0;
  for i in range(N/r, 0, -4):
    lt_dim = i/2 - 1;
    diam_width = i*r;
    diam_height = lt_dim*2*r+1 + nwf-1;

    wf_size = get_blk_size(N, r, diam_width, nwf, nbuf, word_size)
    cache_size_cond = wf_size*ntg < CS*1024;

    cuncurrency_cond = (N/diam_width) >= ntg;
    int_diam_cond = N%diam_width == 0;
    wf_len_cond = diam_height <= N;

#    print("i:%d, diam_width %d,  cuncurrency_cond %d, cache_size_cond %d, int_diam_cond %d, wf_len_cond %d, cache_blk_size: %lu kB" %
#        (i, diam_width, cuncurrency_cond, cache_size_cond, int_diam_cond, wf_len_cond, wf_size*ntg/1024))
    if( (int_diam_cond == 1) and (wf_len_cond == 1) and (cuncurrency_cond == 1)  and (cache_size_cond == 1) ):
      max_diam_width = diam_width
      break

  return max_diam_width


def models(r, diam_width, nbuf, word_size):
  YT_section = float((diam_width)**2)/( 2 * r)

  # no temporal blocking model
  if diam_width == 0:
    BpLUP = (1 + nbuf) * word_size
  else: # temporal blocking model
    BpLUP =( diam_width-2*r + diam_width + nbuf*diam_width + 2*r ) * word_size / YT_section

  return BpLUP


def main():    

  word_size = 8  # DP

  # 7pt const. coeff.
  r = 1
  coeff_type = 'const' 

  if coeff_type == 'const':         nbuf = 2
  elif coeff_type == 'var_all_sym': nbuf = 2 + 1 + r
  elif coeff_type == 'var_ax_sym':  nbuf = 2 + 1 + 3*r
  elif coeff_type == 'var_no_sym':  nbuf = 2 + 1 + 6*r


  print("N     Dw  Bytes/LUP ")

  # spatial blocking case
  BpLUP = models(r, 0, nbuf, word_size)
  print("%04d  %02d  %02d "%(0, 0, BpLUP))

  # loop over stencil type
  for CS in [512]: # loop over desired cache block sizes
    for N in range(32, 1025, 32): # loop over increasing grid size
      
      # find the largets diamond tile fitting in cache
      diam_width = get_max_diam(N, CS, r, nbuf, word_size)
      # compute the corresponding bytes requirement
      BpLUP = models(r, diam_width, nbuf, word_size)
      
      print("%04d  %02d  %02d "%(N, int(diam_width), BpLUP))
  
  return


if __name__ == "__main__":
  main()
