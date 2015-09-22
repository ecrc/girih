#!/usr/bin/env python

def main():
  import sys
  from utils import load_csv, get_stencil_num

  raw_data = load_csv(sys.argv[1])

  k_l = set()
  for k in raw_data:
    k_l.add((get_stencil_num(k), k['Global NX']))
  k_l = list(k_l)

  bsz_l = set()
  for k in raw_data:
    if k['Multi-wavefront updates']=='0': continue
    bsz_l.add(k['Multi-wavefront updates'])
  bsz_l = sorted(list(bsz_l))

  for k, N in k_l:
    for bsz in bsz_l:
      gen_res(raw_data, int(k), int(bsz), int(N))

#    gen_res(raw_data)

def gen_res(raw_data, stencil_kernel, bsz,  N):
  from operator import itemgetter
  import matplotlib.pyplot as plt
  import pylab
  from csv import DictWriter
  from operator import itemgetter
  from utils import get_stencil_num

  #fig_width = 8.588*0.393701 # inches
  fig_width = 4.0*0.393701 # inches
  fig_height = 1.0*fig_width #* 210.0/280.0#433.62/578.16

  fig_size =  [fig_width,fig_height]
  params = {
         'axes.labelsize': 5,
         'axes.linewidth': 0.5,
         'lines.linewidth': 1,
         'font.size': 5,
         'legend.fontsize': 5,
         'xtick.labelsize': 5,
         'ytick.labelsize': 5,
         'lines.markersize': 5,
         'text.usetex': True,
         'figure.figsize': fig_size}
  pylab.rcParams.update(params)


  req_fields = [('Total cache block size (kiB)', int), ('MStencil/s  MAX', float), ('Time stepper orig name', str), ('Stencil Kernel semi-bandwidth', int), ('Stencil Kernel coefficients', str), ('Precision', str), ('Time unroll',int), ('Number of time steps',int), ('Number of tests',int), ('Local NX',int), ('Local NY',int), ('Local NZ',int), ('Total Memory Transfer', float), ('Thread group size' ,int), ('Intra-diamond prologue/epilogue MStencils',int), ('Multi-wavefront updates', int), ('Intra-diamond width', int)]
  data = []
  for k in raw_data:

    # Use single field to represent the performance
    if 'Total RANK0 MStencil/s MAX' in k.keys():
      if(k['Total RANK0 MStencil/s MAX']!=''):
        k['MStencil/s  MAX'] = k['MWD main-loop RANK0 MStencil/s MAX']
    # temporary for deprecated format
    if 'RANK0 MStencil/s  MAX' in k.keys():
      if k['RANK0 MStencil/s  MAX']!='':
        k['MStencil/s  MAX'] = k['RANK0 MStencil/s  MAX']

    tup = dict()
    # add the general fileds
    for f in req_fields:
      try:
        tup[f[0]] = map(f[1], [k[f[0]]] )[0]
      except:
        print f[0]
    # add the stencil operator
    tup['Kernel'] = get_stencil_num(k)
    data.append(tup)

    # add precision information
    tup['word size'] = 8 if k['Precision'] in 'DP' else 4

  #for i in data: print i

  for tup in data:
    tup['Actual Bytes/LUP'] = actual_BpU(tup)
    tup['Model'] = models(tup)
    # model error
    tup['Err %'] = 100 * (tup['Model'] - tup['Actual Bytes/LUP'])/tup['Actual Bytes/LUP']
    tup['D_width'] = tup['Intra-diamond width']
    tup['bsz'] = tup['Multi-wavefront updates']
    tup['Performance'] = tup['MStencil/s  MAX']
#    tup['Cache block'] = get_bs(Dw=tup['D_width'], Nd=get_nd(tup['Kernel']), Nf=tup['Multi-wavefront updates'], Nx=tup['Local NX'], WS=tup['word size'], R=tup['Stencil Kernel semi-bandwidth'])
    tup['Cache block'] =  tup['Total cache block size (kiB)']/1024.0
#    print "cache block C:", tup['Total cache block size (kiB)']/1024.0, " Python:", tup['Cache block']
#    try: print "%6.3f  %6.3f  %6.3f" % (tup['Cache block'], tup['Total cache block size (kiB)']/1024.0,tup['Cache block']- tup['Total cache block size (kiB)']/1024.0)
#    except: pass

  data = sorted(data, key=itemgetter('Kernel', 'Local NX', 'D_width'))
  #for i in data: print i['Kernel'], i['Local NX'], i['D_width'], i['bsz']

  fig, ax = plt.subplots()
  cs=[]
  cb=[]
  cb_meas=[]
  Dw=[]
  for k in data:
    if k['Kernel']==stencil_kernel and k['Local NX']==N and (k['bsz']==bsz or k['Time stepper orig name'] == 'Spatial Blocking'):
      cs.append(k['Cache block'])
      cb.append(k['Model'])
      cb_meas.append(k['Actual Bytes/LUP'])
      Dw.append(k['D_width'])

  #for i in range(len(cs)):
  #  print Dw[i], cs[i], cb_meas[i], cb[i]
  if Dw==[]: return
  ax.plot(cs, cb     , marker='+', linestyle='-', color='k', label="Model")
  ax.plot(cs, cb_meas, marker='x', linestyle='--', color='b', label="Measured")


  # show the usable cache size limits
#  ax.plot([22.5, 22.5], [0, 0.7*cb[0]], linestyle='-', color='r', label="Usable cache size")

  ax.set_ylabel('Code balance (Bytes/LUP)')
  ax.set_xlabel('Cache block size (MiB) PER THREAD')
  ax.set_ylim([0, max(cb_meas+cb)+1])
  ax.set_xlim([0, max(cs)+0.5])
  ax2 = ax.twiny()
  ax2.set_xticks(cs)
  ax2.set_xlabel('Diamond width')
  ax2.set_xlim(ax.get_xlim())

  if stencil_kernel==1:
    Dw = map(str,Dw)
    Dw[1]=''
    Dw[3]=''
    Dw[5]=''
  ax2.set_xticklabels(Dw)

#  for i, d in enumerate(Dw):
    #if ((d+4)%8 == 0):
#    ax.annotate(d, (cs[i], cb[i]))  

  title = '_code_balance_vs_cache_size_N%d_bsz%d'%(N, bsz)
  if stencil_kernel == 0:
      title = '25_pt_const' + title
  elif stencil_kernel == 1:
      title = '7_pt_const' + title
  elif stencil_kernel == 4:
      title = '25_pt_var' + title
  elif stencil_kernel == 5:
      title = '7_pt_var' + title
  elif stencil_kernel == 6:
      title = 'solar' + title


  ax.legend(loc='best')
  ax.grid()
  pylab.savefig(title+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
  plt.clf()



def actual_BpU(tup):
  total_mem = tup['Total Memory Transfer']
  R = tup['Stencil Kernel semi-bandwidth']
  nt = tup['Number of time steps']

  nx = tup['Local NX']
  ny = tup['Local NY']
  nz = tup['Local NZ']

  oh = 0#tup['Intra-diamond prologue/epilogue MStencils']

  stencil_size = nx*ny*nz#2*ny*nz + ny*nz*(nx+2*R) 
  BpU = (total_mem * 10**9) / ( (stencil_size * nt - oh*10**6)*tup['Number of tests'])

  #print BpU, total_mem, stencil_size, nt, oh, tup['Number of tests']
  return BpU

#def get_bs(Dw, Nd, Nf, Nx, WS):
#  Ww = Dw + Nf - 1
#  Bs = WS*Nx*( Nd*(Dw**2/2.0 + Dw*Nf) + 2.0*(Dw+Ww) )
#  if Dw==0: Bs=0
#  return Bs/(1024.0*1024.0)


def get_nd(k):
  if k== 0:
    nd = 2 + 1
  if k== 1:
    nd = 2
  if k== 4:
    nd = 2 + 13
  if k== 5:
    nd = 2 + 7
  if k== 6:
    nd = 40
  return nd


def models(tup):
  if   tup['Precision'] == 'DP': WS = 8
  elif tup['Precision'] == 'SP': WS = 4

  # number of streamed copies of the domain (buffers)
  nb = get_nd(tup['Kernel'])

  # Spatial blocking
  if tup['Time stepper orig name'] == 'Spatial Blocking':
    if 'Solar' in tup['Stencil Kernel coefficients']:
      return (12 + nb) * WS
    else:
      return (1 + nb) * WS
 

  R = tup['Stencil Kernel semi-bandwidth']
  Dw = tup['Intra-diamond width']

  if 'Solar' in tup['Stencil Kernel coefficients']:
    bpu = 32. * (6*(2*Dw-1) + nb*Dw+2*6) / (float(Dw)**2)
  else:
    bpu = 2.*WS * R * (2*Dw-2*R + nb*Dw+2*R) / (float(Dw)**2)

  return bpu



if __name__ == "__main__":
  main()

