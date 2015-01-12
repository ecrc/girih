
def get_bs(Dw, Nd, Nf, Nx, WS):
  Ww = Dw + Nf - 1
  Bs = WS*Nx*( Nd*(Dw**2/2.0 + Dw*Nf) + 2.0*(Dw+Ww) )
  if Dw==0: Bs=0
  return Bs/(1024.0*1024.0)


def get_nd(k):
  if k== 0:
    nd = 2 + 1
  if k== 1:
    nd = 2
  if k== 4:
    nd = 2 + 13
  if k== 5:
    nd = 2 + 7

  return nd


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


def models(tup):
  if   tup['Precision'] == 'DP': word_size = 8
  elif tup['Precision'] == 'SP': word_size = 4

  R = tup['Stencil Kernel semi-bandwidth']
  TB = tup['Time unroll']
  ny = tup['Local NY']

  # number of streamed copies of the domain (buffers)
  nb = get_nd(tup['Kernel'])

  width = tup['Intra-diamond width']
  YT_section = float((TB+1)**2 * 2 * R)

  # temporal blocking model
  if tup['Time stepper orig name'] == 'Spatial Blocking':
    bpu = (1 + nb) * word_size
  else: # no temporal blocking model
    bpu = ( ((width - 2*R) + width) + (nb*width + 2*R) ) * word_size / YT_section

  return bpu


