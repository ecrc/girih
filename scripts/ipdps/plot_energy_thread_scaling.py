#!/usr/bin/env python
def main():
    import sys

    raw_data = load_csv(sys.argv[1])

    k_l = set()
    for k in raw_data:
        k_l.add(get_stencil_num(k))
    k_l = list(k_l)

 #   for ts in ['Naive', 'Dynamic-Intra-Diamond']  
    is_dp = 1
    for k in k_l:
        for val in ['Energy', 'Power']:
            plot_lines(raw_data, k, is_dp, val)

 
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
    
def plot_lines(raw_data, stencil_kernel, is_dp, val):
    from operator import itemgetter
    import matplotlib.pyplot as plt
    import matplotlib
    import pylab
    from pylab import arange,pi,sin,cos,sqrt

    fig_width = 4.0*0.393701 # inches
    fig_height = 1.0*fig_width #* 210.0/280.0#433.62/578.16

    fig_size =  [fig_width,fig_height]
    params = {
           'axes.labelsize': 7,
           'axes.linewidth': 0.5,
           'lines.linewidth': 0.75,
           'text.fontsize': 7,
           'legend.fontsize': 5,
           'xtick.labelsize': 7,
           'ytick.labelsize': 7,
           'lines.markersize': 3,
           'text.usetex': True,
           'figure.figsize': fig_size}
    pylab.rcParams.update(params)

    ts_l = set()
    for k in raw_data:
        ts_l.add(k['Time stepper orig name'])
    ts_l = list(ts_l)

    th = set()
    for k in raw_data:
        th.add(int(k['OpenMP Threads']))
    th = list(th)

    tb_l = set()
    for k in raw_data:
        tb_l.add(k['Time unroll'])
    tb_l = list(tb_l)
    tb_l = map(int,tb_l)
    tb_l.sort()
  
    tgs_l = set()
    for k in raw_data:
        tgs_l.add(k['Thread group size'])
    tgs_l = list(tgs_l)
    tgs_l = map(int,tgs_l)
    tgs_l.sort()

    req_fields = [('Thread group size', int), ('WD main-loop RANK0 MStencil/s  MAX', float), ('Time stepper orig name', str), ('OpenMP Threads', int), ('MStencil/s  MAX', float), ('Time unroll',int), ('Energy', float), ('Energy DRAM', float), ('Power',float), ('Power DRAM', float), ('Intra-diamond prologue/epilogue MStencils', int), ('Global NX', int), ('Number of time steps', int)]
    data = []
    for k in raw_data:
        if k['Intra-diamond prologue/epilogue MStencils'] == '':
            k['Intra-diamond prologue/epilogue MStencils'] = 0
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
        data.append(tup)

    data = sorted(data, key=itemgetter('Time stepper orig name', 'stencil', 'Thread group size', 'OpenMP Threads'))
#    for i in data: print i

    max_single = 0
    cols = {0:'g', 1:'k', 2:'b', 4:'c', 5:'r', 8:'m', 10:'m'}
    markers = {0:'o', 1:'^', 2:'v', 4:'.', 5:'x', 8:'.', 10:'*'}
    for tgs in [ 10, 1, 2, 4, 5, 8]:
        marker = markers[tgs]
        x = []
        y = []
        y1 = []
        y2 = []
        cpuv = val
        memv = val + ' DRAM'
        scale = 1
        for k in data:
            if ( (k['Thread group size'] == tgs)  and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp) ):
                ext_lups = k['Intra-diamond prologue/epilogue MStencils']*1e6
                NT = k['Global NX']**3 * k['Number of time steps']
                Neff = NT - ext_lups
                if 'Energy' in val:
                    scale = Neff / 1e9
                y.append((k[cpuv]+k[memv])/scale)
                y1.append(k[cpuv]/scale)
                y2.append(k[memv]/scale)
                x.append(k['OpenMP Threads'])
        col = cols[tgs]
        ts2 = str(tgs) + 'WD' if tgs !=0 else 'Spt.blk.'
        markersize = 12 if tgs == 10 else 4
        if(x):
            plt.plot(x, y, color=col, marker=marker, markersize=markersize, linestyle='-', label=ts2)

 
    title = '_threadscaling_'+val.lower()
    if stencil_kernel == 0:
        title = '25_pt_const' + title
    elif stencil_kernel == 1:
        title = '7_pt_const' + title
    elif stencil_kernel == 4:
        title = '25_pt_var' + title
    elif stencil_kernel == 5:
        title = '7_pt_var' + title

    if 'Energy' in val:
        plt.ylabel('pJ/LUP')
    elif 'Power' in val:
        plt.ylabel('W')

    f_name = title

    plt.xlabel('Threads')
 
    plt.legend(loc='best')
    plt.grid()
    pylab.savefig(f_name+'.png', bbox_inches="tight", pad_inches=0.04)
    pylab.savefig(f_name+'.pdf', format='pdf', bbox_inches="tight", pad_inches=0)
    #plt.show()     
    plt.clf()
        
def load_csv(data_file):
    from csv import DictReader
    with open(data_file, 'rb') as output_file:
        data = DictReader(output_file)
        data = [k for k in data]
    return data
    
    
if __name__ == "__main__":
    main()
