#!/usr/bin/env python
def main():
    import sys

    raw_data = load_csv(sys.argv[1])

    k_l = set()
    for k in raw_data:
        k_l.add(get_stencil_num(k))
    k_l = list(k_l)

    for k in k_l:
        for is_dp in [1]:
            plot_lines(raw_data, k, is_dp)


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

        
def plot_lines(raw_data, stencil_kernel, is_dp):
    from operator import itemgetter
    import matplotlib.pyplot as plt
    import matplotlib
    import pylab
    import numpy as np
    from pylab import arange,pi,sin,cos,sqrt

    fig_width = 5.2*0.393701 # inches
    fig_height = 0.95*fig_width #* 210.0/280.0#433.62/578.16

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

    req_fields = [('Thread group size', int), ('WD main-loop RANK0 MStencil/s  MAX', float), ('Time stepper orig name', str), ('MPI size', int), ('MStencil/s  MAX', float)]
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
        data.append(tup)

    data = sorted(data, key=itemgetter('Time stepper orig name', 'Thread group size', 'MPI size'))
    #for i in data: print i

    x = []
    y = []
    max_x = 0
    max_single = 0;
    for k in data:
        if ( ('Naive' in k['Time stepper orig name']) and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp)):
            if k['MPI size'] == 1 and max_single < k['MStencil/s  MAX']: max_single = k['MStencil/s  MAX']
            x.append(k['MPI size'])
            y.append(k['MStencil/s  MAX']*k['MPI size']/10**3)
    marker = 'o'
    col = 'g'
    ts2 = 'Spt.blk.'
    if(x):
        plt.plot(x, y, color=col, marker=marker, linestyle='-', label=ts2)
        max_y = max(y)
        max_x = max(x)
 

    perf_str = 'WD main-loop RANK0 MStencil/s  MAX'
    cols = {0:'y', 1:'k', 2:'b', 5:'r', 10:'m'}
    markers = {0:'.', 1:'^', 2:'v', 4:'.', 5:'x', 10:'*'}
    for tgs in [1,2,5, 10]:
        marker = markers[tgs]
        x = []
        y = []
        y_m = []
        for k in data:
            if ( ('Diamond' in k['Time stepper orig name']) and (k['Thread group size'] == tgs)  and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp) ):
                if k['MPI size'] == 1 and max_single < k[perf_str]: max_single = k[perf_str]
                x.append(k['MPI size'])
                y.append(k[perf_str]*k['MPI size']/10**3)
        col = cols[tgs]
        ts2 = str(tgs) + 'WD'
        if(x):
            plt.plot(x, y, color=col, marker=marker, linestyle='-', label=ts2)
            max_y = max(max(y), max_y)
            max_x = max(max(x), max_x)

    #title = 'Distributed memory strong scaling performance'
    plt.ylabel('Aggregate GLUP/s')
    #plt.title(title)
    title = ''

    if stencil_kernel == 0:
        title = title + ' 25_pt constant coeff. star stencil'
    elif stencil_kernel == 1:
        title = title + ' 7_pt constant coeff. star stencil'
    elif stencil_kernel == 2:
        title = title + ' 7_pt varialbe. coeff. star stencil'
    elif stencil_kernel == 3:
        title = title + ' 7_pt varialbe. coeff. axis symm star stencil'
    elif stencil_kernel == 4:
    #    title = title + ' 25_pt varialbe. coeff. axis symm star stencil'
        title = '25_pt_var_dist'
    elif stencil_kernel == 5:
    #    title = title + ' 7_pt varialbe. coeff. no symm star stencil'
        title = '7_pt_var_dist'

    np_l = [k['MPI size'] for k in data if k['stencil'] == stencil_kernel]
    np_set = sorted(list(set(np_l)))
 
    # add ideal scaling
    max_single = max_single/10**3
    np_l = [k['MPI size'] for k in data if k['stencil'] == stencil_kernel]
    plt.plot([1, max(np_l)], [max_single, max_single*max(np_l)], color='c', linestyle='--', label='Ideal scaling')

    f_name = title.replace(' ', '_')

#    plt.xlim([0.75, len(np_set)+0.25])
#    plt.ylim([0, max_single*max(np_l)])

    plt.xlabel('Processors \#')

    plt.xticks(np_set[1:], np_set[1:])

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
