#!/usr/bin/env python
def main():
    import sys

    raw_data = load_csv(sys.argv[1])

    k_l = set()
    for k in raw_data:
        k_l.add(get_stencil_num(k))
    k_l = list(k_l)

 #   for ts in ['Naive', 'Dynamic-Intra-Diamond']  
    for k in k_l:
        for is_dp in [1]:
            for t in [0, 1]:
                plot_lines(raw_data, k, is_dp, t)

 
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
    
def plot_lines(raw_data, stencil_kernel, is_dp, t):
    from operator import itemgetter
    import matplotlib.pyplot as plt
    import matplotlib
    import pylab
    from pylab import arange,pi,sin,cos,sqrt

    fig_width = 3.8*0.393701 # inches
    fig_height = 1.0*fig_width #* 210.0/280.0#433.62/578.16

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

    req_fields = [('Thread group size', int), ('WD main-loop RANK0 MStencil/s  MAX', float), ('Time stepper orig name', str), ('OpenMP Threads', int), ('MStencil/s  MAX', float), ('Time unroll',int), ('Sustained Memory BW', float)]
    data = []
    for k in raw_data:
        tup = {}
        # add the general fileds
        for f in req_fields:
            tup[f[0]] = map(f[1], [k[f[0]]] )[0]

        # add the stencil operator
#        if  k['Stencil Kernel coefficients'] in 'constant':
#            if  int(k['Stencil Kernel semi-bandwidth'])==4:
#                stencil = 0
#            else:
#                stencil = 1
#        elif  'no-symmetry' in k['Stencil Kernel coefficients']:
#            stencil = 5
#        elif  'sym' in k['Stencil Kernel coefficients']:
#            if int(k['Stencil Kernel semi-bandwidth'])==1:
#                stencil = 3
#            else:
#                stencil = 4
#        else:
#            stencil = 2
#        tup['stencil'] = stencil
        tup['stencil'] = get_stencil_num(k)

        # add the precision information
        if k['Precision'] in 'DP':
            p = 1
        else:
            p = 0
        tup['Precision'] = p
        data.append(tup)

    data = sorted(data, key=itemgetter('Time stepper orig name', 'Time unroll', 'Thread group size', 'OpenMP Threads'))
#    for i in data: print i

    max_single = 0
#    fig, ax1 = plt.subplots()
#    lns = []

    marker = 'o'
    x = []
    y = []
    y_m = []
    for k in data:
        if ( ('Naive' in k['Time stepper orig name']) and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp)):
            if k['OpenMP Threads'] == 1 and max_single < k['MStencil/s  MAX']/10**3: max_single = k['MStencil/s  MAX']/10**3
            y_m.append(k['Sustained Memory BW']/10**3)
            x.append(k['OpenMP Threads'])
            y.append(k['MStencil/s  MAX']/10**3)
    marker = 'o'
    col = 'g'
    ts2 = 'Spt.blk.'
    if(x) and t==0:
        plt.plot(x, y, color=col, marker=marker, linestyle='-', label=ts2)
    if(y_m) and t==1:
        plt.plot(x, y_m, color=col, marker=marker, linestyle='-', label=ts2)

    x = []
    y = []
    y_m = []
    perf_str = 'WD main-loop RANK0 MStencil/s  MAX'
    for k in data:
        if ( ('Diamond' in k['Time stepper orig name']) and (k['Thread group size'] == 10) and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp)):
            y_m.append(k['Sustained Memory BW']/10**3)
            x.append(k['OpenMP Threads'])
            y.append(k[perf_str]/10**3)
    marker = '*'
    markersize = 12
    col = 'm'
    ts2 = str(10) + 'WD'
    if(x) and t==0:
        plt.plot(x, y, color=col, marker=marker, markersize=markersize,linestyle='', label=ts2)
    if(y_m) and t==1:
        plt.plot(x, y_m, color=col, marker=marker, markersize=markersize,linestyle='', label=ts2)




    cols = {0:'y', 1:'k', 2:'b', 4:'c', 5:'r', 8:'m'}
    markers = {0:'.', 1:'^', 2:'v', 4:'.', 5:'x', 8:'.'}
    for tgs in [1,2, 4, 8, 5]:
        marker = markers[tgs]
        x = []
        y = []
        y_m = []
        for k in data:
            if ( ('Diamond' in k['Time stepper orig name']) and (k['Thread group size'] == tgs)  and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp) ):
                if k['OpenMP Threads'] == 1 and max_single < k[perf_str]/10**3: max_single = k[perf_str]/10**3
                y_m.append(k['Sustained Memory BW']/10**3)
                x.append(k['OpenMP Threads'])
                y.append(k[perf_str]/10**3)
        col = cols[tgs]
        ts2 = str(tgs) + 'WD'
        if(x) and t==0:
            plt.plot(x, y, color=col, marker=marker, linestyle='-', label=ts2)
        if(y_m) and t==1:
            plt.plot(x, y_m, color=col, marker=marker, linestyle='-', label=ts2)



    # add limits
    mem_limit=0
#    sus_mem_bw = 36500 #SB
    sus_mem_bw = 40 #IB

    if stencil_kernel == 0:
        mem_limit = sus_mem_bw/16
    elif stencil_kernel == 1:
        mem_limit = sus_mem_bw/12
    elif stencil_kernel == 2:
        mem_limit = sus_mem_bw/20
    if is_dp == 1: mem_limit = mem_limit / 2
    if t == 0:
        #plt.plot([1, len(th)], [mem_limit, mem_limit], color='g', linestyle='--', label='Spatial blk. limit')
        pass

    # add ideal scaling
    ideal = [i*max_single for i in th]
    if t == 0:
        plt.plot(th, ideal, color='k', linestyle='--', label='Ideal scaling')

    if t == 0:
        title = '25_pt_var_all_methods_perf'
#        plt.ylabel('GLUP/s')
    else:
        title = '25_pt_var_all_methods_bw'
#        plt.ylabel('GBytes/s')

    #plt.title(title)

    f_name = title.replace(' ', '_')

    plt.xlabel('Threads')
 
#    if t == 0: plt.legend(loc='best')
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
