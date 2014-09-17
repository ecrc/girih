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

    m = 3.0
  
    fig_width = 4.0*0.393701*m # inches
    fig_height = 1.0*fig_width #* 210.0/280.0#433.62/578.16

    fig_size =  [fig_width,fig_height]
    params = {
           'axes.labelsize': 6*m,
           'axes.linewidth': 0.25*m,
           'lines.linewidth': 0.75*m,
           'text.fontsize': 7*m,
           'legend.fontsize': 5*m,
           'xtick.labelsize': 6*m,
           'ytick.labelsize': 6*m,
           'lines.markersize': 1,
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

    req_fields = [('Thread group size', int), ('WD main-loop RANK0 MStencil/s  MAX', float), ('Time stepper orig name', str), ('OpenMP Threads', int), ('MStencil/s  MAX', float), ('Time unroll',int)]
    req_fields = req_fields + [ ('Sustained Memory BW', float)]
    req_fields = req_fields + [('Local NX', int), ('Local NY', int), ('Local NZ', int)]
    data = []
    for k in raw_data:
        tup = {}
        # add the general fileds
        for f in req_fields:
            try:
                tup[f[0]] = map(f[1], [k[f[0]]] )[0]
            except:
                print("ERROR: results entry missing essential data")
                print k
                return

        # add the stencil operator
        tup['stencil'] = get_stencil_num(k)

        # add the precision information
        if k['Precision'] in 'DP':
            p = 1
        else:
            p = 0
        tup['Precision'] = p
        data.append(tup)

    data = sorted(data, key=itemgetter('Time stepper orig name', 'Thread group size', 'Local NY'))
    #for i in data: print i

    max_single = 0
#    fig, ax1 = plt.subplots()
#    lns = []

    markers = 'o^v*x'
    marker_s = 7
    line_w = 1
    line_s = '-' 
 
    x = []
    y = []
    y_m = []
    max_x = 0
    max_y = 0
    for k in data:
        if ( ('Naive' in k['Time stepper orig name']) and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp)):
            y_m.append(k['Sustained Memory BW']/10**3)
            x.append(k['Local NY'])
            y.append(k['MStencil/s  MAX']/10**3)
    marker = markers[0]
    col = 'g'
    ts2 = 'Spt.blk.'
    if(x) and t==0:
        plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=ts2)
        max_y = max(y)
        max_x = max(x)
    if(y_m) and t==1:
        plt.plot(x, y_m, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=ts2)
 

    perf_str = 'WD main-loop RANK0 MStencil/s  MAX'
    cols = {0:'y', 1:'k', 2:'b', 5:'c', 10:'m'}
    for idx, tgs in enumerate([1,2,5, 10]):
        marker = markers[1+idx]
        x = []
        y = []
        y_m = []
        for k in data:
            if ( ('Diamond' in k['Time stepper orig name']) and (k['Thread group size'] == tgs)  and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp) ):
                y_m.append(k['Sustained Memory BW']/10**3)
                x.append(k['Local NY'])
                y.append(k[perf_str]/10**3)
        col = cols[tgs]
        ts2 = str(tgs) + 'WD'
        if(x) and t==0:
            plt.plot(x, y, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=ts2)
            max_y = max(max(y), max_y)
            max_x = max(max(x), max_x)
        if(y_m) and t==1:
            plt.plot(x, y_m, color=col, marker=marker, markersize=marker_s, linestyle=line_s, linewidth=line_w, label=ts2)

    # add limits
    sus_mem_bw = 40.0 #IB

    # divide by streams number
    if stencil_kernel == 0:
        mem_limit = sus_mem_bw/4.0
    elif stencil_kernel == 1:
        mem_limit = sus_mem_bw/3.0
    elif stencil_kernel == 2:
        mem_limit = sus_mem_bw/5.0
    elif stencil_kernel == 3:
        mem_limit = sus_mem_bw/6.0
    elif stencil_kernel == 4:
        mem_limit = sus_mem_bw/16.0
    elif stencil_kernel == 5:
        mem_limit = sus_mem_bw/10.0

    # divide by word size
    if is_dp == 1: 
        mem_limit = mem_limit/8.0
    else:
        mem_limit = mem_limit/4.0

    if t == 0:
        plt.plot([1, max_x], [mem_limit, mem_limit], color='.2', linestyle='--', label='Spt.lim.')


    if t == 0:
        title = '_strongscaling_perf' #'Strong scaling performance'
        if stencil_kernel == 1: plt.ylabel('GLUP/s')
    else:
        title = '_strongscaling_bw' #'Strong scaling main memory BW usage'
        if stencil_kernel == 1: plt.ylabel('GBytes/s')

    if stencil_kernel == 0:
        title = '25_pt_const' + title
    elif stencil_kernel == 1:
        title = '7_pt_const' + title
    elif stencil_kernel == 4:
        title = '25_pt_var' + title
    elif stencil_kernel == 5:
        title = '7_pt_var' + title

    f_name = title.replace(' ', '_')

    if t==0: plt.ylim([0,max_y*1.1])

    plt.xlabel('Size in each dimension')
 
    if t == 0 and stencil_kernel == 1: plt.legend(loc='best')
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
