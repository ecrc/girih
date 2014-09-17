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
    from pylab import arange,pi,sin,cos,sqrt

    if stencil_kernel == 1: 
        fig_width = 3.4*0.393701 # inches
    else:
        fig_width = 3.0*0.393701 # inches

    fig_height = 1.8*fig_width #* 210.0/280.0#433.62/578.16

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

    #tb_l = [3, 7]
    tb_l = set()
    for k in raw_data:
        tb_l.add(k['Time unroll'])
    tb_l = list(tb_l)
    tb_l = map(int,tb_l)
    tb_l.sort()
    #print tb_l

    req_fields = [('Thread group size', int), ('WD main-loop RANK0 MStencil/s  MAX', float), ('Time stepper orig name', str), ('OpenMP Threads', int), ('MStencil/s  MAX', float), ('Time unroll',int), ('Sustained Memory BW', float)]
    data = []
    for k in raw_data:
        tup = {}
        # add the general fileds
        for f in req_fields:
            tup[f[0]] = map(f[1], [k[f[0]]] )[0]
    
        # add the stencil operato
        tup['stencil'] = get_stencil_num(k)

        # add the precision information
        if k['Precision'] in 'DP':
            p = 1
        else:
            p = 0
        tup['Precision'] = p
        data.append(tup)

 


    data = sorted(data, key=itemgetter('Time stepper orig name', 'Time unroll', 'Thread group size', 'OpenMP Threads')) 
    #for i in data: print i

 
    max_single = 0
    fig, ax1 = plt.subplots()
    lns = []

    col = 'g'
    ts2 = 'Spatial blk.'
    x = []
    y = []
    y_m = []
    for k in data:
        if ( ('Naive' in k['Time stepper orig name']) and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp)):
            if k['OpenMP Threads'] == 1 and max_single < k['MStencil/s  MAX']/10**3: max_single = k['MStencil/s  MAX']/10**3
            y_m.append(k['Sustained Memory BW']/10**3)
            x.append(k['OpenMP Threads'])
            y.append(k['MStencil/s  MAX']/10**3)
    if(x):
        lns = lns + ax1.plot(x, y, color=col, marker='o', linestyle='-', label='Perf. (GLUP/s)')
    if(y_m):
        ax2 = ax1.twinx()
        lns = lns + ax2.plot(x, y_m, color='r', marker='^', linestyle='-', label='BW (GBytes/s)')


    # add ideal scaling
    ideal = [i*max_single for i in th]
    lns2 =  ax1.plot(th, ideal, color='g', linestyle='--', label='Ideal scaling')


    # add limits
    mem_limit=0
#    sus_mem_bw = 36500 #SB
    sus_mem_bw = 40.0 #IB

    if stencil_kernel == 0:
        mem_limit = sus_mem_bw/(4*4)
    elif stencil_kernel == 1:
        mem_limit = sus_mem_bw/(4*3)
    elif stencil_kernel == 4:
        mem_limit = sus_mem_bw/(4*(3+13))
    elif stencil_kernel == 5:
        mem_limit = sus_mem_bw/(4*(3+7))
    if is_dp == 1: mem_limit = mem_limit / 2

    #if stencil_kernel ==1: 
    lns3 = ax1.plot([5, len(th)], [mem_limit, mem_limit], color='k', linestyle='-', label='Perf. limit')
#        lns = lns + ax2.plot([6, len(th)], [sus_mem_bw, sus_mem_bw], color='m', linestyle='-', label='STREAM BW')


    #title = 'Thread scaling of'
    title = ''
    if stencil_kernel == 0:
        title = '25_pt_const_naive_large'
#        title = title + ' 25_pt constant coeff. star stencil'
    elif stencil_kernel == 1:
        title = '7_pt_const_naive_large'
#        title = title + ' 7_pt constant coeff. star stencil'
    elif stencil_kernel == 2:
        title = title + ' 7_pt varialbe. coeff. star stencil'
    elif stencil_kernel == 3:
        title = title + ' 7_pt varialbe. coeff. axis symm star stencil'
    elif stencil_kernel == 4:
        title = '25_pt_var_naive_large'
#        title = title + ' 25_pt varialbe. coeff. axis symm star stencil'
    elif stencil_kernel == 5:
        title = '7_pt_var_naive_large'
#        title = title + ' 7_pt varialbe. coeff. no symm star stencil'

    #if is_dp == 1:
    #    title = title + ' in double precision'
    #else:
    #    title = title + ' in single precision'

    f_name = title.replace(' ', '_') 

    ax1.set_xlabel('Threads')
    plt.xticks(range(2,11,2))

    if stencil_kernel == 1: 
        ax1.set_ylabel('GLUP/s', color='g', labelpad=0)
        ax2.set_ylabel('GBytes/s', color='r', labelpad=0)

 #   plt.title(title)
 
    lns = lns2 + lns3 + lns
    labs = [l.get_label() for l in lns]
    if stencil_kernel == 1: 
        ax1.legend(lns, labs, loc='lower right')

    ax2.set_ylim(ymin=0, ymax=42)

    ax1.grid()
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
