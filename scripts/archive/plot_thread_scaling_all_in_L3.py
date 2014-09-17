#!/usr/bin/env python
def main():
    import sys

    raw_data = load_csv(sys.argv[1])

 #   for ts in ['Naive', 'Dynamic-Intra-Diamond']  
    #for k in [0, 1, 2]:
    for k in [1]:
        #for is_dp in [0,1]:
        for is_dp in [0]:
            plot_lines(raw_data, k, is_dp)

        
def plot_lines(raw_data, stencil_kernel, is_dp):
    from operator import itemgetter
    import matplotlib.pyplot as plt
    import pylab

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
    
    req_fields = [('Time stepper orig name', 0), ('OpenMP Threads', 1), ('MStencil/s  MAX', 2), ('Time unroll',1), ('Thread group size',1)]
    data = []
    for k in raw_data:
        tup = {}
        # add the general fileds
        for f in req_fields:
            v = k[f[0]]
            if f[1]==1: v = int(k[f[0]]) 
            if f[1]==2: v = float(k[f[0]]) 
            tup[f[0]] = v

        # add the stencil operator
        if  int(k['Stencil Kernel semi-bandwidth'])==4:
            stencil = 0
        elif  k['Stencil Kernel coefficients'] in 'constant':
            stencil = 1
        else:
            stencil = 2
        tup['kernel'] = stencil

        # add the precision information
        if k['Precision'] in 'DP':
            tup['Precision'] = 1
        else:
            tup['Precision'] = 0

        data.append(tup)

    
    data = sorted(data, key=itemgetter('Time stepper orig name', 'Thread group size', 'Time unroll', 'OpenMP Threads'))
    #for i in data: print i

    max_single = 0
   
    # Plot naive results
    x = []
    y = []
    for k in data:
        if (k['Time stepper orig name'] == 'Naive') and (k['kernel']==1) and (k['Precision']==0):
            if k['OpenMP Threads'] == 1 and max_single < k['MStencil/s  MAX']: max_single = k['MStencil/s  MAX']
            x.append(k['OpenMP Threads'])
            y.append(k['MStencil/s  MAX'])
    plt.plot(x, y, color='g', marker='o', linestyle='-', label='Naive')
   
    # Plot all core WD results
    x = []
    y = []
    for k in data:
        if (k['Time stepper orig name'] == 'Intra-Diamond') and (k['kernel']==1) and (k['Precision']==0):
            if k['OpenMP Threads'] == 1 and max_single < k['MStencil/s  MAX']: max_single = k['MStencil/s  MAX']
            x.append(k['OpenMP Threads'])
            y.append(k['MStencil/s  MAX'])
    plt.plot(x, y, color='r', marker='o', linestyle='-', label='All-core-WD')

    # Plot MWD results
    for tgs, col in [(1, 'b'), (2, 'c'), (5, 'm')]:
        x = []
        y = []
        for k in data:
            if (k['Time stepper orig name'] == 'Dynamic-Intra-Diamond') and (k['Thread group size'] == tgs) and (k['kernel']==1) and (k['Precision']==0):
                if k['OpenMP Threads'] == 1 and max_single < k['MStencil/s  MAX']: max_single = k['MStencil/s  MAX']
                x.append(k['OpenMP Threads'])
                y.append(k['MStencil/s  MAX'])
        plt.plot(x, y, color=col, marker='o', linestyle='-', label=str(tgs)+'-core-WD')
   
 
    mem_limit=0
    sus_mem_bw = 36500
    if stencil_kernel == 0:
        mem_limit = sus_mem_bw/16
    elif stencil_kernel == 1:
        mem_limit = sus_mem_bw/12
    elif stencil_kernel == 2:
        mem_limit = sus_mem_bw/20
    if is_dp == 1: mem_limit = mem_limit / 2
    #plt.plot([1,8], [mem_limit, mem_limit], color='g', linestyle='--', label='Naive mem BW limit')


    # add ideal scaling
    ideal = [i*max_single for i in th]
    plt.plot(th, ideal, color='r', linestyle='--', label='Ideal scaling')

    title = 'Thread scaling of'
    if stencil_kernel == 0:
        title = title + ' 25_pt constant coeff. star stencil'
    elif stencil_kernel == 1:
        title = title + ' 7_pt constant coeff. star stencil'
    elif stencil_kernel == 2:
        title = title + ' 7_pt varialbe. coeff. star stencil'

    if is_dp == 1:
        title = title + ' in double precision'
    else:
        title = title + ' in single precision'

    f_name = title.replace(' ', '_') 

    plt.xlabel('Threads #')
    plt.ylabel('MStencils/s')
    plt.title(title)
    plt.legend(loc='best')
    plt.grid()
    pylab.savefig(f_name+'.png')
    pylab.savefig(title+'.pdf', format='pdf')

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
