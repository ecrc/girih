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

    #matplotlib.rcParams.update({'font.size': 14})
    matplotlib.rcParams.update({'font.size': 16})


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
            if k['OpenMP Threads'] == 1 and max_single < k['MStencil/s  MAX']: max_single = k['MStencil/s  MAX']/1024.0
            y_m.append(k['Sustained Memory BW']/1024.0)
            x.append(k['OpenMP Threads'])
            y.append(k['MStencil/s  MAX']/1024.0)
    marker = 'o'
    col = 'g'
    ts2 = 'Spatial blk.'
    if(x) and t==0:
        plt.plot(x, y, color=col, marker=marker, linestyle='-', label=ts2)
    if(y_m) and t==1:
        plt.plot(x, y_m, color=col, marker=marker, linestyle='-', label=ts2)



    perf_str = 'WD main-loop RANK0 MStencil/s  MAX'
    cols = {0:'y', 1:'k', 2:'b', 4:'c', 5:'c', 8:'m'}
    for tgs in [1,2, 4, 8, 5]:
        x = []
        y = []
        y_m = []
        for k in data:
            if ( ('Diamond' in k['Time stepper orig name']) and (k['Thread group size'] == tgs)  and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp) ):
                if k['OpenMP Threads'] == 1 and max_single < k[perf_str]: max_single = k[perf_str]/1024.0
                y_m.append(k['Sustained Memory BW']/1024.0)
                x.append(k['OpenMP Threads'])
                y.append(k[perf_str]/1024.0)
        marker = 'o'
        col = cols[tgs]
        ts2 = str(tgs) + 'core-WD'
        if(x) and t==0:
            plt.plot(x, y, color=col, marker=marker, linestyle='-', label=ts2)
        if(y_m) and t==1:
            plt.plot(x, y_m, color=col, marker=marker, linestyle='-', label=ts2)

    x = []
    y = []
    y_m = []
    for k in data:
        if ( ('Diamond' in k['Time stepper orig name']) and (k['Thread group size'] == 10) and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp)):
            if k['OpenMP Threads'] == 1 and max_single < k[perf_str]: max_single = k[perf_str]/1024.0
            y_m.append(k['Sustained Memory BW']/1024.0)
            x.append(k['OpenMP Threads'])
            y.append(k[perf_str]/1024.0)
    marker = 'o'
    col = 'm'
    ts2 = str(10) + 'core-WD'
    if(x) and t==0:
        plt.plot(x, y, color=col, marker=marker, linestyle='-', label=ts2)
    if(y_m) and t==1:
        plt.plot(x, y_m, color=col, marker=marker, linestyle='-', label=ts2)




    # add limits
    mem_limit=0
#    sus_mem_bw = 36500 #SB
    sus_mem_bw = 40000 #IB

    if stencil_kernel == 0:
        mem_limit = sus_mem_bw/16 /1024.0
    elif stencil_kernel == 1:
        mem_limit = sus_mem_bw/12 /1024.0
    elif stencil_kernel == 2:
        mem_limit = sus_mem_bw/20 / 1024.0
    if is_dp == 1: mem_limit = mem_limit / 2
    if t == 0:
        #plt.plot([1, len(th)], [mem_limit, mem_limit], color='g', linestyle='--', label='Spatial blk. limit')
        pass

    # add ideal scaling
    ideal = [i*max_single for i in th]
    if t == 0:
        plt.plot(th, ideal, color='k', linestyle='--', label='Ideal scaling')

    if t == 0:
        title = 'Thread scaling performance'
        plt.ylabel('GStencil update/s')
    else:
        title = 'Thread scaling main memory BW usage'
        plt.ylabel('Transferred GBytes/s')

    plt.title(title)

    if stencil_kernel == 0:
        title = title + ' 25_pt constant coeff. star stencil'
    elif stencil_kernel == 1:
        title = title + ' 7_pt constant coeff. star stencil'
    elif stencil_kernel == 2:
        title = title + ' 7_pt varialbe. coeff. star stencil'
    elif stencil_kernel == 3:
        title = title + ' 7_pt varialbe. coeff. axis symm star stencil'
    elif stencil_kernel == 4:
        title = title + ' 25_pt varialbe. coeff. axis symm star stencil'
    elif stencil_kernel == 5:
        title = title + ' 7_pt varialbe. coeff. no symm star stencil'




    if is_dp == 1:
        title = title + ' in double precision'
    else:
        title = title + ' in single precision'

    f_name = title.replace(' ', '_')

    plt.xlabel('Cores #')
 
    plt.legend(loc='best')
    plt.grid()
    pylab.savefig(f_name+'.png')
    pylab.savefig(f_name+'.pdf', format='pdf')

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
