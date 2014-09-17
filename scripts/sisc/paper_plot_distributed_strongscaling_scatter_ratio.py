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

   # matplotlib.rcParams.update({'font.size': 18})


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

    req_fields = [('OpenMP Threads', int), ('Thread group size', int), ('WD main-loop RANK0 MStencil/s  MAX', float), ('Time stepper orig name', str), ('MPI size', int), ('MStencil/s  MAX', float)]
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

        # add MWD profiling
        measures = ['Wavefront barrier wait [s]:', 'Wavefront barrier wait [%]:',
                   'Wavefront steady state [s]:', 'Wavefront steady state [%]:',
                   'Wavefront startup/end [s]:', 'Wavefront startup/end [%]:',
                   'Wavefront communication [s]:', 'Wavefront communication [%]:',
                   'Wavefront others [s]:', 'Wavefront others [%]:',
                   'Group spin-wait [s]:', 'Group spin-wait [%]:',
                   'Resolved diamonds:']
        if tup['Time stepper orig name'] == 'Dynamic-Intra-Diamond':
            for i in range(tup['OpenMP Threads']):
                for m in measures:
                    key =m[:-1] + ' group %d'%i
                    if key in k:
                        if k[key] != '':
                            tup[key] = float(k[key]) if '.' in k[key] else int(k[key])

        measures =['RANK0 Total',
                    'RANK0 Computation',
                    'RANK0 Communication',
                    'RANK0 Waiting',
                    'RANK0 Other']
        if 'Naive' in tup['Time stepper orig name']:
            for m in measures:
                tup[m] = float(k[m])
            tup['RANK0 Computation [%]'] = tup['RANK0 Computation']/tup['RANK0 Total'] * 100
        data.append(tup)

    data = sorted(data, key=itemgetter('Time stepper orig name', 'Thread group size', 'MPI size'))
    #for i in data: print i


    linestyle = ''
    markersize = 5 
    x = []
    y = []
    max_x = 0
    for k in data:
        if ( ('Naive' in k['Time stepper orig name']) and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp)):
            x.append(k['MPI size']+0.2)
            y.append(k['RANK0 Computation [%]'])
    marker = 'o'
    col = 'g'
    ts2 = 'Spt.blk.'
    if(x):
        plt.plot(x, y, color=col, marker=marker, markersize=markersize, linestyle=linestyle, label=ts2)
        max_y = max(y)
        max_x = max(x)
 

    cols = {0:'y', 1:'k', 2:'b', 5:'r', 10:'m'}
    shift = {0:0, 1:0, 2:-.2, 5:-.1, 10:.1}
    markers = {0:'.', 1:'^', 2:'v', 4:'.', 5:'x', 10:'*'}
    for tgs in [1,2,5, 10]:
        marker = markers[tgs]
        x = []
        y = []
        y_m = []
        for k in data:
            if ( ('Diamond' in k['Time stepper orig name']) and (k['Thread group size'] == tgs)  and (k['stencil']==stencil_kernel) and (k['Precision']==is_dp) ):
                for i in range(k['OpenMP Threads']/k['Thread group size']):
                    x.append(k['MPI size']+ shift[tgs])
                    y.append(k['Wavefront steady state [%] group '+str(i)]+
                             k['Wavefront startup/end [%] group '+str(i)])
        col = cols[tgs]
        ts2 = str(tgs) + 'WD'
        if(x):
            plt.plot(x, y, color=col, marker=marker, markersize=markersize, linestyle=linestyle, label=ts2)
            max_y = max(max(y), max_y)
            max_x = max(max(x), max_x)

    plt.ylabel('Compute time (%)')
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
        title = '25_pt_var_dist_comp_ratio'
    elif stencil_kernel == 5:
    #    title = title + ' 7_pt varialbe. coeff. no symm star stencil'
        title = '7_pt_var_dist_comp_ratio'

   # if is_dp == 1:
   #     title = title + ' in double precision'
   # else:
   #     title = title + ' in single precision'

    f_name = title.replace(' ', '_')

    plt.ylim([0,max_y*1.1])

    plt.xlabel('Processors #')

    np_l = [k['MPI size'] for k in data if k['stencil']==stencil_kernel]
    plt.xticks(np.concatenate((np.array([1,2]), np.arange(4, max(np_l)+5, 4))))
    #plt.xticks(np.arange(0, max(np_l)+3, 2))
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
