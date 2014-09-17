#!/usr/bin/env python
def main():
    import sys

    plot_title='Strong scaling'

    raw_data = load_csv(sys.argv[1])
    plot_stacked_clustered_bars(raw_data, plot_title)
    plot_lines(raw_data, plot_title)

    
def plot_stacked_clustered_bars(raw_data, plot_title):
    from operator import itemgetter
    import matplotlib.pyplot as plt
    import pylab
    import numpy as np
    
    kernels = set()
    for k in raw_data:
        kernels.add(k['Time stepper name'])
    n_ts = len(kernels)
    
    # 1:int, 2:float, 0:str
    req_fields = [('Time stepper name', 0), ('MPI size', 1),
                ('MEAN Computation', 2), ('MEAN Communication', 2),
                ('MAX Computation', 2), ('MAX Communication', 2),
                ('MIN Computation', 2), ('MIN Communication', 2),
                ('Time stepper orig name',0), ('Time unroll', 1)]


    data = []
    for k in raw_data:
        # enforce ploting order and rename kernels
        name = k['Time stepper orig name']
        order = 100
        if 'NB-Naive' in name:
            order = 1
            k['Time stepper name'] = k['Time stepper name'][3:]
        elif 'Intra-Diamond' in name[0:15]: # Blocking diamond
            order = 2
            k['Time stepper name'] = 'i' + k['Time stepper name'][6:]
        elif 'Halo-first' in name:
            order = 3
        elif 'NB-Intra-Diamond' in name: # NB diamond
            order = 4
            k['Time stepper name'] = 'NB-i' + k['Time stepper name'][9:]

        tup = []
        for f in req_fields:
            v = k[f[0]]
            if f[1]==1: v = int(k[f[0]]) 
            if f[1]==2: v = float(k[f[0]]) 
            tup = tup + [v]
        tup = tup + [order]
        data.append(tuple(tup))
        
    data = sorted(data, key=itemgetter(1,10,9))
    #for k in data: print k

    #np_list1 = [d[0]+'-'+str(d[1])+'p' for d in data]
    np_list1 = [str(d[1]) for d in data]
    ts_list1 = [d[0] for d in data]
    
    comp_mean1 = [d[2] for d in data]
    comm_mean1 = [d[3] for d in data]
    comp_max1 = [d[4]-d[2] for d in data]
    comm_max1 = [d[5]-d[3] for d in data]
    comp_min1 = [d[2]-d[6] for d in data]
    comm_min1 = [d[3]-d[7] for d in data]
    
    comp_mean = []
    comm_mean = []
    comp_max = []
    comm_max = []
    comp_min = []
    comm_min = []
    np_list = []
    ts_list = []

    # create a sorted list of kernels and process counts
    procs = map(int,list(set(np_list1)))
    procs = sorted(procs)
    kernels = []
    for i in range(n_ts): kernels.append(ts_list1[i])

    i = 0
    ticks = []
    t = 0.5
    for p in range(len(procs)):
        for j in range(n_ts):
            if j == n_ts/2:
                 if p == len(procs)-1:
                     np_list.append(np_list1[-1])
                 else:
                     np_list.append(np_list1[i])
            else:
                 np_list.append('')

            if i < len(data):
                if p == 0:
                    ts_list.append(kernels[j])
                else:
                    ts_list.append('')

                if kernels[j] == ts_list1[i]:
                    comp_mean.append(comp_mean1[i])
                    comm_mean.append(comm_mean1[i])
                    comp_max.append(comp_max1[i])
                    comm_max.append(comm_max1[i])
                    comp_min.append(comp_min1[i])
                    comm_min.append(comm_min1[i])   
                    i = i + 1
                else:
                    comp_mean.append(0)
                    comm_mean.append(0)
                    comp_max.append(0)
                    comm_max.append(0)
                    comp_min.append(0)
                    comm_min.append(0)
            else:
                comp_mean.append(0)
                comm_mean.append(0)
                comp_max.append(0)
                comm_max.append(0)
                comp_min.append(0)
                comm_min.append(0)
                ts_list.append('')
                np_list.append('')
 
            ticks.append(t)
            t = t + 1
        if p != len(procs)-1  and i < len(data):
            comp_mean.append(0); comp_mean.append(0)
            comm_mean.append(0); comm_mean.append(0)
            comp_max.append(0); comp_max.append(0)
            comm_max.append(0); comm_max.append(0)
            comp_min.append(0); comp_min.append(0)
            comm_min.append(0); comm_min.append(0)
            np_list.append(''); np_list.append('')
            ts_list.append(''); ts_list.append('')
            t = t + 2

    xtk = np.arange(len(comp_mean))    

    fig = plt.figure()
    
    fig.subplots_adjust(bottom=0.25)  # space on the bottom for second time axis

    host = fig.add_subplot(111)  # setup plot

        
    width =1.0
    p1 = host.bar(xtk, comp_mean,   width, color='r', 
                 yerr=[comp_min, comp_max], ecolor='k')
    p2 = host.bar(xtk, comm_mean, width, color='y',
             bottom=comp_mean, yerr=[comm_min, comm_max], ecolor='k')

    host.set_xticks(xtk+width/2)
    host.set_xticklabels(np_list)

    # calculate the height of the plot
    bar_h = []
    for d in range(len(comp_mean)):
        bar_h.append(comp_mean[d] + comm_mean[d] + comm_max[d])
    m_bar = max(bar_h)*1.05

    host.axis((0.0, float(len(np_list)), 0.0, m_bar))

    host.set_ylabel('Time (sec)')
    host.set_xlabel('Procs #')
    plt.title(plot_title)
    plt.legend( (p2[0], p1[0]), ('Communication','Computation'), loc='best')
    plt.grid() 
   

    # Insert the time steppers names at the X-axis 
    newax = host.twiny()  # create new axis
    newax.set_frame_on(False)
    newax.patch.set_visible(False)
    newax.xaxis.set_ticks_position('bottom')
    newax.xaxis.set_label_position('bottom')
    newax.spines['bottom'].set_position(('outward', 15))
#    newax.tick_params('x', width=0)
    newax.set_xticks(ticks)
    newax.set_xticklabels(ts_list, rotation=90, size='small')
    
    newax.axis((0.0, float(len(np_list)), 0.0, m_bar))

    pylab.savefig('scaling_bars.png')

    plt.show() 

        
def plot_lines(raw_data, plot_title):
    from operator import itemgetter
    import matplotlib.pyplot as plt
    import pylab

    kernels = set()
    for k in raw_data:
        kernels.add(k['Time stepper name'])
    kernels = list(kernels)
    
    req_fields = [('Time stepper name', 0), ('MPI size', 1), ('MAX Total', 2)]
    data = []
    for k in raw_data:
        tup = []
        for f in req_fields:
            v = k[f[0]]
            if f[1]==1: v = int(k[f[0]]) 
            if f[1]==2: v = float(k[f[0]]) 
            tup = tup + [v]
        data.append(tuple(tup))
        
    data = sorted(data, key=itemgetter(0,1))
    #for i in data: print i

    for kernel in kernels:
        x = []
        y = []
        for k in data:
            if k[0] == kernel:
                x.append(k[1])
                y.append(k[2])
        plt.plot(x, y, marker='o', linestyle='-', label=kernel)
    plt.xlabel('Procs #')
    plt.ylabel('time (sec)')
    plt.title(plot_title)
    plt.legend(loc='best')
    plt.grid()
    pylab.savefig('scaling_lines.png')
    #plt.show()    

        
def load_csv(data_file):
    from csv import DictReader
    with open(data_file, 'rb') as output_file:
        data = DictReader(output_file)
        data = [k for k in data]
    return data
    
    
if __name__ == "__main__":
    main()
