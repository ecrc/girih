#!/usr/bin/env python
def main():
    import sys

    raw_data = load_csv(sys.argv[1])
    create_table(raw_data)


def create_table(raw_data):
    from operator import itemgetter
    import matplotlib.pyplot as plt
    import pylab
    from csv import DictWriter

    ts_l = set()
    for k in raw_data:
        ts_l.add(k['Time stepper orig name'])
    ts_l = list(ts_l)

    #tb_l = [3, 7]
    tb_l = set()
    for k in raw_data:
        tb_l.add(k['Time unroll'])
    tb_l = list(tb_l)
    tb_l = map(int,tb_l)
    tb_l.sort()
    #print tb_l

    req_fields = [('Time stepper orig name', 0), ('Stencil Kernel semi-bandwidth', 1), ('Stencil Kernel coefficients', 0), ('Precision', 0), ('Time unroll',1), ('Number of time steps',1), ('Number of tests',1), ('Local NX',1), ('Local NY',1), ('Local NZ',1), ('Total Memory Transfer', 2)]
    data = []
    for k in raw_data:
        tup = dict()
        # add the general fileds
        for f in req_fields:
            v = k[f[0]]
            if f[1]==1: v = int(k[f[0]]) 
            if f[1]==2: v = float(k[f[0]]) 
            tup[f[0]] = v

        # add the stencil operator
        if  int(k['Stencil Kernel semi-bandwidth'])==4:
            stencil = '25-pt constant'
        elif  k['Stencil Kernel coefficients'] in 'constant':
            stencil = '7-pt constant'
        else:
            stencil = '7-pt variable'
        tup['Kernel'] = stencil

        data.append(tup)

#    data = sorted(data, key=itemgetter(0, 1, 2, 3,4))
#    for i in data: print i

    data2 = []
    for tup in data:
        tup['Actual Bytes/LUP'] = actual_BpU(tup)
        tup = models(tup)
        data2.append(tup)

    #for i in data2: print i
    from operator import itemgetter
    sorted(data2, key=itemgetter('Time stepper orig name', 'Kernel', 'Precision', 'Time unroll'))

    fields = ['Time stepper orig name'] + ['Kernel'] + ['Precision'] + ['Time unroll'] + ['Actual Bytes/LUP'] + ['No TS'] + ['Interior TS'] + ['All TS'] + ['Interior TS 2']

    with open('Arithmetic_intensity_model.csv', 'w') as output_file:
        r = DictWriter(output_file,fieldnames=fields)
        r.writeheader()
        for k in data2:
            k2 = dict()
            for f in k.keys():
                for f2 in fields:
                    if f == f2:
                        k2[f] = k[f]
            r.writerow(k2)

   

def actual_BpU(tup):
    total_mem = tup['Total Memory Transfer']
    R = tup['Stencil Kernel semi-bandwidth']
    nt = tup['Number of time steps'] * tup['Number of tests']

    nx = tup['Local NX']
    ny = tup['Local NY']
    nz = tup['Local NZ']

    stencil_size = 2*ny*nz + ny*nz*(nx+2*R)

    BpU = (total_mem * 1000**3) / ( stencil_size * nt)
    return BpU


def models(tup):
    if tup['Precision'] == 'DP':
        word_size = 8
    elif tup['Precision'] == 'SP':
        word_size = 4

    R = tup['Stencil Kernel semi-bandwidth']
    TB = tup['Time unroll']
    ny = tup['Local NY']

    if  tup['Stencil Kernel semi-bandwidth']==4:
        t_order = 2
    else:
        t_order = 1
  
    if tup['Time stepper orig name'] == 'Naive':
        if tup['Stencil Kernel coefficients'] in 'constant':
            BpU = (3 + t_order-1) * word_size
        elif tup['Stencil Kernel coefficients'] in 'variable':
            BpU = (3 + t_order-1 + R+1) * word_size
        tup['No TS'] = BpU
        return tup

    width = (TB+1)*2*R
    height = TB*2 + 1
    YT_section = float((TB+1)**2 * 2 * R)
    n_diams = ny/width

    if tup['Stencil Kernel coefficients'] in 'constant':
        BpU1 = (YT_section + width + 2*R + (t_order-1)*width) * word_size / YT_section
        tup['No TS'] = BpU1
       
        BpU2 = (width + 2*R + (t_order-1)*width) * word_size / YT_section
        tup['All TS'] = BpU2

        BpU3 = ((width - 2*R) + 2*width + 2*R + (t_order-1)*width) * word_size / YT_section
        tup['Interior TS'] = BpU3

        BpU4 = ( ((width - 2*R) + width) + ((t_order+1)*width + 2*R) ) * word_size / YT_section
        tup['Interior TS 2'] = BpU4

    elif tup['Stencil Kernel coefficients'] in 'variable':
        BpU1 = (YT_section + width + 2*R + width*(R+1) + (t_order-1)*width) * word_size / YT_section
        tup['No TS'] = BpU1

        BpU2 = (width + 2*R + width*(R+1) + (t_order-1)*width) * word_size / YT_section
        tup['All TS'] = BpU2

        BpU3 = ((width - 2*R) + 2*width + 2*R + width*(R+1) + (t_order-1)*width) * word_size / YT_section
        tup['Interior TS'] = BpU3

        BpU4 = ( ((width - 2*R) + width) + ((t_order+1)*width + 2*R) + ((R+1)*width) ) * word_size / YT_section
        tup['Interior TS 2'] = BpU4


#    print TB, R, width, YT_section, word_size, BpU1, BpU2, BpU3
    return tup


def load_csv(data_file):
    from csv import DictReader
    with open(data_file, 'rb') as output_file:
        data = DictReader(output_file)
        data = [k for k in data]
    return data
    
    
if __name__ == "__main__":
    main()
