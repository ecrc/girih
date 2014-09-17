#!/usr/bin/env python
def main():
    import sys

    raw_data = load_csv(sys.argv[1])
    create_table(raw_data)

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

    req_fields = [('WD main-loop RANK0 MStencil/s  MAX', 2), ('Time stepper orig name', 0), ('Stencil Kernel semi-bandwidth', 1), ('Stencil Kernel coefficients', 0), ('Precision', 0), ('Time unroll',1), ('Number of time steps',1), ('Number of tests',1), ('Local NX',1), ('Local NY',1), ('Local NZ',1), ('Total Memory Transfer', 2), ('Thread group size' ,1), ('Intra-diamond prologue/epilogue MStencils',1), ('Total cache block size (kB):',1)]
    data = []
    for k in raw_data:
        tup = dict()
        # defaults
        if k['Intra-diamond prologue/epilogue MStencils'] == '':
            k['Intra-diamond prologue/epilogue MStencils'] = 0
        if k['Total cache block size (kB):'] == '':
            k['Total cache block size (kB):'] = 0
        # add the general fileds
        for f in req_fields:
            try:
                v = k[f[0]]
                if f[1]==1: v = int(k[f[0]]) 
                if f[1]==2: v = float(k[f[0]]) 
            except:
                print f[0]
                    
            tup[f[0]] = v

        # add the stencil operator
        tup['Kernel'] = get_stencil_num(k)
        data.append(tup)

#    data = sorted(data, key=itemgetter(0, 1, 2, 3,4))
#    for i in data: print i

    data2 = []
    for tup in data:
        if tup['Local NX'] > 96:
            tup['Actual Bytes/LUP'] = actual_BpU(tup)
            tup['Model'] = models(tup)
            # model error
            tup['Err %'] = 100 * (tup['Model'] - tup['Actual Bytes/LUP'])/tup['Actual Bytes/LUP']
            tup['D_width'] = (tup['Time unroll']+1)*2*tup['Stencil Kernel semi-bandwidth']
            tup['Performance'] = tup['WD main-loop RANK0 MStencil/s  MAX']
            data2.append(tup)

    #for i in data2: print i
    from operator import itemgetter
    data2 = sorted(data2, key=itemgetter('Time stepper orig name', 'Kernel', 'Thread group size', 'Local NX', 'D_width'))

    fields = ['Time stepper orig name', 'Kernel', 'Thread group size', 'Local NX', 'Precision', 'D_width', 'Total cache block size (kB):', 'Actual Bytes/LUP', 'Model', 'Err %', 'Performance']

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

    oh = tup['Intra-diamond prologue/epilogue MStencils']

    stencil_size = 2*ny*nz + ny*nz*(nx+2*R) 
    BpU = (total_mem * 10**9) / ( stencil_size * nt - oh*10**6*tup['Number of tests'])
    return BpU


def models(tup):
    if   tup['Precision'] == 'DP': word_size = 8
    elif tup['Precision'] == 'SP': word_size = 4

    R = tup['Stencil Kernel semi-bandwidth']
    TB = tup['Time unroll']
    ny = tup['Local NY']

    # number of streamed copies of the domain (buffers)
    if   tup['Kernel'] == 0: nb = 3
    elif tup['Kernel'] == 1: nb = 2
    elif tup['Kernel'] == 4: nb = 2+13
    elif tup['Kernel'] == 5: nb = 2+7



    width = (TB+1)*2*R
    YT_section = float((TB+1)**2 * 2 * R)

    # no temporal blocking model
    if tup['Time stepper orig name'] == 'Naive':
        bpu = (1 + nb) * word_size
    else: # temporal blocking model
        bpu = ( ((width - 2*R) + width) + (nb*width + 2*R) ) * word_size / YT_section
   
    return bpu


def load_csv(data_file):
    from csv import DictReader
    with open(data_file, 'rb') as output_file:
        data = DictReader(output_file)
        data = [k for k in data]
    return data
    
    
if __name__ == "__main__":
    main()

#    if 'constant' in tup['Stencil Kernel coefficients']:
#        BpU1 = (YT_section + width + 2*R + (t_order-1)*width) * word_size / YT_section
#        tup['No TS'] = BpU1
#       
#        BpU2 = (width + 2*R + (t_order-1)*width) * word_size / YT_section
#        tup['All TS'] = BpU2
#
#        BpU3 = ((width - 2*R) + 2*width + 2*R + (t_order-1)*width) * word_size / YT_section
#        tup['Interior TS'] = BpU3
#
#        BpU4 = ( ((width - 2*R) + width) + ((t_order+1)*width + 2*R) ) * word_size / YT_section
#        tup['Interior TS 2'] = BpU4
#
#    elif 'variable' in  tup['Stencil Kernel coefficients']:
#        BpU1 = (YT_section + width + 2*R + width*(6*R+1) + (t_order-1)*width) * word_size / YT_section
#        tup['No TS'] = BpU1
#
#        BpU2 = (width + 2*R + width*(6*R+1) + (t_order-1)*width) * word_size / YT_section
#        tup['All TS'] = BpU2
#
#        BpU3 = ((width - 2*R) + 2*width + 2*R + width*(6*R+1) + (t_order-1)*width) * word_size / YT_section
#        tup['Interior TS'] = BpU3
#
#        BpU4 = ( ((width - 2*R) + width) + ((t_order+1)*width + 2*R) + ((R*6+1)*width) ) * word_size / YT_section
#        tup['Interior TS 2'] = BpU4


