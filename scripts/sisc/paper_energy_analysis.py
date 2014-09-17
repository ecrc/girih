#!/usr/bin/env python
def main():
    import sys

    raw_data = load_csv(sys.argv[1])
    create_table(raw_data)

def get_stencil_num(k):
    # add the stencil operator
    if  k['Stencil Kernel coefficients'] in 'constant':
        if  int(k['Stencil Kernel semi-bandwidth'])==4:
            stencil = '25pt_const'
        else:
            stencil = '7pt_const'
    elif  'no-symmetry' in k['Stencil Kernel coefficients']:
        stencil = '7pt_var'
    elif  'sym' in k['Stencil Kernel coefficients']:
        if int(k['Stencil Kernel semi-bandwidth'])==1:
            stencil = '7pt_var_ax_sym'
        else:
            stencil = '25pt_var'
    else:
        stencil = '7pt_var_all_sym'
    return stencil

def create_table(raw_data):
    from operator import itemgetter
    import matplotlib.pyplot as plt
    import pylab
    from csv import DictWriter

    req_fields = [('Time stepper orig name', 0), ('Stencil Kernel semi-bandwidth', 1), ('Stencil Kernel coefficients', 0), ('Precision', 0), ('Number of time steps',1), ('Number of tests',1), ('Global NX',1), ('Global NY',1), ('Global NZ',1), ('Thread group size' ,1), ('Intra-diamond prologue/epilogue MStencils',1), ('Energy', 2), ('Energy DRAM', 2), ('Power',2), ('Power DRAM',2), ('WD main-loop RANK0 MStencil/s  MAX', 2),('MStencil/s  MAX', 2), ('OpenMP Threads',1)]
    data = []
    for k in raw_data:
        tup = dict()
        # defaults
        if k['Intra-diamond prologue/epilogue MStencils'] == '':
            k['Intra-diamond prologue/epilogue MStencils'] = 0
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
        tup['Stencil'] = get_stencil_num(k)
        data.append(tup)

#    for i in data: print i

    data2 = []
    for tup in data:
        glups = (tup['Number of time steps'] * tup['Global NX']*tup['Global NY']*tup['Global NZ'] - tup['Intra-diamond prologue/epilogue MStencils']*10**6 ) * tup['Number of tests'] / 10**9
        tup['Total pJoul/LUP'] = (tup['Energy'] + tup['Energy DRAM'])/glups
        tup['DRAM pJoul/LUP'] = (tup['Energy DRAM'])/glups
        tup['CPU pJoul/LUP'] = (tup['Energy'])/glups
        if 'Dynamic' in tup['Time stepper orig name']:
            tup['Time stepper orig name'] = 'MWD'

        if 'Dynamic' in tup['Time stepper orig name']:
            tup['Performance'] = tup['WD main-loop RANK0 MStencil/s  MAX']
        else:
            tup['Performance'] = tup['MStencil/s  MAX']
        tup['Threads'] = tup['OpenMP Threads']
        tup['Method'] = tup['Time stepper orig name']
        data2.append(tup)


    #for i in data2: print i
    from operator import itemgetter
    data2 = sorted(data2, key=itemgetter('Stencil', 'Thread group size', 'Time stepper orig name', 'Global NX', 'Global NY','Global NZ'))

    fields = ['Method', 'Stencil', 'Threads', 'Thread group size', 'Global NX', 'Global NY','Global NZ', 'Precision', 'Power', 'Power DRAM', 'CPU pJoul/LUP', 'DRAM pJoul/LUP', 'Total pJoul/LUP', 'Performance']

    with open('energy_consumption.csv', 'w') as output_file:
        r = DictWriter(output_file,fieldnames=fields)
        r.writeheader()
        for k in data2:
            k2 = dict()
            for f in k.keys():
                for f2 in fields:
                    if f == f2:
                        k2[f] = k[f]
            r.writerow(k2)


def load_csv(data_file):
    from csv import DictReader
    with open(data_file, 'rb') as output_file:
        data = DictReader(output_file)
        data = [k for k in data]
    return data
    
    
if __name__ == "__main__":
    main()

