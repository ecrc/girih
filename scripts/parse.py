#!/usr/bin/env python


def sort_fields(data,items):
    for i in items:
        data.remove(i)
    data = items + list(data)
    return data

def main():
    #from itertools import izip
    #from csv import reader, writer
    from csv import DictWriter
    import sys
    import os

    output_name = "summary.csv"
    
    all_fields = set()
    data = []
    for f in sys.argv[1:]:
        with open(f, 'r') as fp:
            k = get_summary(fp)
        k['file_name'] = os.path.split(f)[-1]
        all_fields.update(set(k.keys()))
        if (k.has_key('Time stepper name')) and (k.has_key('MStencil/s  MAX')):
            data.append(k)
            print('Parsed the file: '+ f)
        else:
            print('Rejected the file: '+ f)
            

#    order=['Time stepper name', 'Time unroll', 'npx', 'npy', 'npz', 'Local NX', 'Local NY', 'Local NZ']
    order = ['Stencil Kernel semi-bandwidth', 'Time stepper orig name', 'Thread group size', 'MPI size', 'MStencil/s  MAX','WD main-loop RANK0 MStencil/s  MAX', 'Intra-diamond width']
    fields = sort_fields(all_fields, order) 
    

    with open(output_name, 'w') as output_file:
        r = DictWriter(output_file,fieldnames=fields)
        r.writeheader()
        for k in data:
            r.writerow(k)

# transpose the file
    #with open(output_name, 'rb') as output_file:
        #data = izip(*reader(output_file))
    #with open(output_name, 'w') as output_file:
        #writer(output_file).writerows(data)
    
    
def get_summary(f):
    float_fields =   ('MStencil/s  MIN',
                'MStencil/s  AVG',
                'MStencil/s  MAX',
                'RANK0 Total',
                'RANK0 Computation',
                'RANK0 Communication',
                'RANK0 Waiting',
                'RANK0 Other',
                'MAX Total',
                'MAX Computation',
                'MAX Communication',
                'MAX Waiting',
                'MAX Other',
                'MIN Total',
                'MIN Computation',
                'MIN Communication',
                'MIN Waiting',
                'MIN Other',
                'MEAN Total',
                'MEAN Computation',
                'MEAN Communication',
                'MEAN Waiting',
                'MEAN Other',
                'WD main-loop RANK0 MStencil/s  MAX', 
                'RANK0 ts main loop')

    str_fields = ('Time stepper name',
                'Stencil Kernel name',
                'Stencil Kernel coefficients',
                'Precision')
                
    int_fields = ('Number of time steps',
                'Alignment size',
                'Number of tests',
                'Verify',
                'Time unroll',
                'Intra-diamond width',
                'Halo concatenation:',
                'OpenMP Threads',
                'MPI size',
                'Stencil Kernel semi-bandwidth',
                'number of idiamond wavefronts',
                'Time parallel muti-core wavefront',
                'Using separate call to 1-stride loop',
                'Thread group size',
                'Intra-diamond prologue/epilogue MStencils',
                'Cache block size/wf (kB):',
                'Total cache block size (kB):')

    mlist = []
    # default not cotiguous MPI datatype
    z_contiguous = 0
    mlist.append(('Halo concatenation', 0))
    mlist.append(('Intra-diamond width', 0))

    mlist.append(('Thread group size', 0))
    mlist.append(('WD main-loop RANK0 MStencil/s  MAX', 0))

    for line in f:
        # General float cases
        for field in float_fields:
            if field in line:
                val = float(line.split(':')[1].split()[0])
                mlist.append((field,val))

        # General Integer cases
        for field in int_fields:
            if field in line:
                val = int(line.split(':')[1].split()[0])
                mlist.append((field,val))

        # special case for modified fields
        if 'Intra-diamond halo concatenation' in line:
            val = int(line.split(':')[1].split()[0])
            mlist.append(('Halo concatenation',val))

        # Global size
#Global domain    size:262144    nx:64    ny:64    nz:64
        if 'Global domain' in line:
            dims = line.split()[3:6]
            nx = dims[0].split(':')[1]
            ny = dims[1].split(':')[1]
            nz = dims[2].split(':')[1]
            mlist.append(('Global NX',nx))
            mlist.append(('Global NY',ny))
            mlist.append(('Global NZ',nz))

        # Local size
#Rank 0 domain    size:131072    nx:32    ny:64    nz:64
        if 'Rank 0 domain' in line:
            dims = line.split()[4:7]
            nx = dims[0].split(':')[1]
            ny = dims[1].split(':')[1]
            nz = dims[2].split(':')[1]
            mlist.append(('Local NX',nx))
            mlist.append(('Local NY',ny))
            mlist.append(('Local NZ',nz))

        # topology information
        if 'Processors topology' in line:
            dims = line.split(':')[1].strip()
            npx,npy,npz = dims.split(',')
            mlist.append(('npx',npx))
            mlist.append(('npy',npy))
            mlist.append(('npz',npz))

        # General string cases
        for field in str_fields:
            if field in line:
                val = line.split(':')[1].strip()
                mlist.append((field,val))

        # MWD statistics information
        measures = ['Wavefront barrier wait [s]:', 'Wavefront barrier wait [%]:',
                   'Wavefront steady state [s]:', 'Wavefront steady state [%]:',
                   'Wavefront startup/end [s]:', 'Wavefront startup/end [%]:',
                   'Wavefront communication [s]:', 'Wavefront communication [%]:',
                   'Wavefront others [s]:', 'Wavefront others [%]:',
                   'Group spin-wait [s]:', 'Group spinn-wait [%]:',
                   'Resolved diamonds:']
        for m in measures:
            if m in line:
                vals = [i.strip() for i in line.split(':')[1].split(' ') if len(i.strip()) !=0]
                for i in range(len(vals)):
                    mlist.append((m[:-1] + ' group %d'%i ,vals[i])) 



        # LIKWID performance results
        if 'Measuring group' in line:
            val = line.split(' ')[2].strip()
            mlist.append(('LIKWID performance counter',val))

        if '|    Memory BW [MBytes/s]     |' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Sustained Memory BW',val))

        if '| Memory data volume [GBytes] |' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Total Memory Transfer',val))


        if '| L3 bandwidth [MBytes/s] |' in line:
            vals = [i.strip() for i in line.split('|')[2:-1]]
            if len(vals) == 1:
                mlist.append(('Sustained L3 BW',vals[0]))
            else:
                for i in range(len(vals)):
                    mlist.append(('Sustained L3 BW core %d'%i ,vals[i])) 

        if '| L3 data volume [GBytes] |' in line:
            vals = [i.strip() for i in line.split('|')[2:-1]]
            if len(vals) == 1:
                mlist.append(('Total L3 Transfer',vals[0]))
            else:
                for i in range(len(vals)):
                    mlist.append(('Total L3 Transfer core %d'%i ,vals[i]))
                    
        if '| L3 bandwidth [MBytes/s] STAT |' in line:
            vals = [i.strip() for i in line.split('|')[2:6]]
            mlist.append(('Sustained L3 BW sum',vals[0]))
            mlist.append(('Sustained L3 BW max',vals[1]))
            mlist.append(('Sustained L3 BW min',vals[2]))
            mlist.append(('Sustained L3 BW avg',vals[3]))

        if '| L3 data volume [GBytes] STAT |' in line:
            vals = [i.strip() for i in line.split('|')[2:6]]
            mlist.append(('Total L3 Transfer sum',vals[0]))
            mlist.append(('Total L3 Transfer max',vals[1]))
            mlist.append(('Total L3 Transfer min',vals[2]))
            mlist.append(('Total L3 Transfer avg',vals[3]))

        if '|  L1 DTLB miss rate STAT   |' in line:
            vals = [i.strip() for i in line.split('|')[2:6]]
            mlist.append(('L1 DTLB miss rate sum',vals[0]))
            mlist.append(('L1 DTLB miss rate max',vals[1]))
            mlist.append(('L1 DTLB miss rate min',vals[2]))
            mlist.append(('L1 DTLB miss rate avg',vals[3]))

        if '| L2 bandwidth [MBytes/s] |' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Sustained L2 BW',val))

        if '| L2 data volume [GBytes] |' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Total L2 Transfer',val))

        if '|      Energy [J]      |' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Energy',val))
 
        if '|      Power [W]       |' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Power',val))
 
        if '|   Energy DRAM [J]    |' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Energy DRAM',val))
 
        if '|    Power DRAM [W]    |' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Power DRAM',val))
         
        # kernel name
#        field='Time stepper name'
#        if field in line:
#            val = line.split(':')[1].strip()
#            mlist.append((field,val))
        
        # cotiguous MPI datatype?
        field='contiguous across the Z direction'
        if field in line:
            z_contiguous = 1


    if(z_contiguous == 1):
        mlist.append(('contig-z', 1))
    else:
        mlist.append(('contig-z', 0))

    mlist = dict(mlist)
    
    if mlist.has_key('Time stepper name'):  
        mlist['Time stepper orig name']= mlist['Time stepper name']
        # modify the time stepper name when the MPI datatype is contiguous
        if 'Diamond' not in mlist['Time stepper name'] and mlist['contig-z'] == 1:
            new_ts_name = mlist['Time stepper name'] + '-zcontig'
            mlist['Time stepper name'] = new_ts_name

        if ('Halo-first' == mlist['Time stepper name'] or 'NB-Naive' == mlist['Time stepper name'] or 'Diamond' in mlist['Time stepper name']) and mlist['Halo concatenation'] == 0:
                new_ts_name = mlist['Time stepper name'] + '-no-concat'
                mlist['Time stepper name'] = new_ts_name

        # add the time unrolling in the time stepper name
        if 'Diamond' in mlist['Time stepper name']:
            new_ts_name = mlist['Time stepper name'] + '-' + str(mlist['Time unroll']) + 'ur'
            mlist['Time stepper name'] = new_ts_name


    return mlist

if __name__ == "__main__":
    main()
