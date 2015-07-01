#!/usr/bin/env python


def sort_cols(data,items):
    items2 = []
    for i in items:
        if i in data:
            data.remove(i)
            items2.append(i)
    data = items2 + list(data)
    return data

def main():
    #from itertools import izip
    #from csv import reader, writer
    from operator import itemgetter
    from csv import DictWriter
    import sys, os, re, subprocess

    output_name = "summary.csv"

    all_fields = set()
    data = []
    for f in sys.argv[1:]:
        with open(f, 'r') as fp:
            k = get_summary(fp)
        k['file_name'] = os.path.split(f)[-1]
        all_fields.update(set(k.keys()))
        has_res=0
        for key in k.keys():
          if re.search('MStencil/s *MAX', key): has_res=1
        if (k.has_key('Time stepper name')) and has_res==1:
            data.append(k)
            print('Parsed the file: '+ f)
        else:
#            for i,v in k.iteritems(): print i,' :',v
            print('Rejected the file: '+ f)
#            subprocess.call('rm '+ f, shell=True)


    cols_order = ['Thread group size', 'Wavefront parallel strategy', 'LIKWID performance counter', 'Global NX']

    cols_order.append('Sustained Memory BW')
    cols_order.append('Total Memory Transfer')
    for i, st in enumerate(['sum', 'max', 'min', 'avg']):
        cols_order.append('L1 DTLB miss rate %s'%(st))
    snames = [('| %s bandwidth [MBytes/s] ','BW'),
              ('| %s data volume [GBytes] ', 'data volume'),
              ('|   %s Evict [MBytes/s] ', 'evict'),
              ('|   %s Load [MBytes/s] ', 'load')]
    for cache in ['L2', 'L3']:
        for sn in snames:
            for i, st in enumerate(['sum', 'max', 'min', 'avg']):
                cols_order.append(cache+' '+sn[1]+' '+st)
    for cache in ['L2', 'L3']:
        for sn in snames:
            for i in range(1000):
                cols_order.append(cache+' '+sn[1]+' c'+str(i))
    snames = [('|         CPI ', 'CPI'),
              ('| Load to Store ratio ', 'Load to Store ratio')]
    for sn in snames:
        for i, st in enumerate(['sum', 'max', 'min', 'avg']):
            cols_order.append(sn[1]+' '+st)
    for sn in snames:
        for i in range(1000):
            cols_order.append(sn[1]+' c'+str(i))

    for i in range(1000):
        cols_order.append('Wavefront barrier wait [%] group '+str(i))


    try:
      fields = sort_cols(all_fields, cols_order)
      data = sorted(data, key=itemgetter('Stencil Kernel semi-bandwidth', 'Stencil Kernel coefficients', 'Time stepper orig name', 'Thread group size', 'Wavefront parallel strategy', 'LIKWID performance counter', 'OpenMP Threads', 'MPI size', 'Global NX'))
    except:
      pass


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
    import re

    float_fields =   (
                'RANK0 MStencil/s  MIN',
                'RANK0 MStencil/s  AVG',
                'RANK0 MStencil/s  MAX',
                'MWD main-loop RANK0 MStencil/s MIN',
                'MWD main-loop RANK0 MStencil/s MAX',
#                'MWD main-loop RANK0 MStencil/s  MAX', # temporary for depricated format
                'Total RANK0 MStencil/s MIN',
                'Total RANK0 MStencil/s MAX',
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
                'RANK0 ts main loop')

    str_fields = ('Time stepper name',
                'Stencil Kernel name',
                'Stencil Kernel coefficients',
                'Precision',
                'Wavefront parallel strategy')

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
                'Multi-wavefront updates',
                'Time parallel muti-core wavefront',
                'Using separate call to central line update',
                'Thread group size',
                'Intra-diamond prologue/epilogue MStencils',
                'Cache block size/wf (kiB)',
                'Block size in X',
                'Total cache block size (kiB)',
                'Threads per core',
                'Threads along x-axis',
                'Threads along y-axis',
                'User set thread group size',
                'Threads along z-axis')

    mlist = []
    # default not cotiguous MPI datatype
    z_contiguous = 0
    mlist.append(('Halo concatenation', 0))
    mlist.append(('Intra-diamond width', 0))

    mlist.append(('Cache block size/wf (kiB)', 0))
    mlist.append(('Total cache block size (kiB)', 0))
    mlist.append(('Thread group size', 0))
    mlist.append(('Multi-wavefront updates', 0))
    mlist.append(('Intra-diamond prologue/epilogue MStencils', 0))
    mlist.append(('Wavefront parallel strategy',0))
    mlist.append(('Threads along x-axis',0))
    mlist.append(('Threads along y-axis',0))
    mlist.append(('Threads along z-axis',0))
    mlist.append(('User set thread group size', -1))

    for line in f:
        try:
            # General float cases
            for field in float_fields:
                if re.match('^'+re.escape(field)+':', line)!=None:
                    val = float(line.split(':')[1].split()[0])
                    mlist.append((field,val))

            # General Integer cases
            for field in int_fields:
                if re.match('^'+re.escape(field)+':', line)!=None:
                    val = int(line.split(':')[1].split()[0])
                    mlist.append((field,val))

            # General string cases
            for field in str_fields:
                if re.match('^'+re.escape(field)+':', line)!=None:
                    val = line.split(':')[1].strip()
                    mlist.append((field,val))

            # cache size filed
            if 'Assumed usable cache size' in line:
                cs = line.split(':')[1].strip()
                cs = cs.split('K')[0]
                mlist.append(('cache size',int(cs)))

            # special case for modified fields
            if 'Intra-diamond halo concatenation' in line:
                val = int(line.split(':')[1].split()[0])
                mlist.append(('Halo concatenation',val))
        except:
            print(field, line)
            raise

        # Global size
#Global domain    size:262144    nx:64    ny:64    nz:64
        if 'Global domain' in line:
            dims = line.split()[3:6]
            nx = dims[0].split(':')[1]
            ny = dims[1].split(':')[1]
            nz = dims[2].split(':')[1]
            mlist.append(('Global NX',int(nx)))
            mlist.append(('Global NY',int(ny)))
            mlist.append(('Global NZ',int(nz)))

        # Local size
#Rank 0 domain    size:131072    nx:32    ny:64    nz:64
        if 'Rank 0 domain' in line:
            dims = line.split()[4:7]
            nx = dims[0].split(':')[1]
            ny = dims[1].split(':')[1]
            nz = dims[2].split(':')[1]
            mlist.append(('Local NX',int(nx)))
            mlist.append(('Local NY',int(ny)))
            mlist.append(('Local NZ',int(nz)))

        # topology information
        if 'Processors topology' in line:
            dims = line.split(':')[1].strip()
            npx,npy,npz = dims.split(',')
            mlist.append(('npx',npx))
            mlist.append(('npy',npy))
            mlist.append(('npz',npz))

        # MWD statistics information
        measures =[
                   'Wavefront barrier wait [s]:', 'Wavefront barrier wait [%]:',
                   'Wavefront steady state [s]:', 'Wavefront steady state [%]:',
                   'Wavefront startup/end [s]:', 'Wavefront startup/end [%]:',
                   'Wavefront communication [s]:', 'Wavefront communication [%]:',
                   'Wavefront others [s]:', 'Wavefront others [%]:',
                   'Group spin-wait [s]:', 'Group spinn-wait [%]:',
                   'Resolved diamonds:'
                   ]
        for m in measures:
            if m in line:
                vals = [i.strip() for i in line.split(':')[1].split(' ') if len(i.strip()) !=0]
                for i in range(len(vals)):
                    mlist.append((m[:-1] + ' group %d'%i ,vals[i]))

        # PLUTO parameters
        if 'PLUTO tile size of loop' in line:
            mlist.append((line.split(':')[0].strip(), line.split(':')[1].strip()))


        # LIKWID performance results

        if 'Measuring group' in line:
            val = line.split(' ')[2].strip()
            mlist.append(('LIKWID performance counter',val))

        if '|    Memory BW [MBytes/s] STAT' in line or '|    Memory bandwidth [MBytes/s] STAT' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Sustained Memory BW',val))
        if 'Memory data volume [GBytes] STAT' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Total Memory Transfer',val))


        snames = [('%s bandwidth [MBytes/s] STAT ','BW'),
                  ('%s data volume [GBytes] STAT', 'data volume'),
                  ('%s Evict [MBytes/s] STAT', 'evict'),
                  ('%s Load [MBytes/s] STAT', 'load')]

        for cache in ['L2', 'L3']:
            for sn in snames:
                sn0 = sn[0]%(cache)
                if (sn0 in line):
                    vals = [i.strip() for i in line.split('|')[2:6]]
                    for i, st in enumerate(['sum', 'max', 'min', 'avg']):
                        mlist.append((cache+' '+sn[1]+' '+st, vals[i]))
                elif sn0 in line:
                    vals = [i.strip() for i in line.split('|')[2:-1]]
                    for i in range(len(vals)):
                        mlist.append((cache+' '+sn[1]+' c'+str(i), vals[i]))

        snames = [('CPI STAT', 'CPI'),
                  ('Load to Store ratio STAT', 'Load to Store ratio'),
                  ('L2_DATA_WRITE_MISS_MEM_FILL ', 'L2_DATA_WRITE_MISS_MEM_FILL'),
                  ('L2_DATA_READ_MISS_MEM_FILL ', 'L2_DATA_READ_MISS_MEM_FILL'),
                  ('HWP_L2MISS ', 'HWP_L2MISS'),
                  ('L2_VICTIM_REQ_WITH_DATA ', 'L2_VICTIM_REQ_WITH_DATA'),
                  ('SNP_HITM_L2 ', 'SNP_HITM_L2'),
                  ('CPU_CLK_UNHALTED ', 'CPU_CLK_UNHALTED')]
        for sn in snames:
            if ((sn[0].lower() in line.lower()) and ('STAT' in line)):
                vals = [i.strip() for i in line.split('|')[2:6]]
                for i, st in enumerate(['sum', 'max', 'min', 'avg']):
                    mlist.append((sn[1]+' '+st, vals[i]))
            elif sn[0] in line:
                vals = [i.strip() for i in line.split('|')[2:-1]]
                for i in range(len(vals)):
                    mlist.append((sn[1]+' c'+str(i), vals[i]))


        if 'L1 DTLB miss rate STAT' in line:
            vals = [i.strip() for i in line.split('|')[2:6]]
            for i, st in enumerate(['sum', 'max', 'min', 'avg']):
                mlist.append(('L1 DTLB miss rate %s'%(st),vals[i]))
                mlist.append(('L1 DTLB load miss rate %s'%(st),''))

        if 'L1 DTLB load miss rate STAT' in line:
            vals = [i.strip() for i in line.split('|')[2:6]]
            for i, st in enumerate(['sum', 'max', 'min', 'avg']):
                mlist.append(('L1 DTLB load miss rate %s'%(st),vals[i]))
                mlist.append(('L1 DTLB miss rate %s'%(st),''))

        if 'Energy [J]' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Energy',val))

        if 'Power [W]' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Power',val))

        if 'Energy DRAM [J]' in line:
            val = line.split('|')[2].strip()
            mlist.append(('Energy DRAM',val))

        if 'Power DRAM [W]' in line:
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


    # LIKWID 4 compatibility
    hw_ctr_fields = {
                    '':'',
                    'TLB': 'L1 DTLB load miss rate avg',
                    'DATA': 'Load to Store ratio avg',
                    'L2': 'L2 data volume sum',
                    'L3': 'L3 data volume sum',
                    'MEM': 'Total Memory Transfer',
                    'ENERGY': 'Energy'}
    if 'LIKWID performance counter' not in mlist.keys():
        mlist['LIKWID performance counter'] = ''
        for ctr in hw_ctr_fields:
            if(hw_ctr_fields[ctr]==[]): continue
            field = hw_ctr_fields[ctr]
            if field in mlist.keys():
                if mlist[field] !='':
                    mlist['LIKWID performance counter'] = ctr



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
