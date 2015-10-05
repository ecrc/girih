#!/usr/bin/env python

def submit_emmy_experiment(nt, N, tgs, is_dp, procs, ts, kernel, outfile, target_dir):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=10; export I_MPI_PIN_DOMAIN=socket; export KMP_AFFINITY=granularity=fine,scatter; mpirun_rrze -np $procs -npernode 2 -- $exec_path --n-tests 2 --disable-source-point --npx 1 --npy $procs --npz 1 --nx $N --ny $N --nz $N  --verbose 1 --halo-concatenate 1 --target-ts $ts --wavefront 1 --nt $nt --target-kernel $kernel --thread-group-size $tgs | tee $outfile""")

    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")

    job_cmd = job_template.substitute(nt=nt, N=N, tgs=tgs, procs=procs, ts=ts, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    print procs, kernel, ts, tgs
    print job_cmd
    sts = subprocess.call(job_cmd, shell=True)


def experiment(target_dir, exp_name, ts, tgs, kernel, nt, N, proc_l, cs):
    is_dp = 1
    for procs in proc_l:
        outfile=(exp_name + 'isDP%d_ts%d_kernel%d_cs%d_tgs%d_procs%d.txt' % (is_dp, ts, kernel, cs, tgs,  procs))
        submit_emmy_experiment(nt=nt*procs, N=N, tgs=tgs, is_dp=is_dp, procs=procs, ts=ts, kernel=kernel, outfile=outfile, target_dir=target_dir) 


def main():
    from utils import create_project_tarball 
    import socket
    import time,datetime

    time_stamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H_%M')

    exp_name = "distributed_strongscaling_at_IVB10_%s" % (time_stamp)  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
    target_dir='results/' + exp_name 


    exp =       [[0, 5, 768, 0, 8192, proc] for proc in [1,2,4,8,12,16,24,32]]
    exp = exp + [[2, 5, 768, 1, 2048, proc] for proc in [1,2,4,8,12,16]]
    exp = exp + [[2, 5, 768,-1, 2048, proc] for proc in [1,2,4,8,12,16,24,32]]

#    exp = exp + [[2, 5, 768, 2, 2048, proc] for proc in [1,2,4,8,12,16,24,32]]
#    exp = exp + [[2, 5, 768, 5, 2048, proc] for proc in [1,2,4,8,12,16,24,32]]
#    exp = exp + [[2, 5, 768, 10,2048, proc] for proc in [1,2,4,8,12,16,24,32]]

    exp = exp + [[0, 4, 512, 0, 8192, proc] for proc in [1,2,4,8,16]]
    exp = exp + [[2, 4, 512, 1, 0, proc] for proc in [1,2]]
    exp = exp + [[2, 4, 512,-1, 0, proc] for proc in [1,2,4,8,16]]

#    exp = exp + [[2, 4, 512, 2, 0, proc] for proc in [1,2,4]]
#    exp = exp + [[2, 4, 512, 5, 0, proc] for proc in [1,2,4,8,16]]
#    exp = exp + [[2, 4, 512, 10, 0, proc] for proc in [1,2,4,8,16,32]]

    exp_f = ['ts', 'kernel', 'N', 'tgs', 'cs', 'procs']
    exp = [dict(zip(exp_f,e)) for e in exp]
  

#    nodes_range = [0, 2]
#    nodes_range = [2, 8]
    nodes_range = [8, 16]

    nt = 50

    for e in exp:
        if (e['procs'] > nodes_range[0]*2) and (e['procs'] <= nodes_range[1]*2):
            experiment(target_dir, exp_name, ts=e['ts'], tgs=e['tgs'], kernel=e['kernel'], nt=nt, N=e['N'], proc_l=[e['procs']], cs=e['cs'])


if __name__ == "__main__":
    main()
