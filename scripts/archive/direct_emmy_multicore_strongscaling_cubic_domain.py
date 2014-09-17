#!/usr/bin/env python

def submit_emmy_experiment(ts, tb, perf_ctr, nx, ny, nz, nt, is_dp, outfile, target_dir):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=10; likwid-perfctr -c S0:0-9 -s 0x03 -g $perf_ctr $exec_path --n-tests 3 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --nt $nt --target-kernel $kernel --t-dim $tb | tee $outfile""")


    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


    job_cmd = job_template.substitute(tb=tb, perf_ctr=perf_ctr, nx=nx, ny=ny, nz=nz, nt=nt, kernel=1, ts=ts, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    #sts = subprocess.call(job_cmd, shell=True)
    print job_cmd

def main():

    from utils import create_project_tarball 

    import socket
    hostname = socket.gethostname()

    exp_name = "emmy_strongscaling_f2p2GHz_at_%s" % (hostname)  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 

    target_dir='results/' + exp_name 

    perf_ctr = 'MEM'


    # single precision work
    is_dp = 0
    tb = 19
    for ts in [0,6]:
        for N in range(40, 521, 20):
            outfile=(exp_name + 'isDP%d_ts%d_N%d_TB%d' % (is_dp, ts,N, tb))
            nt = max(30, 4e9 / N**3)
            submit_emmy_experiment(tb=tb, ts=ts, perf_ctr=perf_ctr, nx=N, ny=N, nz=N, nt=nt, is_dp=is_dp, outfile=outfile, target_dir=target_dir)

    for N, tb in [[40, 1], [80, 3], [120, 5], [160, 7]]:
        for ts in [9]:
            outfile=(exp_name + 'isDP%d_ts%d_N%d_TB%d' % (is_dp, ts,N, tb))
            nt = max(30, 4e9 / N**3)
            submit_emmy_experiment(tb=tb, ts=ts, perf_ctr=perf_ctr, nx=N, ny=N, nz=N, nt=nt, is_dp=is_dp, outfile=outfile, target_dir=target_dir)

    tb = 9
    for N in range(200, 521, 20):
        for ts in [9]:
            outfile=(exp_name + 'isDP%d_ts%d_N%d_TB%d' % (is_dp, ts,N, tb))
            nt = max(30, 4e9 / N**3)
            submit_emmy_experiment(tb=tb, ts=ts, perf_ctr=perf_ctr, nx=N, ny=N, nz=N, nt=nt, is_dp=is_dp, outfile=outfile, target_dir=target_dir)




if __name__ == "__main__":
    main()
