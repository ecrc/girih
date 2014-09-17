#!/usr/bin/env python

def submit_emmy_experiment(ts, kernel, tb, nx, ny, nz, nt, is_dp, outfile, target_dir):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=10; likwid-pin -s 0x03 -c S0:0-9 $exec_path --n-tests 3 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --nt $nt --target-kernel $kernel --t-dim $tb | tee $outfile""")


    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


    job_cmd = job_template.substitute(ts=ts, tb=tb, nx=nx, ny=ny, nz=nz, nt=nt, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    sts = subprocess.call(job_cmd, shell=True)
    #print job_cmd

def main():

    from utils import create_project_tarball 

    import socket
    hostname = socket.gethostname()

    exp_name = "emmy_testing_separating_stride1_impact_no_seperation_f2p2GHz_at_%s" % (hostname)  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 

    target_dir='results/' + exp_name 


#    ts = 0
#    tb = 0
#    for kernel in [0, 1, 2]:
#        for N in [8, 16, 32, 64, 256]:
#            outfile=(exp_name + '_ts%d_kernel%d_TB%d_%d' % (ts, kernel, tb, N))
#            nt = max(800, int(8e9 / N**3))
#            submit_emmy_experiment(ts=ts, kernel=kernel, tb=tb, nx=N, ny=N, nz=N, nt=nt, is_dp=0, outfile=outfile, target_dir=target_dir)

    ts = 6
    tb = 9
    for kernel in [1, 2]:
        for N in [20, 40, 80, 100, 260]:
            outfile=(exp_name + '_ts%d_kernel%d_TB%d_%d' % (ts, kernel, tb, N))
            nt = max(800, int(8e9 / N**3))
            submit_emmy_experiment(ts=ts, kernel=kernel, tb=tb, nx=N, ny=N, nz=N, nt=nt, is_dp=0, outfile=outfile, target_dir=target_dir)


if __name__ == "__main__":
    main()
