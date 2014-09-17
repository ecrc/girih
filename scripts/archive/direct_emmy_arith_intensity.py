#!/usr/bin/env python

def submit_emmy_experiment(is_dp, th, ts, tb, kernel, outfile, target_dir):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=$th; likwid-perfctr -g MEM -C S0:0-0 $exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx 512 --ny 512 --nz 512  --verbose 1 --target-ts $ts --wavefront 1 --nt 100 --target-kernel $kernel --t-dim $tb --thread-group-size 1 | tee $outfile""")

    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


    job_cmd = job_template.substitute(th=th, ts=ts, tb=tb, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    sts = subprocess.call(job_cmd, shell=True)
    #print job_cmd

def main():

    from utils import create_project_tarball 

    import socket
    hostname = socket.gethostname()

    exp_name = "emmy_socket_arithmetic_intensity_512cube_at_%s" % hostname  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 

    target_dir='results/' + exp_name 
    tb = 1

    exp = [[0, 0, 0], # ts, kernel, tb
           [0, 1, 0],
           [0, 2, 0],

           [9, 0, 1],
           [9, 0, 3],
#           [9, 0, 7],
#           [9, 0,15],

           [9, 1, 1],
           [9, 1, 3],
           [9, 1, 7],
           [9, 1,15],

           [9, 2, 1],
           [9, 2, 3],
           [9, 2, 7],
           [9, 2,15]]

    th = 1
    for is_dp in [0,1]:
        for ts, kernel, tb in exp:
            outfile=(exp_name + '_ts%d_kernel%d_isDP%d_tb%d_threads%d' % (ts, kernel, is_dp, tb, th))
            submit_emmy_experiment(is_dp=is_dp, th=th, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)

if __name__ == "__main__":
    main()
