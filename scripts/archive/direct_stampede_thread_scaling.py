#!/usr/bin/env python

def submit_stampede_experiment(is_dp, th, ts, tb, kernel, outfile, target_dir):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=$th;  ibrun -o 0 -n 1  numactl -N 0 $exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx 512 --ny 512 --nz 512  --verbose 1 --target-ts $ts --wavefront 1 --nt 500 --target-kernel $kernel --t-dim $tb --thread-group-size 1 | tee $outfile""")

    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


    job_cmd = job_template.substitute(th=th, ts=ts, tb=tb, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    sts = subprocess.call(job_cmd, shell=True)
   # print job_cmd

def main():

    from utils import create_project_tarball 
    exp_name = "stampede_socket_thread_scaling_512cube_diamo_TB1"  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 

    target_dir='results/' + exp_name 
    tb = 1

    for th in list(range(8, 0, -1)):
        for is_dp in [0,1]:
            for ts in [9]:
                for kernel in [0]:
                    outfile=(exp_name + '_ts%d_kernel%d_isDP%d_threads%d' % (ts, kernel, is_dp, th))
                    submit_stampede_experiment(is_dp=is_dp, th=th, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)

if __name__ == "__main__":
    main()
