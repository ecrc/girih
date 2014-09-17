#!/usr/bin/env python

def submit_emmy_experiment(wf, is_dp, th, ts, tb, kernel, outfile, target_dir):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=$thn; likwid-perfctr -s 0x03 -g MEM_SP -C S0:0-$th $exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx 512 --ny 512 --nz 512  --verbose 1 --target-ts $ts --wavefront $wf --nt 200 --target-kernel $kernel --t-dim $tb | tee $outfile""")


    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


    job_cmd = job_template.substitute(wf=wf, th=(th-1), thn=th, ts=ts, tb=tb, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    print job_cmd
    sts = subprocess.call(job_cmd, shell=True)


def main():
    from utils import create_project_tarball 

    import socket
    hostname = socket.gethostname()

    #exp_name = "emmy_socket_thread_scaling_520cube_1c_ac_wf_naive_1to10_threads_at_%s" % hostname  
    #exp_name = "emmy_socket_thread_scaling_520cube_ac_wf_1to10_threads_at_%s" % hostname  
    exp_name = "emmy_socket_thread_scaling_1cwd_large_1to10_threads_f2p2GHz_at_%s" % hostname  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 

    target_dir='results/' + exp_name 

    wf = 1 
    tb = 3
    for th in list(range(10, 0, -1)):
        for is_dp in [0]:
            for ts in [9]:
                for kernel in [1]:
                    outfile=(exp_name + '_ts%d_kernel%d_isDP%d_threads%d' % (ts, kernel, is_dp, th))
                    submit_emmy_experiment(wf=wf, is_dp=is_dp, th=th, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)

if __name__ == "__main__":
    main()
