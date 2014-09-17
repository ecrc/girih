#!/usr/bin/env python

def submit_emmy_experiment(nwf, tgs, wf, is_dp, th, ts, tb, kernel, outfile, target_dir):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=$thn; likwid-perfctr -s 0x03 -g MEM -C S0:0-$th $exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx 960 --ny 960 --nz 960  --verbose 1 --target-ts $ts --wavefront $wf --nt 100 --target-kernel $kernel --t-dim $tb --thread-group-size $tgs --num-wavefronts $nwf | tee $outfile""")
#"""export OMP_NUM_THREADS=$thn; numactl -N 0 $exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx 480 --ny 480 --nz 480  --verbose 1 --target-ts $ts --wavefront $wf --nt 100 --target-kernel $kernel --t-dim $tb --thread-group-size $tgs | tee $outfile""")


    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


    job_cmd = job_template.substitute(nwf=nwf, tgs=tgs, wf=wf, th=(th-1), thn=th, ts=ts, tb=tb, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    print job_cmd
    sts = subprocess.call(job_cmd, shell=True)

def main():

    from utils import create_project_tarball 

    import socket
    hostname = socket.gethostname()

    exp_name = "emmy_socket_thread_scaling_960cube_f2.2GHz_mwd_swd_mwd_groups_1to10_threads_7pt_const_at_%s" % hostname  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 

    target_dir='results/' + exp_name 

    kl = [1]

    wf = 1
    for th in list(range(10, 0, -2)):
        for is_dp in [1]:
            for ts, tgs, tb, nwf  in [(9, 2, 7, 8)]:  # DP
                for kernel in kl:
                    outfile=(exp_name + '_ts%d_tb%d_tgs%d_nwf%d_kernel%d_isDP%d_threads%d' % (ts, tb, tgs, nwf, kernel, is_dp, th))
                    submit_emmy_experiment(nwf=nwf, tgs=tgs, wf=wf, is_dp=is_dp, th=th, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)

    for th in[5, 10]:
        for is_dp in [1]:
            for ts, tgs, tb, nwf  in [(9, 5, 11, 20)]:  #DP
                for kernel in kl:
                    outfile=(exp_name + '_ts%d_tb%d_tgs%d_nwf%d_kernel%d_isDP%d_threads%d' % (ts, tb, tgs, nwf, kernel, is_dp, th))
                    submit_emmy_experiment(nwf=nwf, tgs=tgs, wf=wf, is_dp=is_dp, th=th, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)

    for th in[10]:
        for is_dp in [1]:
            for ts, tgs, tb, nwf  in [(9, 10, 15, 40)]:  #DP
                for kernel in kl:
                    outfile=(exp_name + '_ts%d_tb%d_tgs%d_nwf%d_kernel%d_isDP%d_threads%d' % (ts, tb, tgs, nwf, kernel, is_dp, th))
                    submit_emmy_experiment(nwf=nwf, tgs=tgs, wf=wf, is_dp=is_dp, th=th, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)

    for th in list(range(10, 0, -1)):
        for is_dp in [1]:
            for ts, tgs, tb, nwf  in [(0,0,0,0), (9, 1, 5, 4)]:
                for kernel in kl:
                    outfile=(exp_name + '_ts%d_tb%d_tgs%d_nwf%d_kernel%d_isDP%d_threads%d' % (ts, tb, tgs, nwf, kernel, is_dp, th))
                    submit_emmy_experiment(nwf=nwf, tgs=tgs, wf=wf, is_dp=is_dp, th=th, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)



if __name__ == "__main__":
    main()
