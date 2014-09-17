#!/usr/bin/env python

def submit_emmy_experiment(N, nt, nwf, tgs, wf, is_dp, th, ts, tb, kernel, outfile, target_dir):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=$thn; likwid-perfctr -m -s 0x03 -g ENERGY -C S0:0-$th $exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $N --ny $N --nz $N  --verbose 1 --target-ts $ts --wavefront $wf --nt $nt --target-kernel $kernel --t-dim $tb --thread-group-size $tgs --num-wavefronts $nwf | tee $outfile""")


    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


    job_cmd = job_template.substitute(N=N, nt=nt, nwf=nwf, tgs=tgs, wf=wf, th=(th-1), thn=th, ts=ts, tb=tb, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    print job_cmd
    sts = subprocess.call(job_cmd, shell=True)

def main():

    from utils import create_project_tarball 

    import socket
    hostname = socket.gethostname()

    exp_name = "emmy_energy_f2.2GHz_at_%s" % hostname  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 

    target_dir='results/' + exp_name 
    wf = 1
    is_dp = 1
    th = 10

    # (stencil, N, nt, experiments)
    exp = [(1, 960, 100, [(0,0,0,0)])];  th = 6 # special case for spatial blocking
    #exp = [(1, 960, 100, [(2, 2, 5, 12), (2, 5, 11, 20), (2, 10, 11, 50), (0,0,0,0), (2, 1, 5, 1)])]
    #exp = exp + [(5, 680, 50, [(2, 2, 3, 2), (2, 5, 3, 15), (2, 10, 9, 10), (0,0,0,0), (2, 1, 3, 1)])]
    #exp = exp + [(4, 480, 100, [(2, 2, 1, 2), (2, 5, 1, 5), (2, 10, 3, 10), (0,0,0,0), (2, 1, 1, 9)])]
    for kernel, N, nt, exp_l in exp:
        for ts, tgs, tb, nwf in exp_l:
	    outfile=(exp_name + '_ts%d_tgs%d_tb%d_nwf%d_kernel%d_isDP%d_threads%d.txt' % (ts, tgs, tb, nwf, kernel, is_dp, th))
	    submit_emmy_experiment(N=N, nt=nt, nwf=nwf, tgs=tgs, wf=wf, is_dp=is_dp, th=th, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)


if __name__ == "__main__":
    main()
