#!/usr/bin/env python

def submit_emmy_experiment(nx, ny, nz, nt, is_dp, ts, tb, kernel, outfile, target_dir, tgs=1, nwf=4, mwdt=2):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=10; likwid-perfctr  -m -s 0x03 -g MEM -C S0:0-9 -- $exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --wavefront 1 --nt $nt --target-kernel $kernel --t-dim $tb --num-wavefronts $nwf --thread-group-size $tgs --mwd-type $mwdt | tee $outfile""")

    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")

    job_cmd = job_template.substitute(mwdt=mwdt, nwf=nwf, tgs=tgs, nx=nx, ny=ny, nz=nz, nt=nt, ts=ts, tb=tb, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    print job_cmd
    sts = subprocess.call(job_cmd, shell=True)


def main():
    from utils import create_project_tarball 

    import socket
#    hostname = socket.gethostname()

    exp_name = "emmy_cache_size_vs_code_balace_2p2GHz"#_at_%s" % hostname  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 
    target_dir='results/' + exp_name 


    kernel = 1
    exp = []
    exp = exp + [[0, 0, 960, 30,  0]]
#    exp = exp + [[0, 0, 480, 400, 0]]
#    exp = exp + [[2, i, 480, 400, 10] for i in range(1, 30, 2)] 
    exp = exp + [[2, i, 960, 50,  10] for i in range(1, 30, 2)] 
    for ts, tb, N, nt, tgs in exp:
        outfile=('ts%d_kernel%d_tb%d_tgs%d_N%d.txt' % (ts, kernel, tb, tgs, N))
        submit_emmy_experiment(N, N, N, nt, is_dp=1, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir, tgs=tgs, nwf=tgs)


    kernel = 5
    exp = []
#    exp = exp + [[0, 0, 384, 1000, 0]]
    exp = exp + [[0, 0, 680, 30,   0]]
#    exp = exp + [[2, i, 384, 1000, 10] for i in range(1, 20, 2)] 
    exp = exp + [[2, i, 680, 50,   10] for i in range(1, 20, 2)] 
    for ts, tb, N, nt, tgs in exp:
        outfile=('ts%d_kernel%d_tb%d_tgs%d_N%d.txt' % (ts, kernel, tb, tgs, N))
        submit_emmy_experiment(N, N, N, nt, is_dp=1, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir, tgs=tgs, nwf=tgs)


    kernel = 4
    exp = []
#    exp = exp + [[0, 0, 240, 500,  0]]
    exp = exp + [[0, 0, 480, 30,   0]]
#    exp = exp + [[2, i, 240, 500, 10] for i in range(1, 8, 2)] 
    exp = exp + [[2, i, 480, 50,  10] for i in range(1, 6, 2)] 
    for ts, tb, N, nt, tgs in exp:
        outfile=('ts%d_kernel%d_tb%d_tgs%d_N%d.txt' % (ts, kernel, tb, tgs, N))
        submit_emmy_experiment(N, N, N, nt, is_dp=1, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir, tgs=tgs, nwf=4*tgs, mwdt=0)



if __name__ == "__main__":
    main()
