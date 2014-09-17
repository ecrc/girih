#!/usr/bin/env python

def submit_emmy_experiment(nx, ny, nz, nt, is_dp, ts, tb, kernel, outfile, target_dir, tgs=1):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=10; likwid-perfctr  -m -s 0x03 -g MEM -C S0:0-9 -- $exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts $ts --wavefront 1 --nt $nt --target-kernel $kernel --t-dim $tb --num-wavefronts $tgs --thread-group-size $tgs | tee $outfile""")

    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


    job_cmd = job_template.substitute(tgs=tgs, nx=nx, ny=ny, nz=nz, nt=nt, ts=ts, tb=tb, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    print job_cmd
    sts = subprocess.call(job_cmd, shell=True)

def main():

    from utils import create_project_tarball 

    import socket
    hostname = socket.gethostname()

    exp_name = "emmy_socket_arithmetic_intensity_2p2GHz_at_%s" % hostname  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 

    target_dir='results/' + exp_name 
    is_dp = 1
    
    exp = [[0, 0, 480, 400]]
    exp = exp + [[2, i, 480, 400] for i in range(1,10,2)] 
    exp = exp + [[0, 0, 960, 50]]
    exp = exp + [[2, i, 960, 50] for i in range(1,8,2)] 
    kernel = 1
    for ts, tb, N, nt in exp:
        outfile=(exp_name + '_ts%d_kernel%d_tb%d_N%d' % (ts, kernel,tb,N))
#        submit_emmy_experiment(N, N, N, nt, is_dp=is_dp, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)



#    exp = [[0, 0, 240, 1000]]
#    exp = exp + [[2, i, 240, 1000] for i in range(1,8,2)] 
#    exp = exp + [[0, 0, 760, 50]]
    exp = [[0, 0, 760, 50]]
    exp = exp + [[2, i, 760, 50] for i in range(1,6,2)] 
    kernel = 5
    for ts, tb, N, nt in exp:
        outfile=(exp_name + '_ts%d_kernel%d_tb%d_N%d' % (ts, kernel,tb,N))
        submit_emmy_experiment(N, N, N, nt, is_dp=is_dp, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir)





    kernel = 1
    exp = [[2, i, 480, 400, 2] for i in range(1,10,2)] 
    exp = exp + [[2, i, 960, 50, 2] for i in range(1,8,2)] 

    exp = exp + [[2, i, 480, 400, 5] for i in range(1,12,2)] 
    exp = exp + [[2, i, 960, 50, 5] for i in range(1,8,2)] 

    exp = exp + [[2, i, 480, 400, 10] for i in range(1,18,2)] 
    exp = exp + [[2, i, 960, 50, 10] for i in range(1,14,2)] 
    for ts, tb, N, nt, tgs in exp:
        outfile=(exp_name + '_ts%d_kernel%d_tb%d_tgs%d_N%d' % (ts, kernel, tb, tgs, N))
        submit_emmy_experiment(N, N, N, nt, is_dp=is_dp, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir, tgs=tgs)



    kernel = 5
    exp = [[2, i, 240, 1000, 2] for i in range(1,10,2)] 
    exp = exp + [[2, i, 760, 50, 2] for i in range(1,6,2)] 
 
    exp = exp + [[2, i, 240, 1000, 5] for i in range(1,12,2)] 
    exp = exp + [[2, i, 760, 50, 5] for i in range(1,6,2)] 
  
    exp = exp + [[2, i, 240, 1000, 10] for i in range(1,14,2)] 
    exp = exp + [[2, i, 760, 50, 10] for i in range(1, 12, 2)] 
    for ts, tb, N, nt, tgs in exp:
        outfile=(exp_name + '_ts%d_kernel%d_tb%d_tgs%d_N%d' % (ts, kernel, tb, tgs, N))
        submit_emmy_experiment(N, N, N, nt, is_dp=is_dp, ts=ts, tb=tb, kernel=kernel, outfile=outfile, target_dir=target_dir, tgs=tgs)


if __name__ == "__main__":
    main()
