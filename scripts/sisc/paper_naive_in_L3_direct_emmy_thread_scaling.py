#!/usr/bin/env python

def submit_emmy_experiment(nt, nx, ny, nz, is_dp, th, kernel, outfile, target_dir):
    import os
    import subprocess
    from string import Template
    from utils import ensure_dir    

    job_template=Template(
"""export OMP_NUM_THREADS=$thn; likwid-perfctr -s 0x03 -g MEM -C S0:0-$th $exec_path --n-tests 2 --disable-source-point --npx 1 --npy 1 --npz 1 --nx $nx --ny $ny --nz $nz  --verbose 1 --target-ts 0 --nt $nt --target-kernel $kernel | tee $outfile""")


    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    if(is_dp==1):
        exec_path = os.path.join(os.path.abspath("."),"build_dp/mwd_kernel")
    else:
        exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")


    job_cmd = job_template.substitute(nt=nt, nx=nx, ny=ny, nz=nz, th=(th-1), thn=th, kernel=kernel, outfile=outpath, exec_path=exec_path, target_dir=target_dir)

    print '# ' + job_cmd
    sts = subprocess.call(job_cmd, shell=True)


def main():
    from utils import create_project_tarball 

    import socket
    hostname = socket.gethostname()

    exp_name = "emmy_socket_thread_scaling_naive_inL3_1to10_threads_f2p2GHz_at_%s" % hostname  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 

    target_dir='results/' + exp_name 

    size_l = {0:(128, 64, 64, 1e4), 1:(96, 96, 96, 4e4), 4:(64, 32, 32, 5e4), 5:(64, 48,48, 1e5)}

    for th in list(range(10, 0, -1)):
        for is_dp in [1]:
            for kernel in [0, 1, 4, 5]:
                nx, ny, nz, nt = size_l[kernel]
                outfile=(exp_name + 'naive_kernel%d_isDP%d_threads%d' % (kernel, is_dp, th))
                submit_emmy_experiment(nt=nt, nx=nx, ny=ny, nz=nz, is_dp=is_dp, th=th, kernel=kernel, outfile=outfile, target_dir=target_dir)

if __name__ == "__main__":
    main()
