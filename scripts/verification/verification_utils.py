#!/usr/bin/env python

def run_verification(nx, ny, nz, ts, t_dim=3, npx=1, npy=1, npz=1, kernel=0, dp=0, concat=0, num_threads=0, tgs=1, nwf=-1, bs_x=1000000):
    import os
    import subprocess
    from string import Template

    # Compute the required number of time steps for the wave to propogate from to the whole domain
    stencil_radius=[4,1,1,1,4,1] # IMPORTANT: this list has to match the stencil kernels list of the code
    nt = max(10, nx/stencil_radius[kernel]/2)

    binary="./build"    
    if(dp == 1):
        binary = binary + "_dp"
    
    binary=binary + "/mwd_kernel"

    if(nwf == -1): nwf = 4*tgs

    cmd_template=Template('mpirun_rrze -block -np $np $binary --nt $nt --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz --verbose 0 --verify 1 --target-ts $ts --target-kernel $kernel --t-dim $t_dim --halo-concatenate $concat --thread-group-size $tgs --num-wavefronts $nwf --bsx $bs_x')

    np = npx*npy*npz
    cmd_str = cmd_template.substitute(np=np, nt=nt, npx=npx, npy=npy, npz=npz, nx=nx, ny=ny, nz=nz, t_dim=t_dim, ts=ts, kernel=kernel, binary=binary, concat=concat, tgs=tgs, nwf=nwf, bs_x=bs_x)
    if num_threads == 0:
        num_threads = max(1, min(8,16/np))
    omp_threads = "export OMP_NUM_THREADS=%d" % num_threads
    #print(omp_threads + ";" + cmd_str)
    sts = subprocess.call(omp_threads + ";" + cmd_str, shell=True)
    #print('')
    #os.system(omp_threads + ";" + cmd_str)

