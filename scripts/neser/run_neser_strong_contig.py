#!/usr/bin/env python

def llsub_neser_experiment(npx, npy, npz, nx, ny, nz, nt, ts, t_dim, outfile, target_dir):
    import os
    from string import Template
    from utils import ensure_dir    

    job_template=Template("""#!/bin/sh
#@ job_name         = stencil_scaling
#@ output           = $outfile.out
#@ error            = $outfile.err
#@ job_type         = parallel
#@ node             = $nn
#@ tasks_per_node   = $ppn
#@ wall_clock_limit = 29:00
#@ account_no = k156
#@ queue
source /etc/profile.d/modules.sh
module unload openmpi
module load openmpi/1.5.3-intel intel-compilers/11.1

$$MPIEXEC -np $np $exec_path --n-tests 2 --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz --nt $nt --target-ts $ts --t-dim $t_dim --z-mpi-contig 1""")

    ensure_dir(target_dir)
    outpath = os.path.join(os.path.abspath("."),target_dir,outfile)

    np = npx*npy*npz
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")
    ppn = 8
    nn = np/ppn
    job_script = job_template.substitute(nn=nn, ppn=ppn, np=np, npx=npx, npy=npy, npz=npz, nx=nx, ny=ny, nz=nz, nt=nt, ts=ts, t_dim=t_dim, outfile=outpath, exec_path=exec_path)
    outscript = outpath+".ll"
    with open(outscript, 'w') as f:
        f.write(job_script)

    cmd_str = "llsubmit "+ outscript
    os.system(cmd_str)

def main():

    from utils import create_project_tarball   
 
    # Strong scaling study
    target_dir='results/neser_strong_scaling_contig'
    create_project_tarball(target_dir, "project_neser_strongs_contig")
    min_p = 16
    nnz = min_p*256
    npx=1
    npy=1
    nx=256
    ny=256
    nt=802

    # loop over standard time steppers
    for k in [0,1,2,3]:
        for i in [16,32,64,128]:
            npz=i
            outfile=('neser_strongscaling_ts%d_%d_%d_%dp_contig' % (k,npx,npy,npz))
            llsub_neser_experiment(
                   npx = npx,
                   npy = npy,
                   npz = npz,
                   nx = nx,
                   ny = ny,
                   nz = nnz,
                   nt = nt,
                   ts=k,
                   t_dim=0,
                   outfile=outfile,
                   target_dir=target_dir)

if __name__ == "__main__":
    main()
