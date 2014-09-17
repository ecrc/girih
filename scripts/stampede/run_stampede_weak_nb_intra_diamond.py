#!/usr/bin/env python

def submit_stampede_experiment(npx, npy, npz,  nx, ny, nz, nt, ts, t_dim, is_contig, outfile, target_dir, jobname):
    import os
    from string import Template
    from utils import ensure_dir    

    job_template=Template("""#!/bin/bash

#SBATCH -J $jobname           # Job name
#SBATCH -o $outfile.out       # Name of stdout output file
#SBATCH -e $outfile.err       # Name of stderr output file
#SBATCH -p development        # Submit to the 'normal' or 'development' queue
#SBATCH -N $nn                # Total number of nodes
#SBATCH -n $np                # Total number of mpi tasks requested
#SBATCH -t 00:59:00           # Run time (hh:mm:ss)
#SBATCH -A TG-ASC130040

# Set the number of threads per task(Default=1)
export OMP_NUM_THREADS=7

export MPICH_ASYNC_PROGRESS=1 

cd $target_dir

# Run the hybrid application 
ibrun tacc_affinity $exec_path --n-tests 2 --disable-source-point --halo-concatenate 1 --target-ts $ts --t-dim $t_dim --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz --nt $nt --z-mpi-contig $is_contig

""")

    target_dir = os.path.join(os.path.abspath("."),target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    np = npx*npy*npz
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")
    ppn=2
    nn = np/ppn
    if np==1: nn=1
    job_script = job_template.substitute(nn=nn, np=np, npx=npx, npy=npy, npz=npz, nx=nx, ny=ny, nz=nz, nt=nt, ts=ts, t_dim=t_dim, is_contig=is_contig, outfile=outpath, exec_path=exec_path, target_dir=target_dir, jobname=jobname)
    outscript = outpath+".ll"
    with open(outscript, 'w') as f:
        f.write(job_script)

    cmd_str = "sbatch "+ outscript
    os.system(cmd_str)
    #print cmd_str

def main():

    from utils import create_project_tarball 
    exp_name = "stampede_weakscaling_nb_intra_diamond"  
    tarball_dir='results/'+exp_name
    create_project_tarball(tarball_dir, "project_"+exp_name)
 
    npx=1
    npz=1
    nx=512
    ny_b=512
    nz=512
    nt=2402
    k = 7 # NB-intra-diamond

    target_dir='results/' + exp_name 
    is_contig = 0

    t_dim = 15 #unroll size in time

    for npy in [2, 4, 8, 16, 32]:
        ny = ny_b*npy
        outfile=(exp_name + '_ts%d_npy%d_tdim%d' % (k,npy,t_dim))
        submit_stampede_experiment(
                   npx = npx,
                   npy = npy,
                   npz = npz,
                   nx = nx,
                   ny = ny,
                   nz = nz,
                   nt = nt,
                   ts=k,
                   t_dim=t_dim,
                   is_contig=is_contig,
                   outfile=outfile,
                   target_dir=target_dir,
                   jobname=exp_name)

if __name__ == "__main__":
    main()
