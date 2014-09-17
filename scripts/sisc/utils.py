def run_experiment(np, nx, ny, nz, nt, ts, outfile, target_dir):
    import os
    from string import Template
    
    ensure_dir(target_dir)
    outpath = os.path.join(os.path.curdir,target_dir,outfile)
    
    cmd_template=Template('mpirun -np $np ./build/mwd_kernel --npx $np --verbose 1 --verify 0 --n-tests 2 --nx $nx --ny $ny --nz $nz --nt $nt --target-ts $ts | tee $outfile')

    cmd_str = cmd_template.substitute(np=np, nx=nx, ny=ny, nz=nz, nt=nt, ts=ts, 
                                      outfile=outpath)    
    os.system(cmd_str)


def llsub_shaheen_experiment(npx, npy, npz, nx, ny, nz, nt, ts, t_dim, outfile, target_dir):
    import os
    from string import Template

    job_template=Template("""#!/usr/bin/env bash
#
# @ job_name            = Stencil_scaling
# @ job_type            = bluegene
# @ output              = $outfile.out
# @ error               = $outfile.err
# @ environment         = COPY_ALL; 
# @ wall_clock_limit    = 59:00,59:00
# @ notification        = always
# @ bg_size             = $np
# @ account_no          = k156

# @ queue

export LD_LIBRARY_PATH=$${KSL_PPC450D_LD_LIBRARY_PATH}
export PYTHONPATH=$${KSL_PPC450D_PYTHONPATH}
/bgsys/drivers/ppcfloor/bin/mpirun -exp_env LD_LIBRARY_PATH -exp_env PYTHONPATH -env BG_MAPPING=TXYZ -np $np -mode VN $exec_path --n-tests 2 --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz --nt $nt --target-ts $ts --t-dim $t_dim""")

    ensure_dir(target_dir)
    outpath = os.path.join(os.path.abspath("."),target_dir,outfile)
    
    np = npx*npy*npz
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")
    job_script = job_template.substitute(np=np, npx=npx, npy=npy, npz=npz, nx=nx, ny=ny, nz=nz, nt=nt, ts=ts, t_dim=t_dim, 
outfile=outpath, exec_path=exec_path)
    outscript = outpath+".ll"
    with open(outscript, 'w') as f:
        f.write(job_script)
    
    cmd_str = "llsubmit "+ outscript
    os.system(cmd_str)


def llsub_neser_experiment(npx, npy, npz, nx, ny, nz, nt, ts, t_dim, outfile, target_dir):
    import os
    from string import Template

    job_template=Template("""#!/bin/sh
#@ job_name         = stencil_scaling
#@ output           = $outfile.out
#@ error            = $outfile.err
#@ job_type         = parallel
#@ node             = $nn
#@ tasks_per_node   = $ppn
#@ wall_clock_limit = 19:00
#@ account_no = k156
#@ queue
source /etc/profile.d/modules.sh

module switch openmpi/1.5.3-intel
module load intel-compilers/11.1 

$$MPIEXEC -np $np $exec_path --n-tests 2 --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz --nt $nt --target-ts $ts --t-dim $t_dim""")

    ensure_dir(target_dir)
    outpath = os.path.join(os.path.abspath("."),target_dir,outfile)
    
    np = npx*npy*npz
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")
    ppn = 4
    nn = np/ppn
    job_script = job_template.substitute(nn=nn, ppn=ppn, np=np, npx=npx, npy=npy, npz=npz, nx=nx, ny=ny, nz=nz, nt=nt, ts=ts, t_dim=t_dim, outfile=outpath, exec_path=exec_path)    
    outscript = outpath+".ll"
    with open(outscript, 'w') as f:
        f.write(job_script)
    
    cmd_str = "llsubmit "+ outscript
    os.system(cmd_str)
#    print cmd_str

def llsub_noor_experiment(npx, npy, npz, nx, ny, nz, nt, ts, t_dim, outfile, target_dir):
    import os
    from string import Template

    job_template=Template("""#!/bin/bash -l
#BSUB -q rh6_q2hr
######BSUB -W 59:00
#BSUB -n $np
#BSUB -e $outfile.err
#BSUB -o $outfile.out
#BSUB -J stencil_scaling
#BSUB -app default
#BSUB -a openmpi

cd $cur_dir
module load gcc/4.6.0 mpi-openmpi/1.4.3-gcc-4.6.0

mpirun.lsf -np $np $exec_path --n-tests 2 --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz --nt $nt --target-ts $ts --t-dim $t_dim""")

    ensure_dir(target_dir)
    outpath = os.path.join(os.path.abspath("."),target_dir,outfile)
    
    cur_dir=os.path.abspath(".")

    np = npx*npy*npz
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")
    
    job_script = job_template.substitute(np=np, npx=npx, npy=npy, npz=npz, nx=nx, ny=ny, nz=nz, nt=nt, ts=ts, t_dim=t_dim, outfile=outpath, exec_path=exec_path, cur_dir=cur_dir)    
    outscript = outpath+".sh"
    with open(outscript, 'w') as f:
        f.write(job_script)
    
    cmd_str = "bsub < "+ outscript
    print cmd_str
    os.system(cmd_str)


def create_project_tarball(dest_dir, fname):
    import tarfile, glob
    import os

    nl = ["src/kernels/*.c", "src/*.c", "src/*.h", "scripts/*.py", "make.inc", "Makefile"]
    nl = [glob.glob(n) for n in nl]
    nl = [n for nn in nl for n in nn]

    out_dir = os.path.join(os.path.abspath("."),dest_dir)
    ensure_dir(out_dir)
    out_name = os.path.join(out_dir, fname+".tar.gz")

    print "Writing project files to:" + out_name
    with tarfile.open(out_name, "w:gz") as tar:
        for n in nl:
            print "Adding to the tar file: " + n
            tar.add(n)


def ensure_dir(d):
    import os, errno
    try:
        os.makedirs(d)
    except OSError as exc:
        if exc.errno == errno.EEXIST:
            pass
        else: raise
