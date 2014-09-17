#!/usr/bin/env python

def llsub_neser_experiment(npx, npy, npz, ppn, nx, ny, nz, nt, ts, t_dim, is_contig, outfile, target_dir):
    import os
    from string import Template
    from utils import ensure_dir    

    job_template=Template("""
$$MPIEXEC -np $np $exec_path --n-tests 2 --npx $npx --npy $npy --npz $npz --nx $nx --ny $ny --nz $nz --nt $nt --target-ts $ts --t-dim $t_dim --z-mpi-contig $is_contig | tee $outfile.out""")

    origin_dir = os.path.abspath(".")
    target_dir = os.path.join(origin_dir, target_dir)
    ensure_dir(target_dir)
    outpath = os.path.join(target_dir,outfile)

    np = npx*npy*npz
    exec_path = os.path.join(os.path.abspath("."),"build/mwd_kernel")
    nn = np/ppn
    cmd_str = job_template.substitute(np=np, npx=npx, npy=npy, npz=npz, nx=nx, ny=ny, nz=nz, nt=nt, ts=ts, t_dim=t_dim, is_contig=is_contig, outfile=outpath, exec_path=exec_path, target_dir=target_dir, origin_dir=origin_dir)

    #os.system(cmd_str)
    #os.system(cmd_str)
    #os.system(cmd_str)
    print '---------------------'
    print 'cd ' + target_dir
    print cmd_str
    print 'cd ' + origin_dir
    print '---------------------'

def main():

    from utils import create_project_tarball   
    tarball_dir='results/neser_ppn_map_tau'
    create_project_tarball(tarball_dir, "project_neser_ppn_map_tau")
 
    # Strong scaling study
    min_p = 8
    nnz = min_p*256
    npx=1
    npy=1
    nx=256
    ny=256
    nt=402
    npz= 8
    k = 2 # halo-first TS

    is_contig = 0
    for ppn in [1, 2, 4, 8]:
        target_dir='results/neser_ppn_map_tau/ts%d_%d_%d_%dp_%dppn' % (k,npx,npy,npz,ppn)
        outfile=('neser_ppn_map_tau_ts%d_%d_%d_%dp_%dppn_tau' % (k,npx,npy,npz,ppn))
        llsub_neser_experiment(
                   npx = npx,
                   npy = npy,
                   npz = npz,
                   ppn = ppn,
                   nx = nx,
                   ny = ny,
                   nz = nnz,
                   nt = nt,
                   ts=k,
                   t_dim=0,
                   is_contig=is_contig,
                   outfile=outfile,
                   target_dir=target_dir)

if __name__ == "__main__":
    main()
