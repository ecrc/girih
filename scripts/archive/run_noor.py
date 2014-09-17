#!/usr/bin/env python

def main():
    from utils import llsub_noor_experiment 
    
    # Strong scaling study
    target_dir='results/noor_strong_scaling'
#    create_project_tarball(target_dir, "project_noor_strongs")
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
            outfile=('noor_strongscaling_ts%d_%d_%d_%dp' % (k,npx,npy,npz))
            llsub_noor_experiment(
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

    #loop over pyramid method
    k = 4
    for p in [2,4]:
        for i in [16,32,64,128]:
            npz=i
            outfile=('noor_strongscaling_ts%d_%d_%d_%dp_%dtunroll' % (k,npx,npy,npz,p))
            llsub_noor_experiment(
                   npx = npx,
                   npy = npy,
                   npz = npz,
                   nx = nx,
                   ny = ny,
                   nz = nnz,
                   nt = nt,
                   ts=k,
                   t_dim=p,
                   outfile=outfile,
                   target_dir=target_dir)


if __name__ == "__main__":
    main()