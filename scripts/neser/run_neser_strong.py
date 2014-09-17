#!/usr/bin/env python

base_name = 'neser_strongscaling'
base_dir = 'results/' + base_name
 
def main():
    from utils import create_project_tarball
    create_project_tarball(base_dir, "project_" + base_name)
 
    #exp(base_size=512, npz_list=[16,32,64,128]) #omitted
    #exp(base_size=512, npz_list=[   32,64,128])
    #exp(base_size=512, npz_list=[      64,128])

    exp(base_size=256, npz_list=[16,32,64,128])
    #exp(base_size=256, npz_list=[   32,64,128])
    #exp(base_size=256, npz_list=[      64,128])

    #exp(base_size=64, npz_list=[16,32       ])
    #exp(base_size=64, npz_list=[   32,64    ])
    #exp(base_size=64, npz_list=[      64,128])


def exp(base_size, npz_list):
    # Strong scaling study
    from utils import llsub_neser_experiment

    # constants
    # create code snapshot
    sub_dir = str(base_size) + '_' + '_'.join(map(str,npz_list))
    target_dir= base_dir + '/' + sub_dir
    npx=1
    npy=1
    nt=362 # nt-2 is a common divisor of most of the unroll sizes below
    # unrolling experiments for a set of local NZ's
    unroll_lists = {512:[2,30,60],
                    256:[2,12,30],
                    128:[2,8,12],
                    64:[2,4,6],
                    32:[2]}


    nx=base_size
    ny=base_size
    nz = npz_list[0] * base_size
    # loop over standard time steppers
    for k in [1,2]:
        for i in npz_list:
            npz=i
            out_file =base_name+'_ts%d_np%dx%dx%d_n%dx%dx%d' % (k,npx,npy,npz,nx,ny,nz)
            outfile=(out_file)
            llsub_neser_experiment(npx=npx, npy=npy, npz=npz, nx=nx, ny = ny, nz = nz,
                   nt = nt, ts=k, t_dim=0, outfile=out_file, target_dir=target_dir)
    #loop over pyramid methods
    for k in [4,5]:
        for i in [16,32,64,128]:
            local_nz = nz/i
            for p in unroll_lists[local_nz]:
                npz=i
                out_file= base_name + '_ts%d_np%dx%dx%d_n%dx%dx%d_%dtunroll' % (k,npx,npy,npz,nx,ny,nz,p)
                llsub_neser_experiment(npx=npx, npy=npy, npz=npz, nx=nx, ny = ny, nz = nz, 
                       nt = nt, ts=k, t_dim=p, outfile=out_file, target_dir=target_dir)


if __name__ == "__main__":
    main()
