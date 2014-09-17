#!/usr/bin/env python

def main():
    import sys
    from verification_utils import run_verification

    k = int(sys.argv[1])

    wf = 1

    print "Verifying kernel number: %d" % k

    dim_list = [[16, 32,  128,1],
                [16, 64,  32, 3]]
    # Selected MPI topologies
    top_list= [[1,1,1],
               [1,2,1]]
    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [0,1,2,3,4,5]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, t_dim=t, wavefronts=wf)


    dim_list = [[16, 32,  128,1],
                [16, 64,  32, 3]]
    # Selected MPI topologies
    top_list= [[1,2,1],
               [1,4,1]]
    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [1,2]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, t_dim=t, dp=1, concat=1, wavefronts=wf)

    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [1,2]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, t_dim=t, dp=1, wavefronts=wf)


    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [0,2]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, t_dim=t, concat=1, wavefronts=wf)

  
if __name__ == "__main__":
    main()
