#!/usr/bin/env python

def main():
    from verification_utils import run_verification

    dim_list = [[16, 32, 64]]    
    # Selected MPI topologies
    top_list= [[1,1,1],
               [1,1,2],
               [2,2,1],
               [2,2,2],
               ]
    for k in [0,1]:
        for npx,npy,npz in top_list:
            for nlx,nly,nlz in dim_list:
                for kernel in [0,2,3,4,5]:
                    nx = nlx*npx
                    ny = nly*npy
                    nz = nlz*npz
                    run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel)

    dim_list = [[16, 32, 64]]    
    # Selected MPI topologies
    top_list= [[1,1,1],
               [1,1,2],
               ]
    for k in [0, 1]:
        for npx,npy,npz in top_list:
            for nlx,nly,nlz in dim_list:
                for kernel in [1,2,3]:
                    nx = nlx*npx
                    ny = nly*npy
                    nz = nlz*npz
                    run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, dp=1, concat=1)

    for k in [0, 1]:
        for npx,npy,npz in top_list:
            for nlx,nly,nlz in dim_list:
                for kernel in [1]:
                    nx = nlx*npx
                    ny = nly*npy
                    nz = nlz*npz
                    run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, concat=1)

 
    for k in [0, 1]:
        for npx,npy,npz in top_list:
            for nlx,nly,nlz in dim_list:
                for kernel in [2]:
                    nx = nlx*npx
                    ny = nly*npy
                    nz = nlz*npz
                    run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, dp=1)

 
if __name__ == "__main__":
    main()
