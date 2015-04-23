#!/usr/bin/env python

def main():
    import sys
    from verification_utils import run_verification

    dim_list = [[32, 128,  128,1],
                [32, 256,  32, 3]]

    print "Verifying different wavefront paralellization stratigies" 
    # Different MWD approaches

    for kernel in [0,1,5]:
        for nx, ny, nz, t in dim_list:
            for mwd_t in range(3):
                run_verification(nx=nx, ny=ny, nz=nz, ts=2, kernel=kernel, t_dim=t, mwd_type=mwd_t)

    print "Verifying MPI configurations"
    # Selected MPI topologies
    top_list= [[1,1,1],
               [1,2,1]]
    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [0,1,2,3,4,5]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                run_verification(nx=nx, ny=ny, nz=nz, ts=2, npx=npx, npy=npy, npz=npz, kernel=kernel, t_dim=t)


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
                run_verification(nx=nx, ny=ny, nz=nz, ts=2, npx=npx, npy=npy, npz=npz, kernel=kernel, t_dim=t, dp=1, concat=1)

    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [1,2]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                run_verification(nx=nx, ny=ny, nz=nz, ts=2, npx=npx, npy=npy, npz=npz, kernel=kernel, t_dim=t, dp=1)


    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [0,2]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                run_verification(nx=nx, ny=ny, nz=nz, ts=2, npx=npx, npy=npy, npz=npz, kernel=kernel, t_dim=t, concat=1)

  
if __name__ == "__main__":
    main()
