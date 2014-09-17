#!/usr/bin/env python

def verification_std(k):
    from verification_utils import run_verification

    # general test
    print "Verifying general domain size and network topology combinations"            
    dim_list = [[16, 32, 64],
                [32, 64, 16],
                [64, 16, 32]]
    # Selected MPI topologies
    top_list= [[1,1,1],
               #[2,1,1],
               [1,2,1],
               [1,1,2],
               [2,2,1],
               #[2,1,2],
               [1,2,2],
               [2,2,2],
               [4,1,1],
               [1,4,1],
               #[1,1,4],
               [4,2,1],
               #[8,1,1],
               #[1,8,1],
               #[1,1,8]
               ]
    for npx,npy,npz in top_list:
        for nlx,nly,nlz in dim_list:
            for kernel in [0,1,2,3,4,5]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel)

    # testing halo concatenation
    if(k == 1):
        short_dim_list2= [[16, 32, 64],
                          [64, 16, 32]]
        # Selected MPI topologies
        short_top_list2=[[1,2,1],
                        [2,2,1],
                        [2,2,2],
                        [4,1,1],
                        ]
        print "Verifying halo concatenation combinations"            
        for npx,npy,npz in short_top_list2:
            for nlx,nly,nlz in short_dim_list2:
                for kernel in [0,1,2,3,4,5]:
                    nx = nlx*npx
                    ny = nly*npy
                    nz = nlz*npz
                    run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, concat=1)

    # Testing double precision + halo concatenation
    if(k == 1):
        print "Verifying double precision + halo concatenation combinations"            
        short_dim_list1= [[16, 32, 64],
                          [64, 16, 32]]
        short_top_list1=[[1,2,1],
                        [2,2,1],
                        [2,2,2],
                        [4,1,1],
                        ]
        for npx,npy,npz in short_top_list1:
            for nlx,nly,nlz in short_dim_list1:
                for kernel in [1,2,3,4,5]:
                    nx = nlx*npx
                    ny = nly*npy
                    nz = nlz*npz
                    run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, dp=1, concat=1)


    # testing double precision
    print "Verifying double precision combinations"            
    short_dim_list2= [[16, 32, 64],
                      [64, 16, 32]]
    # Selected MPI topologies
    short_top_list2=[[1,2,1],
                    [2,2,1],
                    [2,2,2],
                    [4,1,1],
                    ]
    for npx,npy,npz in short_top_list2:
        for nlx,nly,nlz in short_dim_list2:
            for kernel in [1,2,3]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, kernel=kernel, dp=1)


#    # test large domains to verify the blocking code (assumig usable cache size 8MB)
#    print "Testing large XY slices to test the cache blocking code"
#        for npx,npy,npz in [[1,2,1]]:
#            run_verification(npx=npx, npy=npy, npz=npz, nx=1024, ny=512, nz=16, ts=k)
