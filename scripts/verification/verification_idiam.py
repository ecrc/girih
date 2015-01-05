#!/usr/bin/env python

def verification_idiam(k):
    from verification_utils import run_verification

    npx =1
    npy =1
    npz =1


    dim_list = [[32, 480,  128,5, 8, 2],
                [64, 480,  128,5, 32, 2]]
    print "Verifying different wavefront paralellization stratigies"
    th = 10
    for mwdt in [0,1,2,3]:
        for nlx, nly, nlz, t, bs_x, tgs in dim_list:
            for kernel, R in [(1,1), (0,4)]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                if th <= nly/((t+1)*2*R): # enough concurrency in the wavefront for the provided threads
                    print "MWD type:%d" % mwdt
                    run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t, kernel=kernel,  num_threads=th, tgs=tgs, nwf=tgs*R, bs_x=bs_x, mwd_type=mwdt)


    dim_list = [[16, 128,  128,1, 8, 1],
                [32, 480,  128,5, 8, 2],
                [64, 480,  128,5, 32, 2]]
    print "Verifying different values of cache blocking in X"
    for th in range(2,21,2):
        for nlx, nly, nlz, t, bs_x, tgs in dim_list:
            for kernel, R in [(1,1), (4,4)]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                if th <= nly/((t+1)*2*R): # enough concurrency in the wavefront for the provided threads
                    run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t, kernel=kernel,  num_threads=th, tgs=tgs, nwf=tgs, bs_x=bs_x)


    print "Verifying different multi-core wavefront threads/tile_size combinations"
    dim_list = [[16, 128,  128,1],
                    [16, 480,  128,5],
                    [16, 640,  128,7]]
    tgs=1
    for th in range(1,20):
        for nlx, nly, nlz, t in dim_list:
            for kernel in [1]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                if th <= nly/((t+1)*2): # enough concurrency in the wavefront for the provided threads
                    run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t, kernel=kernel,  num_threads=th, tgs=tgs, nwf=tgs)

    for th, tgs in [(16, 16), (16, 8), (16, 4), (16, 2), (8,8), (8,4), (8,2)]:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [1]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
              #  if th < (t+1)**2: # enough concurrency in the wavefront for the provided threads
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t, kernel=kernel,  num_threads=th, tgs=tgs, nwf=tgs)



    print "Verifying general domain size and network topology combinations"
    dim_list = [[16, 32,  128,1],
                [16, 64,  128, 3],
                [16, 64,  128, 7],
                [32, 128, 128, 7]]
    # Selected MPI topologies
    top_list= [[1,1,1],
               [1,2,1],
               [1,4,1],
               [1,8,1]]   
    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [0,1,2,3,4,5]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                np = npx*npy*npz
                tgs = max(1, min(8,16/np))
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t, kernel=kernel,  tgs=tgs, nwf=tgs)


    print "Verifying temporal blocks with increasing size"
    dim_list = [[16, 96, 128, 5],
                [16, 128, 128, 7],
                [16, 160, 256, 9],
                [16, 224, 512, 13],
                [16, 480, 512, 29]]
    # Selected MPI topologies
    top_list= [[1,1,1],
               [1,4,1]]   
    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [0,1,2,3,4,5]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                np = npx*npy*npz
                tgs = max(1, min(10,20/np))
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t, kernel=kernel,  tgs=tgs, nwf=tgs)


    print "Verifying double precision + halo concatenation combinations"
    dim_list = [[16, 32,  128,3],
                [32, 128, 128, 7]]
    # Selected MPI topologies
    top_list= [[1,2,1],
               [1,4,1]]
    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [1,2,3,4,5]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                np = npx*npy*npz
                tgs = max(1, min(10/20/np))
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t, kernel=kernel, dp=1, concat=1,  tgs=tgs, nwf=tgs)


    print "Verifying double precision combinations"
    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [1,2,3,4,5]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                np = npx*npy*npz
                tgs = max(1, min(10/20/np))
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t, kernel=kernel, dp=1,  tgs=tgs, nwf=tgs)



    print "Verifying halo concatenation combinations"
    for npx,npy,npz in top_list:
        for nlx, nly, nlz, t in dim_list:
            for kernel in [0,1,2,3,4,5]:
                nx = nlx*npx
                ny = nly*npy
                nz = nlz*npz
                np = npx*npy*npz
                tgs = max(1, min(10,20/np))
                run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t, kernel=kernel, concat=1,  tgs=tgs, nwf=tgs)


#    # verify the blocking in Y at large XY plains
#    for npx,npy,npz in [[1,2,1]]:
#        for nlx, nly, nlz, t in [[1024,512,32,63]]:
#            nx = nlx*npx
#            ny = nly*npy
#            nz = nlz*npz
#            run_verification(nx=nx, ny=ny, nz=nz, ts=k, npx=npx, npy=npy, npz=npz, t_dim=t)
