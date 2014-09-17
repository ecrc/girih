#!/usr/bin/env python

def main():
    from utils import run_experiment, create_project_tarball

    # Strong scaling study
    target_dir='results/strong_scaling'
    create_project_tarball(target_dir, "project_ws_strongs")
    for k in [0,1,2]:
        for i in[1,2,4,8]:
            #print("i=%d  n=%d" % (i,n))
            run_experiment(np = i,
                           nx = 256,
                           ny = 128,
                           nz = 128,
                           nt = 100,
                           ts=k,
                           outfile=('strongscaling_ts%d_%dp' % (k,i)),
                           target_dir=target_dir)

if __name__ == "__main__":
    main()
