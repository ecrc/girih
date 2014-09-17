#!/usr/bin/env python

def main():
    from utils import run_experiment

    # Weak scaling study
    for k in [0,1,2]:
        for i in[1,2,4,8]:
            n = 128*i
            #print("i=%d  n=%d" % (i,n))
            run_experiment(np = i,
                           nx = n,
                           ny = 128,
                           nz = 128,
                           nt = 100,
                           ts=k,
                           outfile=('weakscaling_ts%d_%dp' % (k,i)),
                           target_dir='results/weak_scaling')

if __name__ == "__main__":
    main()
