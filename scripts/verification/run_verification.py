#!/usr/bin/env python

def main():
    import sys
    from verification_utils import run_verification
    from verification_std import verification_std
    from verification_idiam import verification_idiam

    kernels = map(int, sys.argv[1:])

    for k in kernels:
        if(k < 2):
            verification_std(k)
        elif(k ==2):
            verification_idiam(k)


if __name__ == "__main__":
    main()
