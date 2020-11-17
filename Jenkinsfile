pipeline {
    agent { label 'knl' }
    triggers {
        pollSCM('H/10 * * * *')
    }
    options {
        disableConcurrentBuilds()
        buildDiscarder(logRotator(numToKeepStr: '50'))
        timestamps()
    }
    stages {
        stage ('build') {
            steps {
            sh '''#!/bin/bash -el
                    # The -x flags indicates to echo all commands, thus knowing exactly what is being executed.
                    # The -e flags indicates to halt on error, so no more processing of this script will be done
                    # if any command exits with value other than 0 (zero)
module purge
module load intel/2018
module load intelmpi/2018-update-1
make clean 
make
'''
    }
           }
               
        stage ('test') {
            steps {
            sh '''#!/bin/bash -el
export OMP_NUM_THREADS=64
export OMP_PROC_BIND=true
export OMP_PLACES=cores
mpirun -np 1 ./build/mwd_kernel --npy 1 --nx 512 --ny 512 --nz 512 --nt 500 --target-kernel 0 --mwd-type 0 --target-ts 2 --verify 0  --threads 64 --thread-group-size 16
'''
}
    }
        }
          }
