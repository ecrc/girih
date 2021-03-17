pipeline {
    agent { label 'jenkinsfile' }
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

nsockets=`lscpu | awk '/^Socket/{print $2}'`
ncorepersocket=`lscpu | awk '/^Core/{print $4}'`
ncores=$(( nsockets * ncorepersocket ))
let nthread_per_group=$ncores/4           

module purge
module load intel/2019
module load intelmpi/2018-update-1
make clean 
make
'''
    }
           }
               
        stage ('test') {
            steps {
            sh '''#!/bin/bash -el
module purge
module load intel/2018
module load intel-mpi/2019.8.254/gcc-7.5.0-mmiar35  
export OMP_NUM_THREADS=${ncores}
export OMP_PROC_BIND=true
export OMP_PLACES=cores
mpirun -np 1 ./build/mwd_kernel --npy 1 --nx 512 --ny 512 --nz 512 --nt 500 --target-kernel 0 --mwd-type 0 --target-ts 2 --verify 0  --threads ${ncores} --thread-group-size ${nthread_per_group}
'''
}
    }
        }
          }
