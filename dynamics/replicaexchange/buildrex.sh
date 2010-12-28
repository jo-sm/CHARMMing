#!/bin/bash

export PATH=/v/apps/compilers/x86_64/gcc-4.3.2/bin:$PATH
export PATH=/v/apps/mpi/openmpi-1.3.1/64bit-gfortran43/bin:$PATH

mpif90 NPTrex.f -o NPTrex.x

#ln -s /v/bigbox4/home/rsingh/replicaexchange/carbtest/AdHocNPTrex/gfortamd-xxl.x11 charmm
