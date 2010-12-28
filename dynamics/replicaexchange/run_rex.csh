#!/bin/csh

# set path for mpirun
#set path = ( /usr/bin/)

# set full path to executable
set exe = ./NPTrex.x

# run the MPI based C program on $n procs, w. config file $c
#mpirun -machinefile $PBS_NODEFILE -np $1 $exe $2
mpirun -np $1 $exe $2

