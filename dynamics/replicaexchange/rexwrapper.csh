#!/bin/csh
#first argument is filepath to chdir to, second is name of config file, third is number of processors
chdir $1
./rex_init.csh $2
#first argument is number of nodes, second is the config file name
./run_rex.csh $3 $2
