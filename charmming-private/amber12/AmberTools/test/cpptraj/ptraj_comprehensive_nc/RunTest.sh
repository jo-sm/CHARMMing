#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles trajectory.netcdf trajectory_test.mdcrd
 
# Test 2
CheckNetcdf
INPUT="ptraj_netcdf.in"
TOP="compound.prmtop"
RunCpptraj "PTRAJ comprehensive mdcrd -> netcdf test."
INPUT="ptraj_mdcrd.in"
RunCpptraj "PTRAJ comprehensive netcdf -> mdcrd test."
DoTest ../../ptraj_comprehensive/trajectory.mdcrd trajectory_test.mdcrd 
CheckTest

EndTest

exit 0
