#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles rms.in rmsd1.dat rmsd2.dat ref.nc rmsd.mass.dat

CheckNetcdf
TOP="../tz2.parm7"
CRD="../tz2.nc"
INPUT="rms.in"

# Test 1 - 2drms
cat > rms.in <<EOF
noprogress
trajin $CRD 1 10 
2drms :3-7 rmsout rmsd2.dat
EOF
RunCpptraj "2D RMSD Test."
DoTest rmsd.dat.save rmsd2.dat
CheckTest

# Test 2 - 2drms, mass-weighted
cat > rms.in <<EOF
noprogress
trajin $CRD 1 10 
2drms :3-7 rmsout rmsd.mass.dat mass
EOF
RunCpptraj "2D RMSD Test, mass-weighted."
DoTest rmsd.mass.dat.save rmsd.mass.dat
CheckTest

# Test 3 - 2drms to reference traj
cat > rms.in <<EOF
trajin $CRD 1 10
trajout ref.nc netcdf
EOF
RunCpptraj "Generating ref traj for 2D RMSD test."
cat > rms.in <<EOF
noprogress
trajin $CRD 1 10 
2drms :3-7 rmsout rmsd1.dat reftraj ref.nc
EOF
RunCpptraj "2D RMSD Test with reference trajectory."
DoTest rmsd.dat.save rmsd1.dat
CheckTest

EndTest

exit 0
