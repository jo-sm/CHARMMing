#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles prec.in prec.dat

CheckNetcdf
TOP="../tz2.truncoct.parm7"

# Test 1
cat > prec.in <<EOF
noprogress
trajin ../tz2.truncoct.nc
rms first :2-11 out prec.dat
datafile precision prec.dat * 8 3
EOF
INPUT="prec.in"
RunCpptraj "Data file output precision test."
DoTest prec.dat.save prec.dat

CheckTest

EndTest

exit 0
