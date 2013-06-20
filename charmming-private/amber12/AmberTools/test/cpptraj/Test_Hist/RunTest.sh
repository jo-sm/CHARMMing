#!/bin/bash

. ../MasterTest.sh

CleanFiles hist.in hist.gnu hist.agr freeE.gnu norm.gnu
CheckNetcdf
cat > hist.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc

distance d1 :1 :10 
distance d2 :10 :13  

dihedral phi6 :5@C :6@N :6@CA :6@C 
dihedral psi6 :6@N :6@CA :6@C :7@N 

hist d1 d2:8:10:0.1:* min 9.0 max 26.0 step 0.5 out hist.gnu
hist d1:9:26:0.5 out hist.agr
hist phi6 psi6 min -180 max 180 bins 72 out freeE.gnu free 300.0
hist phi6 psi6 min -180 max 180 bins 72 out norm.gnu norm 
EOF

INPUT="-i hist.in"
RunCpptraj "Histogram Analysis Test"
DoTest hist.gnu.save hist.gnu
DoTest hist.agr.save hist.agr
DoTest freeE.gnu.save freeE.gnu
DoTest norm.gnu.save norm.gnu

EndTest

exit 0