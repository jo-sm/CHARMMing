#!/bin/bash

. ../MasterTest.sh

CleanFiles corr.in corr.dat 
INPUT="-i corr.in"
cat > corr.in <<EOF
parm ../tz2.parm7
trajin ../tz2.nc

distance d1 :2 :12 

corr d1 d1 out corr.dat
EOF
RunCpptraj "Correlation Test."
DoTest corr.dat.save corr.dat
CheckTest
EndTest

exit 0
