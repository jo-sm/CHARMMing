#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles nastruct.in BP.nastruct.dat BPstep.nastruct.dat bases.pdb baseaxes.pdb basepairaxes.pdb Helix.nastruct.dat

# Test 1
cat > nastruct.in <<EOF
noprogress
#parm 3OIE.pdb
#trajin 3OIE.pdb
parm ../adh026.pdb
trajin ../adh026.pdb 
nastruct naout nastruct.dat
#nastruct resrange 25-30
#nastruct resrange 1-12
#nastruct resrange 1-500
EOF
INPUT="-i nastruct.in"
RunCpptraj "NAstruct command test."
#DoTest test.out.save test.out
#DoTest bases.pdb.save bases.pdb
#DoTest baseaxes.pdb.save baseaxes.pdb
#DoTest basepairaxes.pdb.save basepairaxes.pdb
DoTest BP.nastruct.dat.save BP.nastruct.dat
DoTest BPstep.nastruct.dat.save BPstep.nastruct.dat
DoTest Helix.nastruct.dat.save Helix.nastruct.dat

CheckTest

EndTest

exit 0
