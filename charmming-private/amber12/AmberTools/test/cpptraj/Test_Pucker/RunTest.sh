#!/bin/bash

. ../MasterTest.sh

# Clean
CleanFiles pucker.in pucker.dat pucker2.in pucker2.dat

# Test 1
TOP=../adh026.pdb
INPUT=pucker.in
cat > pucker.in <<EOF
noprogress
trajin ../adh026.pdb
pucker p1 :1@C1' :1@C2' :1@C3' :1@C4' :1@O4' out pucker.dat 
pucker p2 :2@C1' :2@C2' :2@C3' :2@C4' :2@O4' out pucker.dat 
pucker p3 :3@C1' :3@C2' :3@C3' :3@C4' :3@O4' out pucker.dat 
pucker p4 :4@C1' :4@C2' :4@C3' :4@C4' :4@O4' out pucker.dat 
pucker p5 :5@C1' :5@C2' :5@C3' :5@C4' :5@O4' out pucker.dat 
pucker p6 :6@C1' :6@C2' :6@C3' :6@C4' :6@O4' out pucker.dat 
pucker p7 :7@C1' :7@C2' :7@C3' :7@C4' :7@O4' out pucker.dat 
pucker p8 :8@C1' :8@C2' :8@C3' :8@C4' :8@O4' out pucker.dat 
EOF
RunCpptraj "Pucker command test (Altona & Sundaralingam)"
DoTest pucker.dat.save pucker.dat

# Test 2
TOP=3OIE.parm7
INPUT=pucker2.in
cat > pucker2.in <<EOF
noprogress
trajin 3OIE.rst7
pucker p1 :1@C1' :1@C2' :1@C3' :1@C4' :1@O4' out pucker2.dat cremer
pucker p2 :2@C1' :2@C2' :2@C3' :2@C4' :2@O4' out pucker2.dat cremer
pucker p3 :3@C1' :3@C2' :3@C3' :3@C4' :3@O4' out pucker2.dat cremer
pucker p4 :4@C1' :4@C2' :4@C3' :4@C4' :4@O4' out pucker2.dat cremer
pucker p5 :5@C1' :5@C2' :5@C3' :5@C4' :5@O4' out pucker2.dat cremer
pucker p6 :6@C1' :6@C2' :6@C3' :6@C4' :6@O4' out pucker2.dat cremer
pucker p7 :7@C1' :7@C2' :7@C3' :7@C4' :7@O4' out pucker2.dat cremer
pucker p8 :8@C1' :8@C2' :8@C3' :8@C4' :8@O4' out pucker2.dat cremer
EOF
RunCpptraj "Pucker command test (Cremer & Pople)"
DoTest pucker2.dat.save pucker2.dat

CheckTest

EndTest

exit 0
