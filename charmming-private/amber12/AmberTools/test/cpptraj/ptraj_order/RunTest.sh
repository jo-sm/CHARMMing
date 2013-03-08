#!/bin/bash

# Test of ptraj vector / matrix
. ../MasterTest.sh

CleanFiles ptraj.* orderparam ired.vec noe v0.cjt v0.cmt

CheckPtrajAnalyze

TOP="../../ptraj_order/1IEE_A_prot.prmtop"
#INPUT="ptraj_ired_step1.in"
#ptraj $TOP $INPUT > ptraj.out
#INPUT="ptraj_ired_step2.in"
#ptraj $TOP $INPUT >> ptraj.out 
#mv orderparam ptraj.orderparam
#mv ired.vec ptraj.ired.vec
#mv noe ptraj.noe
INPUT="../../ptraj_order/ptraj_ired_step1.in"
RunCpptraj "ptraj vector/matrix test, step 1"
INPUT="../../ptraj_order/ptraj_ired_step2.in"
RunCpptraj "ptraj vector/matrix test, step 2"
DoTest ../../ptraj_order/orderparam.save orderparam
DoTest ../../ptraj_order/noe.save noe
#DoTest ptraj.orderparam orderparam
#DoTest ptraj.ired.vec ired.vec
#DoTest ptraj.noe noe
CheckTest
EndTest

exit 0
