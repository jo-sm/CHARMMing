#!/bin/bash

. ../MasterTest.sh

CleanFiles 1rrb_vac_distmat.dat 1rrb_vac_mwcovarmat.dat 1rrb_vac_mwcovarmat_evecs.dat 1rrb_vac_distcovarmat.dat  1rrb_vac_distcovarmat_evecs.dat

CheckPtrajAnalyze

INPUT="../../ptraj_matrix/ptraj.in"
TOP="../../ptraj_matrix/1rrb_vac.prmtop"

RunCpptraj "matrix analysis test."
DoTest ../../ptraj_matrix/1rrb_vac_distmat.dat.save 1rrb_vac_distmat.dat

CheckTest
EndTest

