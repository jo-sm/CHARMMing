#
#                            PUBLIC DOMAIN NOTICE
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the authors' official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software is freely available
#  to the public for use.  There is no restriction on its use or
#  reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, NIH and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. NIH, NHLBI, and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any
#  particular purpose.

# maximum number of atoms
maxatoms = 50000

# set to 1 if you have Q-Chem enabled in your CHARMM exe and want
# to allow your users to use QM/MM
haveqchem = 1

# root of the CHARMMing install
charmming_root = "/var/www/charmming"

# place where user directories are...
user_home = "/home/schedd"

# place where executables etc. are kept
data_home = "/usr/local/charmming"

# topology file to be used for proteins and nucleic acids ... these
# must live in data_home/toppar.
default_pro_top = "top_all36_prot.rtf"
default_pro_prm = "par_all36_prot.prm"
default_na_top = "top_all36_na.rtf"
default_na_prm = "par_all36_na.prm"

# Where AmberTools lives -- used for Antechamber topology/parameter
# determination.
amber_home = "/usr/local/charmming/amber12"

# methods to use to try to generate topology/parameters for BADHET
# atoms. These methods will be tried in order until one succeeds 
# (or they all fail).
toppar_generators = 'cgenff,match,antechamber,genrtf'
dd_toppar_generators = 'match,cgenff'#,antechamber,genrtf'

# CGenFF host & port
cgenff_host = 'dogmans.umaryland.edu'
cgenff_port = 32108

# Path to all libraries needed by the CHARMM executables
lib_path = ''

# path to the single threaded CHARMM executable
#charmm_exe = "/usr/local/charmming/c37b2-prelease.exe"
#charmm_apbs_exe = "/usr/local/charmming/c37b2-prelease.exe"

charmm_exe = "/usr/local/charmming/c37b2-qc-apbs.one"
charmm_apbs_exe = "/usr/local/charmming/c37b2-qc-apbs.one"
charmm_mscale_exe = "/usr/local/charmming/c37b2-qc-mscale.ompi"

# path to the MPI-enabled CHARMM executable
charmm_mpi_exe = "/usr/local/charmming/gfortran-xxlg-qc.ompi"

# path to mpirun (not used by default -- testing parallel)
mpirun_exe = "/bin/false"

# path to stride (needed for CG stuff)
stride_bin = "/usr/local/bin/stride"

#Default number of processors for MSCALE jobs. Replace once we have proper multi-proc intefaces.
default_mscale_nprocs = 2

# dd targets
user_dd_targets_home = user_home + "/dd/targets"

# dd ligands
user_dd_ligands_home = user_home + "/dd/ligands"

# dd jobs
user_dd_jobs_home = user_home + "/dd/jobs"

# path to the DAIM executable and properties files
daim_exe = "/usr/local/charmming/drug_design/DAIM/bin/daim"
daim_param = user_home + "/dd/app_files/daim/daim.param"
daim_prop = user_home + "/dd/app_files/daim/daim.prop"
daim_weight = user_home + "/dd/app_files/daim/daim.prop"

# path to the SEED executable and parameter files
seed_exe = "/usr/local/charmming/drug_design/seed_send_3.3.5/seed_3.3.5"
seed_param = user_home + "/dd/app_files/seed/seed.par"

# path to the FFLD executable
ffld_exe = "/usr/local/charmming/drug_design/ffld_send_3.3/ffld_3.3_gcc4.0_static"
ffld_param = user_home + "/dd/app_files/ffld/FFLD_param_cgenff"

#path to VMD executable
vmd_exe = "/usr/local/charmming/drug_design/bin/vmd"

# path to the FLEA exeutable and parameter files
flea_param = user_home + "/dd/app_files/flea/PARM.flea"
flea_exe = "/usr/local/charmming/drug_design/FLEA/flea_1.0_32b_static"

# path to the CHARMM files
charmm_files = user_home + "/dd/app_files/charmm/"
charmm_param = user_home + "/dd/app_files/charmm/top_all36_prot.rtf"

# dd scripts
dd_scripts_home = user_home + "/dd/scripts"

# dd shell launch command
dd_submit_parallel="qsub -cwd -l h_rt=00:20:00" #to enable parallel execution of docking commands
dd_submit_serial = "" #for serial execution of docking commands

## Limits on various calculations

# max atoms for all atom normal modes
max_nma_atoms = 1500

# step limits for minimization and dynamics
minimize_steplimit = 1000
dyna_steplimit = 1000

# docking parameters
docking_iterations = 1
ffld_energy_evaluations=2000
ffld_generations = 100
clustering_energy_cutoff = 10

# qsar
# user_qsar_home = user_home + "/qsar"

# qsar jobs
# user_qsar_jobs_home = user_qsar_home + "/jobs"

# qsar jobs
# user_qsar_models_home = user_qsar_home + "/models"
