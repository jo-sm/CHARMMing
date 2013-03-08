"""

 This is a module that contains functions generally useful for the
 MMPBSA.py script. A full list of functions/subroutines is shown below.
 It must be included to insure proper functioning of MMPBSA.py

          Last updated: 09/04/2011                                     


                           GPL LICENSE INFO                             

Copyright (C) 2010  Dwight McGee, Billy Miller III, and Jason Swails

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
   
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330,
Boston, MA 02111-1307, USA.
"""

################################################################################
#  List of functions and a brief description of their purpose
#
#  remove: Removes temporary work files in this directory. It has a number of 
#          different levels to remove only a small number of files.
#  InterruptHandler: Exit cleanly displaying an "interrupt" message
#  CheckIncomp: Checks incompatibilities in the input files
#  PrintFileInfo: Prints information about files to the final output file
#  Abort: Invokes MPI.Abort for a failed calculation and prints out a useful 
#         message
#  concatenate: combines 2 files into a single, common file
################################################################################

import os, sys
from MMPBSA_mods.exceptions import InputError

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def remove(flag, mpi_size=0):
   """ Removes temporary files. Allows for different levels of cleanliness """

   # A list of all input files that we keep for the flag -use-mdins
   input_files = ['_MMPBSA_gb.mdin', '_MMPBSA_pb.mdin', 
                  '_MMPBSA_gb_decomp_com.mdin', '_MMPBSA_gb_decomp_rec.mdin',
                  '_MMPBSA_gb_decomp_lig.mdin', '_MMPBSA_pb_decomp_com.mdin',
                  '_MMPBSA_pb_decomp_rec.mdin', '_MMPBSA_pb_decomp_lig.mdin',
                  '_MMPBSA_ptrajentropy.in', '_MMPBSA_mutant_ptrajentropy.in']

   # All the extra files we keep for keep_files = 1
   keep_files_1=['_MMPBSA_ligand.mdcrd', '_MMPBSA_ligand.nc', 
                 '_MMPBSA_mutant_ligand.mdcrd', '_MMPBSA_mutant_ligand.nc',
                 '_MMPBSA_complex.mdcrd', '_MMPBSA_mutant_complex.mdcrd',
                 '_MMPBSA_complex.nc', '_MMPBSA_mutant_complex.nc',
                 '_MMPBSA_receptor.mdcrd', '_MMPBSA_mutant_receptor.mdcrd',
                 '_MMPBSA_receptor.nc', '_MMPBSA_mutant_receptor.nc',
                 '_MMPBSA_dummycomplex.inpcrd', '_MMPBSA_complex.pdb',
                 '_MMPBSA_dummyreceptor.inpcrd','_MMPBSA_receptor.pdb',
                 '_MMPBSA_dummyligand.inpcrd','_MMPBSA_ligand.pdb',
                 '_MMPBSA_mutant_dummycomplex.inpcrd',
                 '_MMPBSA_mutant_dummyreceptor.inpcrd',
                 '_MMPBSA_mutant_dummyligand.inpcrd',
                 '_MMPBSA_mutant_complex.pdb', '_MMPBSA_mutant_receptor.pdb',
                 '_MMPBSA_mutant_ligand.pdb', '_MMPBSA_complex_nm.mdcrd',
                 '_MMPBSA_complex_nm.nc', '_MMPBSA_mutant_complex_nm.mdcrd',
                 '_MMPBSA_mutant_complex_nm.nc', '_MMPBSA_receptor_nm.mdcrd',
                 '_MMPBSA_receptor_nm.nc', '_MMPBSA_mutant_receptor_nm.nc',
                 '_MMPBSA_mutant_receptor_nm.mdcrd', '_MMPBSA_ligand_nm.nc',
                 '_MMPBSA_ligand_nm.mdcrd', '_MMPBSA_mutant_ligand_nm.nc',
                 '_MMPBSA_mutant_ligand_nm.mdcrd', '_MMPBSA_avgcomplex.pdb',
                 '_MMPBSA_mutant_avgcomplex.pdb', '_MMPBSA_ligand_entropy.out',
                 '_MMPBSA_complex_entropy.out', '_MMPBSA_receptor_entropy.out',
                 '_MMPBSA_ptraj_entropy.out','_MMPBSA_mutant_ptraj_entropy.out',
                 '_MMPBSA_mutant_complex_entropy.out',
                 '_MMPBSA_mutant_receptor_entropy.out',
                 '_MMPBSA_mutant_ligand_entropy.out',
                 '_MMPBSA_complex_gb.mdout', '_MMPBSA_mutant_complex_gb.mdout',
                 '_MMPBSA_receptor_gb.mdout','_MMPBSA_mutant_receptor_gb.mdout',
                 '_MMPBSA_ligand_gb.mdout', '_MMPBSA_mutant_ligand_gb.mdout',
                 '_MMPBSA_complex_pb.mdout', '_MMPBSA_mutant_complex_pb.mdout',
                 '_MMPBSA_receptor_pb.mdout','_MMPBSA_mutant_receptor_pb.mdout',
                 '_MMPBSA_ligand_pb.mdout', '_MMPBSA_mutant_ligand_pb.mdout',
                 '_MMPBSA_complex_rism.mdout', 
                 '_MMPBSA_mutant_complex_rism.mdout',
                 '_MMPBSA_receptor_rism.mdout',
                 '_MMPBSA_mutant_receptor_rism.mdout',
                 '_MMPBSA_ligand_rism.mdout','_MMPBSA_mutant_ligand_rism.mdout',
                 '_MMPBSA_complex_nm.out', '_MMPBSA_mutant_complex_nm.out',
                 '_MMPBSA_receptor_nm.out','_MMPBSA_mutant_receptor_nm.out',
                 '_MMPBSA_ligand_nm.out', '_MMPBSA_mutant_ligand_nm.out',
                 '_MMPBSA_info']

   # Collect all of the temporary files (those starting with _MMPBSA_)
   allfiles = os.listdir(os.getcwd())
   tempfiles = []
   for fil in allfiles:
      if fil.startswith('_MMPBSA_'): tempfiles.append(fil)

   if flag == -1: # internal -- keep all mdin files
      for fil in tempfiles:
         if not fil in input_files: os.remove(fil)
   elif flag == 0: # remove all temporary files
      for fil in tempfiles: os.remove(fil)
   elif flag == 1: # keep keep mdcrds, mdouts, and other relevant output files
      for fil in tempfiles:
         if fil in keep_files_1: continue # keep this file
         # Now we have to split out this file and analyze the base. If the
         # suffix is just a number (corresponding to a thread-specific output
         # file or trajectory), then we only want to remove it if in the base
         # name is not in keep_files_1
         base, ext = os.path.splitext(fil)
         if ext.strip('.').isdigit() and base in keep_files_1: continue
         # if we've made it this far, remove the file
         os.remove(fil)
      
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def CheckIncomp(INPUT, stability):
   """ Checks for inconsistencies in the input parameters """
   if not INPUT['igb'] in [1, 2, 5, 7, 8]:
      raise InputError('Invalid value for IGB (%s)! ' % INPUT['igb'] +
                       'It must be 1, 2, 5, 7, or 8.')
   if INPUT['saltcon'] < 0: raise InputError('SALTCON must be non-negative!')
   if INPUT['surften'] < 0: raise InputError('SURFTEN must be non-negative!')
   if INPUT['indi'] < 0: raise InputError('INDI must be non-negative!')
   if INPUT['exdi'] < 0: raise InputError('EXDI must be non-negative!')
   if INPUT['scale'] < 0: raise InputError('SCALE must be non-negative!')
   if INPUT['linit'] < 0: raise InputError('LINIT must be a positive integer!')
   if not INPUT['prbrad'] in [1.4, 1.6]:
      raise InputError('PRBRAD (%s) must be 1.4 and 1.6!' % INPUT['prbrad'])
   if INPUT['istrng'] < 0: raise InputError('ISTRNG must be non-negative!')
   if not INPUT['inp'] in [0, 1, 2]:
      raise InputError('INP/NPOPT (%s) must be 0, 1, or 2!' % INPUT['inp'])
   if INPUT['cavity_surften'] < 0:
      raise InputError('CAVITY_SURFTEN must be non-negative!')
   if INPUT['fillratio'] <= 0: raise InputError('FILL_RATIO must be positive!')
   if not INPUT['radiopt'] in [0, 1]:
      raise InputError('RADIOPT (%s) must be 0 or 1!' % INPUT['radiopt'])
   if INPUT['dielc'] <= 0: raise InputError('DIELC must be positive!')
   if INPUT['maxcyc'] < 1:
      raise InputError('MAXCYC must be a positive integer!')
   if not INPUT['idecomp'] in [0, 1, 2, 3, 4]:
      raise InputError('IDECOMP (%s) must be 1, 2, 3, or 4!' % INPUT['idecomp'])
   if INPUT['idecomp'] != 0 and INPUT['sander_apbs'] == 1:
      raise InputError('IDECOMP cannot be used with sander.APBS!')
   if not INPUT['entropy'] in [0, 1]:
      raise InputError('ENTROPY (%s) must be 0 or 1!' % INPUT['entropy'])
   if not INPUT['sander_apbs'] in [0, 1]:
      raise InputError('SANDER_APBS must be 0 or 1!')
   if INPUT['alarun'] and INPUT['netcdf'] != '':
      print 'INPUT[netcdf] = %s' % INPUT['netcdf']
      raise InputError('Alanine scanning is incompatible with NETCDF != 0!')
   if INPUT['decomprun'] and INPUT['idecomp'] == 0:
      raise InputError('IDECOMP cannot be 0 for Decomposition analysis!')
   if not INPUT['use_sander'] in [0,1]:
      raise InputError('USE_SANDER must be set to 0 or 1!')
   if not INPUT['ifqnt'] in [0,1]: raise InputError('QMMM must be 0 or 1!')
   if INPUT['ifqnt'] == 1:
      if not INPUT['qm_theory'] in ['PM3', 'AM1', 'AM1_D*', 'AM1-D*', 'AM1_DH+',
                                    'AM1-DH+', 'AM1D', 'AM1_D', 'AM1/D', 'MNDO',
                                    'MNDOD', 'MNDO_D', 'MNDO/D', 'PM3-PDDG',
                                    'PM3PDDG', 'PM3_PDDG', 'PDDG-PM3',
                                    'PDDGPM3', 'PDDG_PM3', 'MNDO-PDDG',
                                    'MNDOPDDG', 'MNDO_PDDG', 'PDDG-MNDO',
                                    'PDDGMNDO', 'PM3-CARB1', 'PM3CARB1',
                                    'PM3_CARB1', 'PM3ZNB', 'PM3-ZNB', 'PM3_ZNB',
                                    'PM3/ZNB', 'ZNB', 'DFTB', 'SCCDFTB',
                                    'SCC-DFTB', 'SCC_DFTB', 'RM1', 'PM3-PDDG08',
                                    'PM3PDDG08', 'PM3_PDDG08', 'PDDG-PM308',
                                    'PDDGPM308', 'PDDG_PM308', 'PM3-PDDG-08',
                                    'PDDGPM3_08', 'PDDG_PM3_08', 'PM3-PDDG_08',
                                    'PM6', 'PM6_D', 'PM6-D', 'PM6_DH+',
                                    'PM6-DH+', 'PM3-MAIS', 'PM3MAIS',
                                    'PM3_MAIS', 'MAIS']:
         raise InputError('Invalid QM_THEORY (%s)! ' % INPUT['qm_theory'] +
                      'This variable must be set to allowable options.\n' +
                      '       See the AmberTools manual for allowable options.')
      if INPUT['qm_residues'] == '':
         raise InputError('QM_RESIDUES must be specified for IFQNT = 1!')
      if INPUT['decomprun']:
         raise InputError('QM/MM and decomposition analysis are incompatible!')
      if (INPUT['qmcharge_lig'] + INPUT['qmcharge_rec'] != INPUT['qmcharge_com']
                         and not stability):
         raise InputError('The total charge of the ligand and receptor does ' +
                           'not equal the charge of the complex!')
   if INPUT['rismrun']:
      if INPUT['rism_verbose'] > 2 or INPUT['rism_verbose'] < 0:
         raise InputError('RISM_VERBOSE must be 0, 1, or 2!')
      if INPUT['buffer'] < 0 and INPUT['solvcut'] < 0:
         raise InputError('If BUFFER < 0, SOLVCUT must be > 0!')
      if INPUT['tolerance'] < 0:
         raise InputError('TOLERANCE must be positive!')
      if INPUT['buffer'] < 0 and INPUT['ng'] == '':
         raise InputError('You must specify NG if BUFFER < 0!')
      if INPUT['closure'] == 'pse' and INPUT['closureorder'] < 1:
         raise InputError('You must specify CLOSUREORDER if CLOSURE=pse!')
      if not INPUT['polardecomp'] in [0,1]:
         raise InputError('POLARDECOMP must be either 0 or 1!')
      if not INPUT['thermo'] in ['std', 'gf', 'both']:
         raise InputError('THERMO must be "std", "gf", or "both"!')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def Abort(comm, rank, debug):
   """ Prints out useful error message (i.e. which process died) in which rank,
       then calls the MPI Abort routine to kill all processes """
   # Make sure we open any thread back up to printing to stderr
   sys.stderr = os.fdopen(sys.stderr.fileno(),'a',0)
   print >> sys.stderr, ''
   if debug: print >> sys.stderr, "MMPBSA.py failed on rank %d " % rank
   print >> sys.stderr, "Exiting. All files have been retained."
   comm.Abort(1)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def concatenate(file1, file2):
   """ Adds contents of file2 onto beginning of file1 """
   import os
   chunksize = 1048576 # Read the file in 1 MB chunks
   if type(file1).__name__ != 'file' or type(file2).__name__ != 'str':
      raise TypeError('concatenate takes an open file and a file name')
   # Open the 2 files, the first in append mode
   fl2 = open(file2, 'r')
   # Add a newline (make it OS-independent) to the first file if it doesn't
   # already end in one
   file1.write(os.linesep)

   str1 = fl2.read(chunksize)
   while str1:
      file1.write(str1)
      str1 = fl2.read(chunksize)
   
   # Now remove the merged file (file2)
   os.remove(file2)
