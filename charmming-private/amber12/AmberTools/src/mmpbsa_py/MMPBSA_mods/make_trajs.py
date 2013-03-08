"""
This module contains classes and such that are responsible for generating
trajectories for MM/PBSA based on input, starting trajectories. It will
strip trajectories and output the ones that will be used for the calculations.
Needed for proper operation of MMPBSA.py

Last updated:  9/03/2011

                          GPL LICENSE INFO                             

  Copyright (C) 2009 - 2011  Dwight McGee, Billy Miller III, and Jason Swails

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


Methods:
         make_trajectories(INPUT, FILES, size): Makes all of the non-mutant
            trajectories needed for the calculation, as well as all dummy
            files (i.e. restarts and PDBs)

         make_mutant_trajectories(INPUT, FILES, rank): Mutates the trajectories

Classes:
         Trajectory: Class for manipulating Amber trajectories through cpptraj
"""
from MMPBSA_mods.exceptions import TrajError, MMPBSA_Error

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def make_trajectories(INPUT, FILES, size, cpptraj):
   """ 
   This function creates the necessary trajectory files, and creates thread-
   specific trajectories for parallel calculations
   """
   # If we have a solvated_prmtop, set up our trajectory and strip the solvent
   stability = not FILES.receptor_prmtop and not FILES.ligand_prmtop

   # File suffix is dependent on file type
   if INPUT['netcdf']: trj_suffix = 'nc'
   else: trj_suffix = 'mdcrd'

   # Make sure we set up for a solvated topology; image, strip solvent, etc.
   if FILES.solvated_prmtop:
      traj = Trajectory(FILES.solvated_prmtop, FILES.mdcrd, cpptraj)
      traj.Setup(INPUT['startframe'], INPUT['endframe'], INPUT['interval'])
      if not stability: 
         traj.Image(INPUT['receptor_mask'])
      traj.StripSolvent(INPUT['strip_mask'])
   # If we do not have a solvated topology...
   else:
      traj = Trajectory(FILES.complex_prmtop, FILES.mdcrd, cpptraj)
      traj.Setup(INPUT['startframe'], INPUT['endframe'], INPUT['interval'])
   # RMS fit
   traj.rms('!(%s)' % INPUT['strip_mask'])

   num_frames_total = traj.processed_frames
   num_frames_nmode = 0

   # Sanity check
   if traj.processed_frames < size: 
      raise MMPBSA_Error('Must have at least as many frames as processors!')

   # We now know how many frames we have in total, so make a list that lists the
   # number of frames found for each rank, and assign extra frames incrementally
   frames_per_rank = traj.processed_frames // size
   extras = traj.processed_frames - frames_per_rank * size
   frame_count = [frames_per_rank for i in range(size)]
   for i in range(size):
      if i < extras: frame_count[i] += 1

   # Dump our complex trajectories
   if INPUT['full_traj'] or INPUT['entropy']:
      traj.Outtraj('_MMPBSA_complex.%s' % trj_suffix, filetype=INPUT['netcdf'])
   traj.Outtraj('_MMPBSA_complex.pdb', frames='1', filetype='pdb')
   traj.Outtraj('_MMPBSA_dummycomplex.inpcrd',frames='1',filetype='restart')
   
   # Now dump thread-specific trajectories
   last_frame = 1
   for i in range(size):
      frame_string = '%d-%d' % (last_frame, last_frame+frame_count[i]-1)
      traj.Outtraj('_MMPBSA_complex.%s.%d' % (trj_suffix, i), 
                   frames=frame_string, filetype=INPUT['netcdf'])
      last_frame += frame_count[i]
   
   # Now create the receptor/ligand trajectories if we're taking them from
   # the complex trajectory

   if not stability and not FILES.receptor_mdcrd:
      traj.Strip(INPUT['ligand_mask'])
      if INPUT['full_traj'] or INPUT['entropy']:
         traj.Outtraj('_MMPBSA_receptor.%s'%trj_suffix,filetype=INPUT['netcdf'])
      traj.Outtraj('_MMPBSA_receptor.pdb', frames='1', filetype='pdb')
      traj.Outtraj('_MMPBSA_dummyreceptor.inpcrd',frames='1',filetype='restart')
      last_frame = 1
      for i in range(size):
         frame_string = '%d-%d' % (last_frame, last_frame+frame_count[i]-1)
         traj.Outtraj('_MMPBSA_receptor.%s.%d' % (trj_suffix, i),
                      frames=frame_string, filetype=INPUT['netcdf'])
         last_frame += frame_count[i]
      traj.Unstrip(restrip_solvent=True)
      traj.rms('!(%s)' % INPUT['strip_mask'])

   if not stability and not FILES.ligand_mdcrd:
      traj.Strip(INPUT['receptor_mask'])
      if INPUT['full_traj'] or INPUT['entropy']:
         traj.Outtraj('_MMPBSA_ligand.%s' % trj_suffix,filetype=INPUT['netcdf'])
      traj.Outtraj('_MMPBSA_ligand.pdb', frames='1', filetype='pdb')
      traj.Outtraj('_MMPBSA_dummyligand.inpcrd', frames='1', filetype='restart')
      last_frame = 1
      for i in range(size):
         frame_string = '%d-%d' % (last_frame, last_frame+frame_count[i]-1)
         traj.Outtraj('_MMPBSA_ligand.%s.%d' % (trj_suffix, i),
                      frames=frame_string, filetype=INPUT['netcdf'])
         last_frame += frame_count[i]
      traj.Unstrip(restrip_solvent=True)
      traj.rms('!(%s)' % INPUT['strip_mask'])

   # Run cpptraj to get the trajectory
   traj.Run('_MMPBSA_normal_traj_cpptraj.out')

   # Go back and do the receptor and ligand if we used a multiple 
   # trajectory approach
   if not stability and FILES.receptor_mdcrd:
      # Different actions depending on solvation...
      if FILES.solvated_receptor_prmtop:
         rectraj = Trajectory(FILES.solvated_receptor_prmtop, 
                              FILES.receptor_mdcrd, cpptraj)
         rectraj.Setup(INPUT['startframe'],INPUT['endframe'],INPUT['interval'])
         rectraj.StripSolvent(INPUT['strip_mask'])
      else:
         rectraj = Trajectory(FILES.receptor_prmtop, 
                              FILES.receptor_mdcrd, cpptraj)
         rectraj.Setup(INPUT['startframe'],INPUT['endframe'],INPUT['interval'])
      rectraj.rms('!(%s)' % INPUT['strip_mask'])
      if INPUT['full_traj'] or INPUT['entropy']:
         rectraj.Outtraj('_MMPBSA_receptor.%s' % trj_suffix, 
                         filetype=INPUT['netcdf'])
      rectraj.Outtraj('_MMPBSA_receptor.pdb', frames='1', filetype='pdb')
      rectraj.Outtraj('_MMPBSA_dummyreceptor.inpcrd',frames='1',
                      filetype='restart')

      # Now do the same split-up of workload as we did for complex, but don't
      # assume the same number of frames as we had for the complex
      if rectraj.processed_frames < size: 
         raise MMPBSA_Error('Too many procs for receptor snapshots')
      frames_per_rank = rectraj.processed_frames // size
      extras = rectraj.processed_frames - frames_per_rank * size
      frame_count = [frames_per_rank for i in range(size)]
      for i in range(size):
         if i < extras: frame_count[i] += 1
      last_frame = 1
      for i in range(size):
         frame_string = '%d-%d' % (last_frame, last_frame+frame_count[i]-1)
         rectraj.Outtraj('_MMPBSA_receptor.%s.%d' % (trj_suffix, i),
                         frames=frame_string, filetype=INPUT['netcdf'])
         last_frame += frame_count[i]

      rectraj.Run('_MMPBSA_receptor_traj_cpptraj.out')

   # end if not stability and FILES.receptor_mdcrd

   if not stability and FILES.ligand_mdcrd:
      # Different actions depending on solvation...
      if FILES.solvated_ligand_prmtop:
         ligtraj = Trajectory(FILES.solvated_ligand_prmtop, 
                              FILES.ligand_mdcrd, cpptraj)
         ligtraj.Setup(INPUT['startframe'],INPUT['endframe'],INPUT['interval'])
         ligtraj.StripSolvent(INPUT['strip_mask'])
      else:
         ligtraj = Trajectory(FILES.ligand_prmtop, FILES.ligand_mdcrd, cpptraj)
         ligtraj.Setup(INPUT['startframe'],INPUT['endframe'],INPUT['interval'])
      ligtraj.rms('!(%s)' % INPUT['strip_mask'])
      ligtraj.Outtraj('_MMPBSA_ligand.%s'%trj_suffix,filetype=INPUT['netcdf'])
      ligtraj.Outtraj('_MMPBSA_ligand.pdb', frames='1', filetype='pdb')
      ligtraj.Outtraj('_MMPBSA_dummyligand.inpcrd',frames='1',
                      filetype='restart')

      # Now do the same split-up of workload as we did for complex, but don't
      # assume the same number of frames as we had for the complex
      if ligtraj.processed_frames < size: 
         raise MMPBSA_Error('Too many procs for ligand snapshots')
      frames_per_rank = ligtraj.processed_frames // size
      extras = ligtraj.processed_frames - frames_per_rank * size
      frame_count = [frames_per_rank for i in range(size)]
      for i in range(size):
         if i < extras: frame_count[i] += 1
      last_frame = 1
      for i in range(size):
         frame_string = '%d-%d' % (last_frame, last_frame+frame_count[i]-1)
         ligtraj.Outtraj('_MMPBSA_ligand.%s.%d' % (trj_suffix, i),
                         frames=frame_string, filetype=INPUT['netcdf'])
         last_frame += frame_count[i]

      ligtraj.Run('_MMPBSA_ligand_traj_cpptraj.out')

   # end if not stability and FILES.ligand_mdcrd

   # Now make the nmode trajectories
   if INPUT['nmoderun']:
      nmtraj = Trajectory(FILES.complex_prmtop, ['_MMPBSA_complex.%s.%d' %
               (trj_suffix, i) for i in range(size)], cpptraj)
      nmtraj.Setup(INPUT['nmstartframe'], INPUT['nmendframe'], 
                   INPUT['nminterval'])

      num_frames_nmode = nmtraj.processed_frames

      # Now split up the complex trajectory by thread
      if nmtraj.processed_frames < size:
         raise MMPBSA_Error('More processors than complex nmode frames!')

      frames_per_rank = nmtraj.processed_frames // size
      extras = nmtraj.processed_frames - frames_per_rank * size
      frame_count = [frames_per_rank for i in range(size)]
      for i in range(size):
         if i < extras: frame_count[i] += 1
      last_frame = 1
      for i in range(size):
         frame_string = '%d-%d' % (last_frame, last_frame+frame_count[i]-1)
         nmtraj.Outtraj('_MMPBSA_complex_nm.%s.%d' % (trj_suffix, i),
                        frames=frame_string, filetype=INPUT['netcdf'])
         last_frame += frame_count[i]

      nmtraj.Run('_MMPBSA_com_nm_traj_cpptraj.out')

      if not stability:
         nmtraj = Trajectory(FILES.receptor_prmtop, ['_MMPBSA_receptor.%s.%d' %
                  (trj_suffix, i) for i in range(size)], cpptraj)
         nmtraj.Setup(INPUT['nmstartframe'], INPUT['nmendframe'], 
                      INPUT['nminterval'])

         # Now split up the complex trajectory by thread
         if nmtraj.processed_frames < size:
            raise MMPBSA_Error('More processors than receptor nmode frames!')

         frames_per_rank = nmtraj.processed_frames // size
         extras = nmtraj.processed_frames - frames_per_rank * size
         frame_count = [frames_per_rank for i in range(size)]
         for i in range(size):
            if i < extras: frame_count[i] += 1
         last_frame = 1
         for i in range(size):
            frame_string = '%d-%d' % (last_frame, last_frame+frame_count[i]-1)
            nmtraj.Outtraj('_MMPBSA_receptor_nm.%s.%d' % (trj_suffix, i),
                           frames=frame_string, filetype=INPUT['netcdf'])
            last_frame += frame_count[i]

         nmtraj.Run('_MMPBSA_rec_nm_traj_cpptraj.out')

         nmtraj = Trajectory(FILES.ligand_prmtop, ['_MMPBSA_ligand.%s.%d' %
                  (trj_suffix, i) for i in range(size)], cpptraj)
         nmtraj.Setup(INPUT['nmstartframe'], INPUT['nmendframe'], 
                      INPUT['nminterval'])

         # Now split up the complex trajectory by thread
         if nmtraj.processed_frames < size:
            raise MMPBSA_Error('More processors than ligand nmode frames!')

         frames_per_rank = nmtraj.processed_frames // size
         extras = nmtraj.processed_frames - frames_per_rank * size
         frame_count = [frames_per_rank for i in range(size)]
         for i in range(size):
            if i < extras: frame_count[i] += 1
         last_frame = 1
         for i in range(size):
            frame_string = '%d-%d' % (last_frame, last_frame+frame_count[i]-1)
            nmtraj.Outtraj('_MMPBSA_ligand_nm.%s.%d' % (trj_suffix, i),
                           frames=frame_string, filetype=INPUT['netcdf'])
            last_frame += frame_count[i]

         nmtraj.Run('_MMPBSA_lig_nm_traj_cpptraj.out')

      # end if not stability
   
   # end if INPUT['nmoderun']

   return num_frames_total, num_frames_nmode

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def make_mutant_trajectories(INPUT, FILES, rank, cpptraj, norm_sys, mut_sys):
   """ Mutates given trajectories and outputs dummy files for mutants """
   from alamdcrd import MutantMdcrd
   import shutil
   if not INPUT['alarun']: return None, None

   stability = not FILES.receptor_prmtop and not FILES.ligand_prmtop

   if INPUT['netcdf']:
      trj_suffix = 'nc'
      raise TypeError('Alanine scanning requires ASCII trajectories (netcdf=0)')
   else:
      trj_suffix = 'mdcrd'

   if not stability and not ((FILES.ligand_prmtop == FILES.mutant_ligand_prmtop
         and FILES.receptor_prmtop != FILES.mutant_receptor_prmtop) or (
         FILES.ligand_prmtop != FILES.mutant_ligand_prmtop and
         FILES.receptor_prmtop == FILES.mutant_receptor_prmtop)):
      raise MMPBSA_Error('Alanine scanning requires either a mutated ligand ' +
         'or receptor topology file with only 1 mutant residue, but not both')

   master = rank == 0

   # Have each rank mutate our rank's normal complex trajectory
   com_mut = MutantMdcrd('_MMPBSA_complex.%s.%d' % (trj_suffix, rank),
                         norm_sys.complex_prmtop, mut_sys.complex_prmtop)
   com_mut.MutateTraj('_MMPBSA_mutant_complex.%s.%d' % (trj_suffix, rank))

   # Have each rank mutate our rank's normal receptor or ligand trajectory
   # and copy the normal one to the mutant if the mutated residue is *not*
   # present in there
   if not stability:
      if FILES.receptor_prmtop != FILES.mutant_receptor_prmtop:
         rec_mut = MutantMdcrd('_MMPBSA_receptor.%s.%d' % (trj_suffix, rank),
                            norm_sys.receptor_prmtop, mut_sys.receptor_prmtop)
         rec_mut.MutateTraj('_MMPBSA_mutant_receptor.%s.%d' % (trj_suffix,rank))
         shutil.copyfile('_MMPBSA_ligand.%s.%d' % (trj_suffix, rank),
                         '_MMPBSA_mutant_ligand.%s.%d' % (trj_suffix, rank))

      elif FILES.ligand_prmtop != FILES.mutant_ligand_prmtop:
         lig_mut = MutantMdcrd('_MMPBSA_ligand.%s.%d' % (trj_suffix, rank),
                            norm_sys.ligand_prmtop, mut_sys.ligand_prmtop)
         lig_mut.MutateTraj('_MMPBSA_mutant_ligand.%s.%d' % (trj_suffix,rank))
         shutil.copyfile('_MMPBSA_receptor.%s.%d' % (trj_suffix, rank),
                         '_MMPBSA_mutant_receptor.%s.%d' % (trj_suffix, rank))
   
   # Have our master dump out dummy files
   if master:
      com_traj = Trajectory(FILES.mutant_complex_prmtop, 
                            '_MMPBSA_mutant_complex.%s.0' % trj_suffix, cpptraj)
      com_traj.Setup(1,1,1)
      com_traj.Outtraj('_MMPBSA_mutant_complex.pdb', frames='1', filetype='pdb')
      com_traj.Outtraj('_MMPBSA_mutant_dummycomplex.inpcrd' , frames='1',
                       filetype='restart')
      com_traj.Run('_MMPBSA_mutant_complex_cpptraj.out')
      if not stability:
         rec_traj = Trajectory(FILES.mutant_receptor_prmtop,
                           '_MMPBSA_mutant_receptor.%s.0' % trj_suffix, cpptraj)
         rec_traj.Setup(1,1,1)
         rec_traj.Outtraj('_MMPBSA_mutant_receptor.pdb', frames='1', 
                          filetype='pdb')
         rec_traj.Outtraj('_MMPBSA_mutant_dummyreceptor.inpcrd' , frames='1',
                          filetype='restart')
         rec_traj.Run('_MMPBSA_mutant_receptor_cpptraj.out')

         lig_traj = Trajectory(FILES.mutant_ligand_prmtop,
                             '_MMPBSA_mutant_ligand.%s.0' % trj_suffix, cpptraj)
         lig_traj.Setup(1,1,1)
         lig_traj.Outtraj('_MMPBSA_mutant_ligand.pdb', frames='1', 
                          filetype='pdb')
         lig_traj.Outtraj('_MMPBSA_mutant_dummyligand.inpcrd' , frames='1',
                          filetype='restart')
         lig_traj.Run('_MMPBSA_mutant_ligand_cpptraj.out')

   # Mutate our nmode trajectories if need be
   if INPUT['nmoderun']:
      com_mut = MutantMdcrd('_MMPBSA_complex_nm.%s.%d' % (trj_suffix, rank),
                            norm_sys.complex_prmtop, mut_sys.complex_prmtop)
      com_mut.MutateTraj('_MMPBSA_mutant_complex_nm.%s.%d' % (trj_suffix, rank))
      if not stability and FILES.receptor_prmtop !=FILES.mutant_receptor_prmtop:
         rec_mut = MutantMdcrd('_MMPBSA_receptor_nm.%s.%d' % (trj_suffix, rank),
                            norm_sys.receptor_prmtop, mut_sys.receptor_prmtop)
         rec_mut.MutateTraj('_MMPBSA_mutant_receptor_nm.%s.%d' % 
                            (trj_suffix, rank))
         shutil.copyfile('_MMPBSA_ligand_nm.%s.%d' % (trj_suffix, rank),
                         '_MMPBSA_mutant_ligand_nm.%s.%d' % (trj_suffix, rank))

      if not stability and FILES.ligand_prmtop != FILES.mutant_ligand_prmtop:
         lig_mut = MutantMdcrd('_MMPBSA_ligand_nm.%s.%d' % (trj_suffix, rank),
                            norm_sys.ligand_prmtop, mut_sys.ligand_prmtop)
         lig_mut.MutateTraj('_MMPBSA_mutant_ligand_nm.%s.%d' % 
                            (trj_suffix, rank))
         shutil.copyfile('_MMPBSA_ligand_nm.%s.%d' % (trj_suffix, rank),
                         '_MMPBSA_mutant_ligand_nm.%s.%d' % (trj_suffix, rank))

   # If we're doing a quasi-harmonic approximation we need the full com traj
   if (INPUT['full_traj'] or INPUT['entropy']) and master:
      com_mut = MutantMdcrd('_MMPBSA_complex.%s' % trj_suffix,
                            norm_sys.complex_prmtop, mut_sys.complex_prmtop)
      com_mut.MutateTraj('_MMPBSA_mutant_complex.%s' % trj_suffix)

   return str(com_mut), com_mut.mutres

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Trajectory(object):
   """ Base Trajectory class:
       Methods: __init__(prmtop, traj_files)
                Setup(startframe, endframe, interval)
                Query()
                Image(masks)
                Strip(mask)
                StripSolvent(strip_mask)
                Unstrip(restrip_solvent=True)
                Outtraj(fname, startframe, endframe, interval, filetype, nobox)
                Run(output)

       The way this class is typically used is initializing it by giving it a
       list of trajectory files. It will then automatically run Query to make
       sure that all trajectory files are valid, and it counts how many frames
       are present in each trajectory. Then you can strip the solvent (this
       just tags what is the solvent in case you want to re-strip it after every
       Unstrip() command). Then, add masks to strip out. Then, you can run with 
       the current actions via the Run() routine. Typical sequence:

       trajsystem = Trajectory(solvated_prmtop, mdcrds)
       trajsystem.Setup(startframe, endframe, interval)
       trajsystem.Image([receptor_mask])
       trajsystem.StripSolvent(strip_mask)
       trajsystem.Strip(receptor_mask)
       trajsystem.Outtraj('_MMPBSA_ligand.mdcrd', filetype=INPUT['netcdf'])
       trajsystem.Unstrip(restrip_solvent=True)
       trajsystem.Strip(ligand_mask)
       trajsystem.Outtraj('_MMPBSA_receptor.mdcrd', filetype=INPUT['netcdf'])
       trajsystem.Run('_MMPBSA_create_trajectories.out')
   """

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self, prmtop, traj_files, cpptraj='cpptraj'):
      """ Sets up a basic trajectory """

      # Make sure we have a list or string, and force it to be a list
      if type(traj_files).__name__ == 'str':
         self.traj_files = [traj_files]
      elif type(traj_files).__name__ == 'list':
         self.traj_files = traj_files
      else:
         raise TrajError('Trajectory only takes a string or list!')

      if type(prmtop).__name__ != 'str':
         raise TrajError('Trajectory expects the prmtop name only!')
        
      self.prmtop = prmtop

      # Find cpptraj
      self.exe = cpptraj

      self.strip_solvent = False

      self.Query()

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Setup(self, startframe=1, endframe=9999999, interval=1):
      """ Finds the program and clears the queue and other attributes assigned
          during other function calls
      """

      self.strip_solvent = False
      orig_endframe = endframe
      orig_startframe = startframe

      # Delete the stripmask, but make its non-existence non-fatal, since it
      # won't exist unless Setup is called a second time.

      try: del self.stripmask
      except AttributeError: pass

      self.actions = []

      # Catch stupid choices
      if startframe > endframe:
         raise TrajError('Cannot have startframe (%d) > endframe (%d)' %
                         (startframe, endframe))

      if startframe < 0:
         raise TrajError('Startframe (%d) < 0' % startframe)

      if interval <= 0:
         raise TrajError('Interval (%d) <= 0' % interval)

      if startframe > self.total_frames:
         raise TrajError('start frame (%d) > total frames (%d)' % 
                         (startframe, self.total_frames))

      # Bring endframe down if it's > the total number of frames

      endframe = min(endframe, self.total_frames)

      self.analyzed_frames = int((endframe - startframe) / interval) + 1

      # If we have an interval != 1, then it's possible that the endframe that
      # the user set is not the actual one that will be used. To make things
      # simpler, I will adjust the endframe to what is *actually* used

      endframe = startframe + interval*int((endframe - startframe) // interval)

      # Set up start and end arrays (for each trajectory)

      traj_starts = [-1 for i in range(len(self.traj_files))]
      traj_ends   = [0 for i in range(len(self.traj_files))]

      # Now determine where each trajectory starts and ends

      for i in range(len(self.traj_files)):

         # skip over any trajectories that lie entirely before startframe

         if startframe > self.traj_sizes[i]:
            traj_starts[i] = -1 # this will tell us to skip this traj file
            startframe -= self.traj_sizes[i]
            continue

         # Now we start at our startframe

         traj_starts[i] = startframe

         # Now we figure out our last frame, adjusting for interval

         last_frame = startframe + interval * int((self.traj_sizes[i]
                     - startframe) / interval)

         traj_ends[i] = min(endframe, last_frame, self.traj_sizes[i])

         # Now determine our new start frame for the next trajectory. First
         # we need to know how many frames are left over at the end of the
         # last frame, and subtract those from the interval to determine
         # where our next trajectory should start from. Also move our endframe
         # down, but only by our last_frame! (since the last couple frames in
         # our trajectory file could be leftover due to interval > 1)

         startframe = interval + last_frame - self.traj_sizes[i]

         endframe -= last_frame

         if endframe < self.traj_sizes[i]:
            break # we're done now

      # end for i len(traj_files)

      # We now have to trajin each of the analyzed trajectories

      for i in range(len(self.traj_files)):
         if traj_starts[i] < 0: continue # skip -1's

         self.actions.append('trajin %s %d %d %d' % (self.traj_files[i],
                             traj_starts[i], traj_ends[i], interval))

      self.actions.append('noprogress') # quash the progress bar

      self.processed_frames = (min(orig_endframe, self.total_frames) - 
                               orig_startframe) / interval + 1

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def StripSolvent(self, stripmask):
      """ Strips the solvent """
      self.stripmask = stripmask
      self.strip_solvent = True
      self.Strip(stripmask)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Query(self):
      """ Finds out how many frames are in the given traj files """
      from subprocess import Popen, PIPE

      self.traj_sizes = []

      # Determine how many frames are in each trajectory of the list

      for traj in self.traj_files:
         
         process = Popen([self.exe, self.prmtop], stdin=PIPE, stdout=PIPE)
         
         (output, error) = process.communicate('trajin %s\n' % traj)
         
         if process.wait(): # if it quits with return code != 0
            raise TrajError('%s failed when querying %s' % (self.exe, traj))

         # Now parse the output to find out how many frames are there. We are
         # looking for "Coordinate processing will occur on x frames."

         outputlines = output.split('\n')
         for line in outputlines:
            line = line.strip()
            if line.startswith('Coordinate processing will occur'):
               num_frames = int(line.split()[5])
               if num_frames == 0:
                  raise TrajError('Trajectory %s has 0 frames!' % traj)
               self.traj_sizes.append(int(line.split()[5]))
               break

      # end for traj in self.traj_files

      self.total_frames = sum(self.traj_sizes)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Image(self, masks):
      """ Images the trajectories """
      from chemistry.amber.readparm import amberParm

      if type(masks).__name__ == 'str':
         masks = [masks]

      elif type(masks).__name__ == 'tuple':
         masks = list(masks)

      elif type(masks).__name__ != 'list':
         raise TrajError('Image() expects a list of masks or a string!')

      solvated_prmtop = amberParm(self.prmtop)
      if not solvated_prmtop.valid:
         raise TrajError('Invalid solvated prmtop file %s' % self.prmtop)

      ifbox = solvated_prmtop.ptr('ifbox')
      if ifbox == 0:
         raise TrajError('Solvated topology %s has IFBOX == 0' % ifbox)

      for mask in masks:
         self.actions.append('center %s mass origin' % mask)
         if ifbox == 2: # truncated octahedron requires familiar keyword
            self.actions.append('image origin familiar')
         else:
            self.actions.append('image origin')

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def rms(self, mask):
      """ Does an RMS fit around a specific mask """
      self.actions.append('rmsd %s mass first' % mask)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Strip(self, mask):
      """ Strips a mask from the coordinates. """
      self.actions.append('strip %s' % mask)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
   
   def Unstrip(self, restrip_solvent=True):
      """ Returns to unstripped state """
      self.actions.append('unstrip')
      if self.strip_solvent and restrip_solvent:
         self.actions.append('strip %s' % self.stripmask)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Run(self, output_file=None):
      """ Runs cpptraj to actually create the files """
      from sys import stdout as sys_stdout
      from subprocess import Popen, PIPE

      # Accept output_file as both a file object and string object
      if type(output_file).__name__ == 'str':
         stdout = open(output_file, 'w', 0)

      elif type(output_file).__name__ == 'file':
         stdout = output_file

      elif output_file == None:
         stdout = sys_stdout

      else:
         raise TrajError(('%s not allowed type for output  file in' + 
                          'Trajectory.Run()') % type(output_file).__name__)

      # Now it's time to run the program
      
      input_string = ''
      for action in self.actions:
         input_string += action.strip() + '\n'

      process = Popen([self.exe, self.prmtop], stdout=stdout, stdin=PIPE)

      process.communicate(input_string)

      if process.wait():
         raise TrajError('Error running %s' % self.program)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

   def Outtraj(self, filename, frames=None, filetype='', nobox='nobox'):
      """ This adds an outtraj command to the action stack, and you can specify
          the type of trajectory file to output (such as restart/pdb for input
          files, etc.)
      """
      if not frames: frames = '1-%d' % self.total_frames
      self.actions.append('outtraj %s onlyframes %s %s %s' % (filename, frames,
                          nobox, filetype))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
