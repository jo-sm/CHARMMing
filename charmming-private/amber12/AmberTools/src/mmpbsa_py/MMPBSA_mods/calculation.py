"""
This module contains calculation classes that call the necessary programs
for running MM/PBSA calculations.

Last updated: 09/03/2011

######################### GPL LICENSE INFO ############################

  Copyright (C) 2009  Dwight McGee, Billy Miller III, and Jason Swails

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
   run_calculations(FILES, INPUT, rank) : Determines which calculations need to
        be run, then sets up the calculations and runs them

Classes:
   Calculation: Base calculation class
   EnergyCalculation: Typical GB/PB binding FE calculations. Handles all sander
                      and mmpbsa_py_energy program calls
   RISMCalculation: RISM binding FE calculation
   NmodeCalc: normal mode entropy calculation
   QuasiHarmCalc: Quasi-harmonic entropy calculation
"""

from MMPBSA_mods.exceptions import CalcError
import sys, os, shutil

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def run_calculations(FILES, INPUT, rank, MMPBSA_Timer, extern_progs, prefix,
                     parmsystem):
   """ Runs the calculations """

   stability = not parmsystem.receptor_prmtop and not parmsystem.ligand_prmtop

   # Determine if this is a mutant simulation. If it is, and we've also
   # calculated the normal system, then we don't actually have to run all of the
   # calculations, since some of them may be the same as the unmutated form
   mutant = 'mutant' in prefix

   if INPUT['netcdf']: trj_sfx = 'nc'
   else: trj_sfx = 'mdcrd'

   # stores which types of calculations get which program. Start with default
   # and then change them for individual cases
   progs = {'gb' : extern_progs['mmpbsa_py_energy'].full_path, 
            'pb' : extern_progs['mmpbsa_py_energy'].full_path,
            'rism' : extern_progs['rism3d.snglpnt'].full_path,
            'qh' : extern_progs['ptraj'].full_path,
            'nmode' : extern_progs['mmpbsa_py_nabnmode'].full_path}
   progs['gb'] = extern_progs['mmpbsa_py_energy'].full_path

   if INPUT['use_sander']:
      progs['gb'] = extern_progs['sander'].full_path
      progs['pb'] = extern_progs['sander'].full_path
   if INPUT['sander_apbs']:
      progs['pb'] = extern_progs['sander.APBS'].full_path
   if INPUT['ifqnt'] or INPUT['igb'] == 7 or INPUT['igb'] == 8:
      progs['gb'] = extern_progs['sander'].full_path
   if INPUT['decomprun']:
      progs['gb'] = extern_progs['sander'].full_path
      progs['pb'] = extern_progs['sander'].full_path

   master = rank == 0

   # First do the GB simulations
   if INPUT['gbrun']:

      MMPBSA_Timer.StartTimer('gb')

      print 'Beginning GB calculations with %s' % progs['gb']

      if INPUT['decomprun']: mdin = '_MMPBSA_gb_decomp_com.mdin'
      elif INPUT['ifqnt']: mdin = '_MMPBSA_gb_qmmm_com.mdin'
      else: mdin = '_MMPBSA_gb.mdin'

      if 'mmpbsa_py_energy' in progs['gb']: incrd = '%scomplex.pdb' % prefix
      else: incrd = '%sdummycomplex.inpcrd' % prefix

      # The complex is ALWAYS calculated (since mutation must be there if
      # this is the mutant calculation)
      print '  calculating complex contribution...'

      calc = EnergyCalculation(progs['gb'], parmsystem.complex_prmtop,
                   incrd, '%scomplex.%s.%d' % (prefix, trj_sfx, rank), mdin, 
                   '%scomplex_gb.mdout.%d' % (prefix, rank), 
                   '_MMPBSA_restrt.%d' % rank)
      calc.Setup()
      calc.Run()

      if not stability:
         # Don't re-run the receptor if we don't have to
         if mutant and FILES.receptor_prmtop == FILES.mutant_receptor_prmtop \
                   and not INPUT['mutant_only']:
            print '  no mutation found in receptor -- using unmutated files'
            shutil.copy('%sreceptor_gb.mdout.%d' % (
                        prefix[:prefix.index('mutant')], rank), 
                        '%sreceptor_gb.mdout.%d' % (prefix, rank))
         else:
            print '  calculating receptor contribution...'

            if INPUT['decomprun']: mdin = '_MMPBSA_gb_decomp_rec.mdin'
            elif INPUT['ifqnt']: mdin = '_MMPBSA_gb_qmmm_rec.mdin'
            else: mdin = '_MMPBSA_gb.mdin'

            if 'mmpbsa_py_energy' in progs['gb']:
               incrd = '%sreceptor.pdb' % prefix
            else: incrd = '%sdummyreceptor.inpcrd' % prefix

            calc = EnergyCalculation(progs['gb'], parmsystem.receptor_prmtop,
                     incrd, '%sreceptor.%s.%d' % (prefix, trj_sfx, rank), mdin,
                     '%sreceptor_gb.mdout.%d' % (prefix, rank),
                     '_MMPBSA_restrt.%d' % rank)
            calc.Setup()
            calc.Run()

         if mutant and FILES.ligand_prmtop == FILES.mutant_ligand_prmtop \
                   and not INPUT['mutant_only']:
            print '  no mutation found in ligand -- using unmutated files'
            shutil.copy('%sligand_gb.mdout.%d' % (
                        prefix[:prefix.index('mutant')], rank), 
                        '%sligand_gb.mdout.%d' % (prefix, rank))
         else:
            print '  calculating ligand contribution...\n'
            if INPUT['decomprun']: mdin = '_MMPBSA_gb_decomp_lig.mdin'
            elif INPUT['ifqnt']: mdin = '_MMPBSA_gb_qmmm_lig.mdin'
            else: mdin = '_MMPBSA_gb.mdin'
   
            if 'mmpbsa_py_energy' in progs['gb']: incrd = '%sligand.pdb' % prefix
            else: incrd = '%sdummyligand.inpcrd' % prefix
   
            calc = EnergyCalculation(progs['gb'], parmsystem.ligand_prmtop,
                     incrd, '%sligand.%s.%d' % (prefix, trj_sfx, rank), mdin,
                     '%sligand_gb.mdout.%d' % (prefix, rank),
                     '_MMPBSA_restrt.%d' % rank)
            calc.Setup()
            calc.Run()

      MMPBSA_Timer.StopTimer('gb')

   # end if INPUT['gbrun']

   if INPUT['pbrun']:

      MMPBSA_Timer.StartTimer('pb')

      print 'Beginning PB calculations with %s' % progs['pb']

      if INPUT['decomprun']: mdin = '_MMPBSA_pb_decomp_com.mdin'
      else: mdin = '_MMPBSA_pb.mdin'

      if 'mmpbsa_py_energy' in progs['pb']: incrd = '%scomplex.pdb' % prefix
      else: incrd = '%sdummycomplex.inpcrd' % prefix

      print '  calculating complex contribution...'

      calc = PBEnergyCalculation(progs['pb'], parmsystem.complex_prmtop,
                   incrd, '%scomplex.%s.%d' % (prefix, trj_sfx, rank), mdin, 
                   '%scomplex_pb.mdout.%d' % (prefix, rank), 
                   '_MMPBSA_restrt.%d' % rank)
      calc.Setup()
      calc.Run()

      if not stability:
         if mutant and FILES.receptor_prmtop == FILES.mutant_receptor_prmtop \
                   and not INPUT['mutant_only']:
            print '  no mutation found in receptor -- using unmutated files'
            shutil.copy('%sreceptor_pb.mdout.%d' % (
                        prefix[:prefix.index('mutant')], rank), 
                        '%sreceptor_pb.mdout.%d' % (prefix, rank))
         else:
            print '  calculating receptor contribution...'

            if INPUT['decomprun']: mdin = '_MMPBSA_pb_decomp_rec.mdin'
            else: mdin = '_MMPBSA_pb.mdin'
   
            if 'mmpbsa_py_energy' in progs['pb']: incrd = '%sreceptor.pdb' % prefix
            else: incrd = '%sdummyreceptor.inpcrd' % prefix
   
            calc = PBEnergyCalculation(progs['pb'], parmsystem.receptor_prmtop,
                     incrd, '%sreceptor.%s.%d' % (prefix, trj_sfx, rank), mdin,
                     '%sreceptor_pb.mdout.%d' % (prefix, rank),
                     '_MMPBSA_restrt.%d' % rank)
            calc.Setup()
            calc.Run()

         if mutant and FILES.ligand_prmtop == FILES.mutant_ligand_prmtop \
                   and not INPUT['mutant_only']:
            print '  no mutation found in ligand -- using unmutated files'
            shutil.copy('%sligand_pb.mdout.%d' % (
                        prefix[:prefix.index('mutant')], rank), 
                        '%sligand_pb.mdout.%d' % (prefix, rank))
         else:
            print '  calculating ligand contribution...\n'
            if INPUT['decomprun']: mdin = '_MMPBSA_pb_decomp_lig.mdin'
   
            if 'mmpbsa_py_energy' in progs['pb']: incrd = '%sligand.pdb' % prefix
            else: incrd = '%sdummyligand.inpcrd' % prefix
   
            calc = PBEnergyCalculation(progs['pb'], parmsystem.ligand_prmtop,
                     incrd, '%sligand.%s.%d' % (prefix, trj_sfx, rank), mdin,
                     '%sligand_pb.mdout.%d' % (prefix, rank),
                     '_MMPBSA_restrt.%d' % rank)
            calc.Setup()
            calc.Run()

      MMPBSA_Timer.StopTimer('pb')

   # end if INPUT['pbrun']:

   if INPUT['rismrun']:
      
      MMPBSA_Timer.StartTimer('rism')

      print 'Beginning 3D-RISM calculations with %s' % progs['rism']

      print '  calculating complex contribution...'
      calc = RISMCalculation(progs['rism'], parmsystem.complex_prmtop,
          '%scomplex.pdb' % prefix, '%scomplex.%s.%d' % (prefix, trj_sfx, rank),
          FILES.xvvfile, '%scomplex_rism.mdout.%d' % (prefix, rank), INPUT)
      calc.Setup()
      calc.Run()

      if not stability:
         if mutant and FILES.receptor_prmtop == FILES.mutant_receptor_prmtop \
                   and not INPUT['mutant_only']:
            print '  no mutation found in receptor -- using unmutated files'
            shutil.copy('%sreceptor_rism.mdout.%d' % (
                        prefix[:prefix.index('mutant')], rank), 
                        '%sreceptor_rism.mdout.%d' % (prefix, rank))
         else:
            print '  calculating receptor contribution...'
            calc = RISMCalculation(progs['rism'], parmsystem.receptor_prmtop,
                '%sreceptor.pdb' % prefix, '%sreceptor.%s.%d' % (prefix, trj_sfx,
                rank), FILES.xvvfile, '%sreceptor_rism.mdout.%d' % (prefix,rank), 
                INPUT)
            calc.Setup()
            calc.Run()
               
         if mutant and FILES.ligand_prmtop == FILES.mutant_ligand_prmtop \
                   and not INPUT['mutant_only']:
            print '  no mutation found in ligand -- using unmutated files'
            shutil.copy('%sligand_rism.mdout.%d' % (
                        prefix[:prefix.index('mutant')], rank), 
                        '%sligand_rism.mdout.%d' % (prefix, rank))
         else:
            print '  calculating ligand contribution...\n'
            calc = RISMCalculation(progs['rism'], parmsystem.ligand_prmtop,
                '%sligand.pdb' % prefix, '%sligand.%s.%d' % (prefix, trj_sfx,
                rank), FILES.xvvfile, '%sligand_rism.mdout.%d' % (prefix, rank),
                INPUT)
            calc.Setup()
            calc.Run()

      MMPBSA_Timer.StopTimer('rism')
               
   # end if INPUT['rismrun']

   if INPUT['nmoderun']:
      
      MMPBSA_Timer.StartTimer('nmode')

      print 'Beginning nmode calculations with %s' % progs['nmode']

      print '  calculating complex contribution...'
      calc = NmodeCalc(progs['nmode'], parmsystem.complex_prmtop, 
             '%scomplex.pdb' % prefix, '%scomplex_nm.%s.%d' % (prefix, trj_sfx,
             rank), '%scomplex_nm.out.%d' % (prefix, rank), INPUT)
      calc.Setup()
      calc.Run()

      if not stability:
         if mutant and FILES.receptor_prmtop == FILES.mutant_receptor_prmtop \
                   and not INPUT['mutant_only']:
            print '  no mutation found in receptor -- using unmutated files'
            shutil.copy('%sreceptor_nm.out.%d' % (
                        prefix[:prefix.index('mutant')], rank), 
                        '%sreceptor_nm.out.%d' % (prefix, rank))
         else:
            print '  calculating receptor contribution...'
            calc = NmodeCalc(progs['nmode'], parmsystem.receptor_prmtop, 
                   '%sreceptor.pdb' % prefix, '%sreceptor_nm.%s.%d' % (prefix, 
                   trj_sfx, rank), '%sreceptor_nm.out.%d' % (prefix, rank), INPUT)
            calc.Setup()
            calc.Run()

         if mutant and FILES.ligand_prmtop == FILES.mutant_ligand_prmtop \
                   and not INPUT['mutant_only']:
            print '  no mutation found in ligand -- using unmutated files'
            shutil.copy('%sligand_nm.out.%d' % (
                        prefix[:prefix.index('mutant')], rank), 
                        '%sligand_nm.out.%d' % (prefix, rank))
         else:
            print '  calculating ligand contribution...\n'
            calc = NmodeCalc(progs['nmode'], parmsystem.ligand_prmtop, 
                   '%sligand.pdb' % prefix, '%sligand_nm.%s.%d' % (prefix, 
                   trj_sfx, rank), '%sligand_nm.out.%d' % (prefix, rank), INPUT)
            calc.Setup()
            calc.Run()
      
      MMPBSA_Timer.StopTimer('nmode')

   if INPUT['entropy'] and master: # only master does quasi-harmonic calculation

      MMPBSA_Timer.StartTimer('qh')

      print 'Beginning quasi-harmonic calculation...\n'

      calc = QuasiHarmCalc(progs['qh'], parmsystem.complex_prmtop,
               '%scomplex.%s' % (prefix, trj_sfx), '%sptrajentropy.in' % prefix,
               '%sptraj_entropy.out' % prefix, INPUT['receptor_mask'],
               INPUT['ligand_mask'])
      calc.Setup()
      calc.Run()

      MMPBSA_Timer.StopTimer('qh')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Calculation(object):
   """ Base calculation class. All other calculation classes should be inherited
       from this class.
   """

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self, prog, prmtop, incrd, inptraj, input_file, output):
      self.prmtop     = str(prmtop)
      self.incrd      = incrd 
      self.input_file = input_file
      self.inptraj    = inptraj
      self.output     = output
      self.program    = prog

      self.calc_setup = False # This means that the setup has run successfully

      self.command_args = [self.program]

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Run(self, stdout=sys.stdout, stderr=sys.stderr):
      """ Runs the program. All command-line arguments must be set before 
          calling this method. Command-line arguments should be set in Setup()
      """
      from subprocess import Popen

      # If this has not been set up yet
      # then raise a stink
      if not self.calc_setup:
         raise CalcError('Cannot run a calculation without calling its' +
                        ' its Setup() function!') 

      # Here, make sure that we could pass a file *OR* a string as stdout/stderr.
      # If they are strings, then open files up with that name, and make sure to
      # close them afterwards. The Setup() method should make sure that they are
      # either a file or a string!
      if type(stdout).__name__ == 'str':
         stdout_is_string = True
         process_stdout = open(stdout, 'w', 0)
      else:
         stdout_is_string = False
         process_stdout = stdout

      if type(stderr).__name__ == 'str':
         stderr_is_string = True
         process_stderr = open(stderr, 'w', 0)
      else:
         stderr_is_string = False
         process_stderr = stderr

      # The Setup() method sets the command-line arguments and makes sure that
      # all of the CL arguments are set. Now all we have to do is start the 
      # process and monitor it for success.

      # Popen can only take strings as command-line arguments, so convert 
      # everything to a string here
      for i in range(len(self.command_args)):
         self.command_args[i] = str(self.command_args[i])

      process = Popen(self.command_args, stdin=None, stdout=process_stdout,
                      stderr=process_stderr)

      calc_failed = bool(process.wait())

      if stdout_is_string: process_stdout.close()
      if stderr_is_string: process_stderr.close()

      if calc_failed:
         raise CalcError('%s failed with prmtop %s!' % (self.program,
                         self.prmtop))

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Setup(self):
      """ Sets up the Calculation. Finds the program and adds that to the
          first element of the array. Inherited classes should call this
          method first, but then do anything else that is necessary for that
          calculation.
      """
      self.calc_setup= True

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class EnergyCalculation(Calculation):
   """ Uses mmpbsa_py_energy to evaluate energies """
   def __init__(self, prog, prmtop, incrd, inptraj, input_file, output, restrt):
      Calculation.__init__(self, prog, prmtop, incrd, inptraj, 
                           input_file, output)
      self.restrt = restrt

   def Setup(self):
      """ 
      Sets up the command-line arguments. Sander requires a unique restrt file
      for the MPI version (since one is *always* written and you don't want 2
      threads fighting to write the same dumb file)
      """
      self.command_args.append('-O')                    # overwrite flag
      self.command_args.extend(('-i', self.input_file)) # input file flag
      self.command_args.extend(('-p', self.prmtop))     # prmtop flag
      self.command_args.extend(('-c', self.incrd))      # input coordinate flag
      self.command_args.extend(('-y', self.inptraj))    # input trajectory flag
      self.command_args.extend(('-o', self.output))     # output file flag
      if self.restrt:
         self.command_args.extend(('-r', self.restrt))  # restart file flag

      # Now test to make sure that the input file exists, since that's the only
      # one that may be absent (due to the use of -use-mdins)
      if not os.path.exists(self.input_file):
         raise IOError("Input file (%s) doesn't exist" % self.input_file)
      
      self.calc_setup = True

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class RISMCalculation(Calculation):
   """ This class handles RISM calculations """

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self, prog, prmtop, incrd, inptraj, xvvfile, output, INPUT):
      """ Sets up a RISM calculation. It's not as similar to the base class as
          other calculation classes are, but it still inherits useful methods
      """
      # rism3d.snglpnt dumps its output to stdout
      Calculation.__init__(self, prog, prmtop, incrd, inptraj, None, output)

      # Set up instance variables
      self.xvvfile      = xvvfile
      self.closure      = INPUT['closure']
      self.polardecomp  = INPUT['polardecomp']
      self.ng           = INPUT['ng'].replace(' ','') # get rid of spaces
      self.solvbox      = INPUT['solvbox']
      self.closureorder = INPUT['closureorder']
      self.buffer       = INPUT['buffer']
      self.grdspc       = str(INPUT['grdspc']).replace(' ','')
      self.solvcut      = INPUT['solvcut']
      self.tolerance    = INPUT['tolerance']
      self.verbose      = INPUT['rism_verbose']
      self.solvbox      = INPUT['solvbox']

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Setup(self):
      """ Sets up the RISM calculation. All it has to do is fill in the
          necessary command-line arguments
      """
      # Set up some defaults    
      
      if self.ng == "-1,-1,-1":
         ngflag = False
      else:
         ngflag = True
   
      if self.solvbox == "-1,-1,-1":
         solvboxflag = False
      else:
         solvboxflag = True
   
      if self.polardecomp:
         polardecompflag = True
      else:
         polardecompflag = False

      Calculation.Setup(self)
      self.command_args.extend( ('--xvv', self.xvvfile, 
                                 '--closure', self.closure,
                                 '--closureorder', self.closureorder, 
                                 '--buffer', self.buffer, 
                                 '--grdspc', self.grdspc,
                                 '--solvcut', self.solvcut,
                                 '--tolerance', self.tolerance,
                                 '--verbose', self.verbose,
                                 '--prmtop', self.prmtop,
                                 '--pdb', self.incrd,
                                 '--traj', self.inptraj))
      if ngflag: self.command_args.extend(('--ng', self.ng))
      if solvboxflag: self.command_args.extend(('--solvbox', self.solvbox))
      if polardecompflag: self.command_args.extend(['--polarDecomp'])
      
      if not os.path.exists(self.xvvfile):
         raise IOError('XVVFILE (%s) does not exist!' % self.xvvfile)

      self.calc_setup = True

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Run(self):
      Calculation.Run(self, stdout=self.output)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class NmodeCalc(Calculation):
   """ Calculates entropy contribution by normal mode approximation """
   
   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self, prog, prmtop, incrd, inptraj, output, INPUT):
      """ Initializes the nmode calculation. Need to set the options string """
      from math import sqrt
      Calculation.__init__(self, prog, prmtop, incrd, inptraj, None, output)
      
      kappa = sqrt(0.10806 * INPUT['nmode_istrng'])
      if INPUT['nmode_igb']:
         option_string = ('ntpr=10000, diel=C, kappa=%f, cut=1000, gb=1, ' +
                          'dielc=%f, temp0=%f') % (kappa, 
                          INPUT['dielc'], INPUT['temp'])
      else:
         option_string = ('ntpr=10000, diel=R, kappa=%f, cut=1000, gb=0, ' +
                          'dielc=%f, temp0=%f') % (kappa, 
                          INPUT['dielc'], INPUT['temp'])

      self.option_string = option_string
      self.drms = INPUT['drms']
      self.maxcyc = INPUT['maxcyc']

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Setup(self):
      """ Sets up the simulation """

      self.command_args.extend((self.incrd, self.prmtop, self.maxcyc, self.drms,
                                 self.option_string, self.inptraj))
      self.calc_setup = True

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Run(self):
      Calculation.Run(self, stdout=self.output)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class QuasiHarmCalc(Calculation):
   """ Quasi-harmonic entropy calculation class """

   def __init__(self, prog, prmtop, inptraj, input_file, output,
                receptor_mask, ligand_mask):
      """ Initializes the Quasi-harmonic calculation class """
      Calculation.__init__(self, prog, prmtop, None, inptraj, 
                           input_file, output)
      self.stability = not bool(receptor_mask) and not bool(ligand_mask)
      self.receptor_mask, self.ligand_mask = receptor_mask, ligand_mask
      self.calc_setup = False

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
   
   def Setup(self):
      """ Sets up a Quasi-harmonic calculation """
      from subprocess import Popen, PIPE

      # Determine the prefix from our input file... hack way to do this
      if self.input_file.startswith('_MMPBSA_mutant_'): prefix='_MMPBSA_mutant_'
      else: prefix='_MMPBSA_'

      # Make sure masks are a list, and that there are enough masks

      # First thing we need is the average PDB as a reference
      ptraj_str = 'trajin %s\naverage %savgcomplex.pdb pdb' % (self.inptraj,
                  prefix)

      outfile = open('_MMPBSA_create_average.out','w',0)

      process = Popen([self.program, self.prmtop], stdin=PIPE, stdout=outfile)
      process.communicate(ptraj_str)

      if process.wait():
         raise CalcError('Failed creating average PDB')

      outfile.close()

      # Now that we have the PDB file
      
      self.command_args.extend((self.prmtop, self.input_file))

      self.calc_setup = True

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
   
   def Run(self):
      Calculation.Run(self, stdout=self.output)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class PBEnergyCalculation(EnergyCalculation):
   """
   Specially handle the PB calculations to extract warnings and errors PBSA
   prints to stdout and redirect them to the user
   """
   def Run(self, stderr=sys.stderr):
      """ Runs the program. All command-line arguments must be set before 
          calling this method. Command-line arguments should be set in Setup()
      """
      import re
      from subprocess import Popen, PIPE

      # If this has not been set up yet
      # then raise a stink
      if not self.calc_setup:
         raise CalcError('Cannot run a calculation without calling its' +
                        ' its Setup() function!') 

      errorre = re.compile('(pb (?:bomb)|(?:warning))', re.I)
      # Here, make sure that we could pass a file *OR* a string as stderr.
      if type(stderr).__name__ == 'str':
         stderr_is_string = True
         process_stderr = open(stderr, 'w', 0)
      else:
         stderr_is_string = False
         process_stderr = stderr

      # The Setup() method sets the command-line arguments and makes sure that
      # all of the CL arguments are set. Now all we have to do is start the 
      # process and monitor it for success.

      # Popen can only take strings as command-line arguments, so convert 
      # everything to a string here
      for i in range(len(self.command_args)):
         self.command_args[i] = str(self.command_args[i])

      process = Popen(self.command_args, stdin=None, stdout=PIPE,
                      stderr=process_stderr)

      out, err = process.communicate('')
      calc_failed = bool(process.wait())

      if stderr_is_string: process_stderr.close()

      if calc_failed:
         error_list = [s.strip() for s in out.split('\n') 
                                       if errorre.match(s.strip())]
         raise CalcError('%s failed with prmtop %s!\n\t' % (self.program,
                         self.prmtop) + '\n\t'.join(error_list) + '\n')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
