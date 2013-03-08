#! /usr/bin/python2.6

"""
MMPBSA.py is a script for performing (M)olecular (M)echanics          
(P)oisson (B)oltzmann (S)urface (A)rea to find free energies of       
binding. Refer to the AMBER manual and/or relevant literature for a   
more thorough overview of the method. This implementation uses AMBER  
executables to find energies using either Poisson Boltzmann or        
Generalized Born implicit solvent models of a complex of a receptor   
with a bound ligand. This script was written by Dwight McGee,         
Billy Miller III, and Jason Swails in Adrian Roitberg's research      
group at the Quantum Theory Project at the University of Florida.     
                                                                        
Last updated: 10/31/2011                                        

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
"""
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#
# BEGIN import system and MMPBSA-specific modules, unbuffer std output streams,
#       and set up MPI universe
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

# Import system modules
import sys, os, signal, re

# Unbuffer sys.stdout and sys.stderr
sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
sys.stderr = os.fdopen(sys.stderr.fileno(), 'w', 0)

from MMPBSA_mods.fake_mpi import MPI

rank = 0
size = 1
debug_printlevel = 0

# Set up the new excepthook to control how fatal exceptions are handled

def excepthook(exception_type, exception_value, tb):
   """ Replaces sys.excepthook so fatal exceptions kill all MPI threads and
       we can control the printing of tracebacks. Those are helpful for 
       debugging purposes, but may be unsightly to users. debug_printlevel
       set above controls this behavior
   """
   import traceback
   if debug_printlevel > 1: traceback.print_tb(tb)
   sys.stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))
   if size > 1: 
      sys.stderr.write('Error occured on rank %d.' % rank + os.linesep)
   sys.stderr.write('Exiting. All files have been retained.' + os.linesep)
   MPI.COMM_WORLD.Abort(1)

sys.excepthook = excepthook # replace sys.excepthook with my definition above

# Set up the SIGINT handler

def interrupt_handler(signal, frame):
   """ Handles interrupt signals for a clean exit """
   print >> sys.stderr, ''
   print >> sys.stderr, ('%s interrupted! Program terminated. ' +
            'All files are kept.') % os.path.split(sys.argv[0])[1]
   MPI.COMM_WORLD.Abort(1)

signal.signal(signal.SIGINT, interrupt_handler)

# Import the rest of the MMPBSA modules
from MMPBSA_mods import utils
from MMPBSA_mods.calculation import run_calculations
from MMPBSA_mods.commandline_parser import OptionParser
from MMPBSA_mods.createinput import create_inputs
from MMPBSA_mods.exceptions import MMPBSA_Error
from MMPBSA_mods.findprogs import find_progs
from MMPBSA_mods.input_parser import InputFile
from MMPBSA_mods.make_trajs import make_trajectories, make_mutant_trajectories
from MMPBSA_mods.output_file import write_stability_output, write_binding_output
from MMPBSA_mods.output_file import write_decomp_stability_output
from MMPBSA_mods.output_file import write_decomp_binding_output
from MMPBSA_mods.parm_setup import MMPBSA_System
from MMPBSA_mods.timer import Timer

# Now make sure we have AMBERHOME set
if not os.getenv('AMBERHOME'):
   raise MMPBSA_Error('AMBERHOME environment variable is not set!')

rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()
master = rank == 0

# Set up the standard output streams for each thread. Slaves will by default
# have the null device opened in both cases to suppress output
if not master:
   sys.stdout = open(os.devnull, 'w')

# Check that we have AMBERHOME set
if not os.getenv('AMBERHOME'):
   raise MMPBSA_Error('AMBERHOME is not set')

# Give each thread its own timer. For some reason, it seems like each thread
# will work on the same piece of memory. This seems to me to be an mpi4py bug,
# but we'd need an array for a Gather at the end, anyway.

timers = [Timer() for i in range(size)]
MMPBSA_Timer = timers[rank]

# Start the setup timer
MMPBSA_Timer.AddTimer('setup', 'Total setup time:')
MMPBSA_Timer.StartTimer('setup')

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#
# BEGIN set up command-line parser options and and parse the CL arguments
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

# Set up the OptionParser, adding the options and 
parser = OptionParser(sys.stdout, sys.stderr)
parser.addOption('-O', 'overwrite', help='Overwrite existing output files',
                 default=False, num_entries=0)
parser.addOption('-i', 'input_file', help='MM/PBSA input file',
                 default=None)
parser.addOption('-o', 'output_file', default='FINAL_RESULTS_MMPBSA.dat', 
                 help='Final MM/PBSA statistics file. Default ' +
                 'FINAL_RESULTS_MMPBSA.dat',)
parser.addOption('-sp', 'solvated_prmtop',
                 help='Solvated complex topology file')
parser.addOption('-cp', 'complex_prmtop', default='complex_prmtop',
                  help='Complex topology file. Default "complex_prmtop"')
parser.addOption('-rp', 'receptor_prmtop', help='Receptor topology file')
parser.addOption('-lp', 'ligand_prmtop', help='Ligand topology file')
parser.addOption('-y', 'mdcrd', num_entries=-1, default=['mdcrd'], 
                 help='Input trajectories to analyze. Default mdcrd')
parser.addOption('-do', 'decompout', default='FINAL_DECOMP_MMPBSA.dat',
                 help='Decomposition statistics summary file. Default ' +
                 'FINAL_DECOMP_MMPBSA.dat')
parser.addOption('-eo', 'energyout', default=None,
                 help='CSV-format output of all energy terms for every frame ' +
                 'in every calculation. File name forced to end in .csv')
parser.addOption('-deo', 'dec_energies', default=None,
                 help='CSV-format output of all decomposition energy terms ' +
                 'for every frame. File name forced to end in .csv')
parser.addOption('-yr', 'receptor_mdcrd', help='Receptor trajectory file for ' +
                 'multiple trajectory approach', num_entries=-1)
parser.addOption('-yl', 'ligand_mdcrd', help='Ligand trajectory file for ' +
                 'multiple trajectory approach', num_entries=-1)
parser.addOption('-mc', 'mutant_complex_prmtop',
                 help='Alanine scanning mutant complex topology file',
                 default='mutant_complex_prmtop')
parser.addOption('-ml', 'mutant_ligand_prmtop',
                 help='Alanine scanning mutant ligand topology file')
parser.addOption('-mr', 'mutant_receptor_prmtop',
                 help='Alanine scanning mutant receptor topology file')
parser.addOption('-slp', 'solvated_ligand_prmtop',
                 help='Solvated ligand topology file')
parser.addOption('-srp', 'solvated_receptor_prmtop',
                 help='Solvated receptor topology file')
parser.addOption('-xvvfile', 'xvvfile', help='XVV file for 3D-RISM. Default ' +
                 '$AMBERHOME/dat/mmpbsa/spc.xvv', default=os.path.join(
                 os.getenv('AMBERHOME'), 'dat', 'mmpbsa', 'spc.xvv'))
parser.addOption('-make-mdins', 'make_mdins', num_entries=0, default=False,
                 help='Create the Input files for each calculation and quit')
parser.addOption('-use-mdins', 'use_mdins', num_entries=0, default=False,
                 help='Use existing input files for each calculation')
parser.addOption('-rewrite-output', 'rewrite_output', num_entries=0, 
                 default=False, help="Don't rerun any calculations, just parse "
                 + 'existing output files')
parser.addOption('--clean', 'clean', num_entries=0, default=False,
                 help='Clean temporary files from previous run')
parser.SetHelp(['--help', '-h', '--h', '-H'])

FILES = parser.Parse()

if not FILES: sys.exit(0)

if FILES.energyout and not FILES.energyout.endswith('.csv'):
   FILES.energyout += '.csv'

if size > 1: print 'Running MMPBSA.MPI on %s processors' % size

print 'Reading command-line arguments and input files...'

if FILES.clean:
   print 'Cleaning temporary files and quitting'
   utils.remove(0)
   sys.exit(0)

if not FILES.rewrite_output and not FILES.input_file:
   parser.print_help(sys.stderr)
   sys.exit(1)

if not FILES.receptor_prmtop and not FILES.ligand_prmtop: 
   stability = True
else: 
   stability = False

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#
# BEGIN set up input file and parse it. Only master does this, then broadcasts
#       it to the rest of the threads. Then some setup stuff is done immediately
#       afterwards with the input values
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

# We have to declare INPUT as an empty dictionary for the slave nodes bcast call
INPUT = {}

# Now if the user set -rewrite-output, read the _MMPBSA_info file on the master,
# read it, then broadcast it out to everybody else. Make sure it stays in the
# format [ <var>  =  <value> ] to limit what it can do. This can still be
# hijacked since we're using "exec" to load these variables, but it's quite a
# bit safer.
if FILES.rewrite_output:
   if size > 1:
      if master:
         print >> sys.stderr, 'Error: -rewrite-output must be used in serial!'
      sys.exit(1)
   valid_line = re.compile(r'.+ = .+')
   try: infofile = open('_MMPBSA_info', 'r')
   except IOError: 
      raise MMPBSA_Error('_MMPBSA_info is needed for -rewrite-output!')
   restore_info = infofile.readlines()
   infofile.close()
   for line in restore_info:
      # Skip over illegal lines
      if not valid_line.match(line): continue
      if line.startswith('|'): continue # This is the input file text
      else: 
         try: exec(line.strip())
         except: raise MMPBSA_Error('Corrupt _MMPBSA_info file!')

# Now read the input. If we're rewriting output, we still do this to allow
# users to change some settings (like verbose settings, etc.)
if master:
   input_file = InputFile()

   strip_mask = ':WAT,Cl-,CIO,Cs+,IB,K+,Li+,MG2,Na+,Rb+'

# Add namelists with a list of variables. The variables are added to the 
# namelists in lists. The entries are:
# [<variable name> <variable type> <default value> <# of characters to match>]

   input_file.addNamelist('general', 'general', 
                          [ 
                            ['debug_printlevel', 'int', 0, 4],
                            ['endframe', 'int', 9999999, 4],
                            ['entropy', 'int', 0, 4],
                            ['full_traj', 'int', 0, 4],
                            ['interval', 'int', 1, 4],
                            ['keep_files', 'int', 1, 4],
                            ['ligand_mask', 'str', '', 4],
                            ['netcdf', 'int', 0, 4],
                            ['receptor_mask', 'str', '', 4],
                            ['search_path', 'int', 0, 4],
                            ['startframe', 'int', 1, 4],
                            ['strip_mask', 'str', strip_mask, 4],
                            ['use_sander', 'int', 0, 4],
                            ['verbose', 'int', 1, 4]
                          ], trigger = None)

   input_file.addNamelist('gb', 'gb', 
                          [ 
                            ['ifqnt', 'int', 0, 4],
                            ['igb', 'int', 5, 4],
                            ['qm_theory', 'str', '', 4],
                            ['qm_residues', 'str', '', 4],
                            ['qmcharge_com', 'float', 0, 10],
                            ['qmcharge_lig', 'float', 0, 10],
                            ['qmcharge_rec', 'float', 0, 10],
                            ['qmcut', 'float', 9999, 4],
                            ['saltcon', 'float', 0, 4],
                            ['surfoff', 'float', -999999.0, 5],
                            ['surften', 'float', 0.0072, 5]
                          ], trigger = 'gbrun' )

   input_file.addNamelist('pb', 'pb', 
                          [ 
                            ['cavity_offset', 'float', -0.5692, 8],
                            ['cavity_surften', 'float', 0.0378, 8],
                            ['exdi', 'float', 80, 4],
                            ['fillratio', 'float', 4, 4],
                            ['indi', 'float', 1, 4],
                            ['inp', 'int', 2, 3],
                            ['istrng', 'float', 0, 4],
                            ['linit', 'int', 1000, 4],
                            ['prbrad', 'float', 1.4, 4],
                            ['radiopt', 'int', 1, 4],
                            ['sander_apbs', 'int', 0, 4],
                            ['scale', 'float', 2, 4]
                          ], trigger='pbrun')

   input_file.addNamelist('ala', 'alanine_scanning', 
                          [ ['mutant_only', 'int', 0, 4]], 
                          trigger='alarun')

   input_file.addNamelist('nmode', 'nmode', 
                          [ 
                            ['dielc', 'float', 4, 4],
                            ['drms', 'float', 0.001, 4],
                            ['maxcyc', 'int', 10000, 4],
                            ['nminterval', 'int', 1, 4],
                            ['nmendframe', 'int', 1000000, 4],
                            ['nmode_igb', 'int', 1, 8],
                            ['nmode_istrng', 'float', 0, 8],
                            ['nmstartframe', 'int', 1, 4]
                          ], trigger='nmoderun')

   input_file.addNamelist('decomp', 'decomposition', 
                          [ 
                            ['csv_format', 'int', 1, 4],
                            ['dec_verbose', 'int', 0, 4],
                            ['idecomp', 'int', 0, 4],
                            ['print_res', 'str', 'all', 4]
                          ], trigger='decomprun')
   
   input_file.addNamelist('rism', 'rism', 
                          [ 
                            ['buffer', 'float', 14, 4],
                            ['closure', 'str', 'kh', 7],
                            ['closureorder', 'int', 1, 8],
                            ['grdspc', 'float', 0.5, 4],
                            ['ng', 'str', '-1,-1,-1', 2],
                            ['polardecomp', 'int', 0, 4],
                            ['rism_verbose', 'int', 0, 4],
                            ['solvbox', 'str', '-1,-1,-1', 5],
                            ['solvcut', 'float', None, 5],
                            ['thermo', 'str', 'std', 4],
                            ['tolerance', 'float', 1.0e-5, 4]
                          ], trigger='rismrun')

   if FILES.input_file: INPUT = input_file.Parse(FILES.input_file)

# Let everyone know what the INPUT is
INPUT = MPI.COMM_WORLD.bcast(INPUT, root=0)

# Set the debug_printlevel. To trace errors that happen before here, the
# default value will have to be changed at the top of this script, but
# that shouldn't be necessary unless you modify commandline_parser.py or
# input_parser.py, or an unexpected bug in one of those pops up
debug_printlevel = INPUT['debug_printlevel']

# All of this stuff has already been done if we are rewriting output file
if not FILES.rewrite_output:
   # Invert scale
   INPUT['scale'] = 1 / INPUT['scale']

   # Set up netcdf variables and decide trajectory suffix
   if INPUT['netcdf'] == 0:
      INPUT['netcdf'] = ''
      trjsuffix = 'mdcrd'
   else:
      INPUT['netcdf'] = 'netcdf'
      trjsuffix = 'nc'

   # Set default GBSA for Decomp
   if INPUT['decomprun']:
      INPUT['gbsa'] = 2

   # If we're running a stability calculation, only verbose == 2 makes sense,
   # since no terms cancel and there's no DELTA term...
   if stability: INPUT['verbose'] = 2

   # Make thermo option case-insensitive
   INPUT['thermo'] = INPUT['thermo'].lower()

   # RISM stuff
   if not INPUT['solvcut']: INPUT['solvcut'] = INPUT['buffer']
   # We want to trigger different kinds of RISM free energy printouts based on
   # what kinds of thermo output we want
   INPUT['rismrun_std'] = INPUT['rismrun'] and INPUT['thermo'] in ['std','both']
   INPUT['rismrun_gf'] = INPUT['rismrun'] and INPUT['thermo'] in ['gf', 'both']

   # Add the temperature to the INPUT dictionary. Don't add it to any namelist 
   # above, though, since you really have to modify the nmode and ptraj source
   # codes to use a different temperature. Thus, I don't want to make it too
   # easy to get users "thinking" they changed the temperature correctly.
   INPUT['temp'] = 298.15

# Check for incompatibilities
if master: utils.CheckIncomp(INPUT, stability)

MPI.COMM_WORLD.Barrier()

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#
# BEGIN Set up prmtop classes and get the necessary masks
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

print 'Loading and checking parameter files for compatibility...'

normal_system = MMPBSA_System(FILES.complex_prmtop, FILES.receptor_prmtop,
                              FILES.ligand_prmtop)

# If we have a chamber prmtop, make sure we force use_sander and opt out of any
# calculation that requires a NAB program which doesn't have CHAMBER potential
# terms implemented...
if normal_system.complex_prmtop.chamber:
   INPUT['use_sander'] = True
   if INPUT['rismrun']: 
      raise MMPBSA_Error('CHAMBER prmtops cannot be used with 3D-RISM')
   if INPUT['nmoderun']:
      raise MMPBSA_Error('CHAMBER prmtops cannot be used with 3D-RISM')
   print 'CHAMBER prmtops found. Forcing use of sander'

normal_system.Map()
normal_system.CheckConsistency()

# Print warnings if we're overwriting any masks
if not INPUT['ligand_mask'] and INPUT['receptor_mask']:
   sys.stderr.write('receptor_mask overwritten with default\n')

if not INPUT['receptor_mask'] and INPUT['ligand_mask']:
   sys.stderr.write('ligand_mask overwritten with default\n')

# Get default masks
if not INPUT['ligand_mask'] or not INPUT['receptor_mask']:
   com_mask, INPUT['receptor_mask'], INPUT['ligand_mask'] = normal_system.Mask(
      'all', in_complex=True)

if INPUT['alarun']:
   # No need to map this system, as the normal one was mapped, and the 
   # mapping should be the same.
   if not FILES.mutant_receptor_prmtop and not FILES.mutant_ligand_prmtop:
      raise MMPBSA_Error('Alanine scanning requires either a mutated ' +
                         'receptor or mutated ligand topology file!')

   if not FILES.mutant_receptor_prmtop:
      FILES.mutant_receptor_prmtop = FILES.receptor_prmtop

   if not FILES.mutant_ligand_prmtop:
      FILES.mutant_ligand_prmtop = FILES.ligand_prmtop

   mutant_system = MMPBSA_System(FILES.mutant_complex_prmtop,
                                 FILES.mutant_receptor_prmtop,
                                 FILES.mutant_ligand_prmtop)
   mutant_system.Map()
   mutant_system.CheckConsistency()

else: mutant_system = None

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#
# BEGIN Check that we're doing a valid calculation, and remove old files
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

# Check to make sure we're doing a valid calculation
if not (INPUT['gbrun'] or INPUT['pbrun'] or INPUT['rismrun'] or 
        INPUT['nmoderun'] or INPUT['entropy']):
   raise MMPBSA_Error('No calculation type specified!')

# Now we're getting ready, remove existing files, but leave mdin files intact
# if we're using existing mdins
if master and FILES.use_mdins and not FILES.rewrite_output:
   utils.remove(-1)
elif master and not FILES.rewrite_output:
   utils.remove(0)

# Find the executables

# Declare external_progs for bcast
external_progs = {}

# Only need to find external programs if we're doing a real calculation
if not FILES.make_mdins and not FILES.rewrite_output:
   # Find our programs and broadcast them
   if master: external_progs = find_progs(INPUT)
   external_progs = MPI.COMM_WORLD.bcast(external_progs, root=0)

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#
# BEGIN Create input files based on INPUT dictionary
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

if master and not FILES.rewrite_output and not FILES.use_mdins:
   create_inputs(INPUT, normal_system)

MMPBSA_Timer.StopTimer('setup')

# If we only wanted to generate our MDIN files, bail out now
if FILES.make_mdins:
   print 'Created mdin files. Quitting.'
   sys.exit(0)

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#
# BEGIN Create trajectory files necessary for execution
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

if not FILES.rewrite_output:
   MMPBSA_Timer.AddTimer('cpptraj', 'Creating trajectories with cpptraj:')
   MMPBSA_Timer.StartTimer('cpptraj')

   if master:
      print 'Preparing trajectories for simulation...'
      numframes, numframes_nmode = make_trajectories(INPUT, FILES, size, 
             str(external_progs['cpptraj']))
   
   MPI.COMM_WORLD.Barrier()
   
   MMPBSA_Timer.StopTimer('cpptraj')
   
   MMPBSA_Timer.AddTimer('muttraj', 'Mutating trajectories:')
   MMPBSA_Timer.StartTimer('muttraj')
   
   if INPUT['alarun']: print 'Mutating trajectories...'
   mut_str, mutant_residue = make_mutant_trajectories(INPUT, FILES, rank, 
                   str(external_progs['cpptraj']), normal_system, mutant_system)

   MPI.COMM_WORLD.Barrier()
   
   if master:
      print '%s frames were processed by cpptraj ' % numframes + \
            'for use in calculation.'
      if INPUT['nmoderun']: 
         print '%s frames were processed by cpptraj ' % numframes_nmode + \
               'for nmode calculations.'
   
   MMPBSA_Timer.StopTimer('muttraj')

   # Add all of the calculation timers
   MMPBSA_Timer.AddTimer('calc', 'Total calculation time:')
   if INPUT['gbrun']: MMPBSA_Timer.AddTimer('gb', 'Total GB calculation time:')
   if INPUT['pbrun']: MMPBSA_Timer.AddTimer('pb', 'Total PB calculation time:')
   if INPUT['rismrun']: 
      MMPBSA_Timer.AddTimer('rism', 'Total 3D-RISM calculation time:')
   if INPUT['nmoderun']: 
      MMPBSA_Timer.AddTimer('nmode', 'Total normal mode calculation time:')
   if INPUT['entropy']:
      MMPBSA_Timer.AddTimer('qh', 'Total quasi-harmonic calculation time:')

   #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
   #
   # BEGIN Running calculations
   #
   #-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

   MMPBSA_Timer.StartTimer('calc')

   print ''

   if not INPUT['mutant_only']:
      if master and INPUT['alarun']: 
         print 'Running calculations on normal system...\n'
      run_calculations(FILES, INPUT, rank, MMPBSA_Timer, 
                       external_progs, '_MMPBSA_', normal_system)
   if INPUT['alarun']:
      if master:
         print 'Running calculations on alanine mutant...\n'
      run_calculations(FILES, INPUT, rank, MMPBSA_Timer,
                       external_progs, '_MMPBSA_mutant_', mutant_system)

   MMPBSA_Timer.StopTimer('calc')

# We are now done with all calculations, let everybody catch up
MPI.COMM_WORLD.Barrier()

MMPBSA_Timer.AddTimer('output','Statistics calculation & output writing:')
MMPBSA_Timer.StartTimer('output')

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#
# BEGIN writing _MMPBSA_info file for re-writing simulations
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

# Now write the _MMPBSA_info file. This needs to have variable declarations
# for every variable that is necessary to write all of the output files.
# Therefore we back up INPUT, FILES, and a couple other variables
if master and not FILES.rewrite_output:
   infofile = open('_MMPBSA_info', 'w')
   # First back up INPUT
   for key in INPUT.keys():
      if type(INPUT[key]).__name__ == 'str':
         infofile.write("INPUT['%s'] = '%s'" % (key, INPUT[key]) + os.linesep)
      else:
         infofile.write("INPUT['%s'] = %s" % (key, INPUT[key]) + os.linesep)
   # Now back up FILES
   for item in dir(FILES):
      # Ignore the built-in junk
      if item.startswith('__'): continue
      # We also don't want to store input/output file names so they can be reset
      if item in ['input_file', 'output_file', 'decompout', 'energyout',
                  'dec_energies', 'rewrite_output', 'overwrite']: continue
      if type(getattr(FILES,item)).__name__ == 'str':
         infofile.write('FILES.%s = "%s"'%(item,getattr(FILES,item))+os.linesep)
      else:
         infofile.write('FILES.%s = %s'%(item,getattr(FILES,item))+os.linesep)
   # Now we need the MPI size
   infofile.write('size = %d' % size + os.linesep)
   # Now we need the number of frames
   infofile.write('numframes = %d' % numframes + os.linesep)
   infofile.write('numframes_nmode = %d' % numframes_nmode + os.linesep)
   # Now we need the mutant string
   if mut_str: infofile.write('mut_str = "%s"' % mut_str + os.linesep)
   else: infofile.write('mut_str = None' + os.linesep)
   infofile.write('stability = %s' % stability + os.linesep)
   # Now print the original mmpbsa_input file 
   infofile.write(str(input_file))
   infofile.close()

#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#
# BEGIN writing final output files and timing info
#
#-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

# Only the master does output file creation
if master:
   # We want to display the text of the new input file if we used rewrite-output
   # but still supplied a new input file to change some verbose settings, so 
   # load that text here and send it over.
   if FILES.rewrite_output and FILES.input_file:
      infile_string = str(input_file)
   else: infile_string = ''

   if stability: 
      write_stability_output(FILES, INPUT, size, normal_system, infile_string,
                             [numframes, numframes_nmode], mut_str)
   else: 
      write_binding_output(FILES, INPUT, size, normal_system, infile_string,
                           [numframes, numframes_nmode], mut_str)

   if INPUT['decomprun']:
      if stability: 
         write_decomp_stability_output(FILES, INPUT, size, normal_system,
                                       mutant_system, mut_str)
      else: 
         write_decomp_binding_output(FILES, INPUT, size, normal_system, 
                                     mutant_system, mut_str)

MMPBSA_Timer.StopTimer('output')

MMPBSA_Timer.Done()

if not master: sys.exit(0)

# Now print out timing:

if size > 1: print '\nTiming (master node):'
else: print '\nTiming:'

print MMPBSA_Timer.Print('setup')

if not FILES.rewrite_output:
   print MMPBSA_Timer.Print('cpptraj')

   if INPUT['alarun']:
      print MMPBSA_Timer.Print('muttraj')

   print MMPBSA_Timer.Print('calc')
   print ''

   if INPUT['gbrun']:
      print MMPBSA_Timer.Print('gb')

   if INPUT['pbrun']:
      print MMPBSA_Timer.Print('pb')

   if INPUT['nmoderun']:
      print MMPBSA_Timer.Print('nmode')

   if INPUT['entropy']:
      print MMPBSA_Timer.Print('qh')

   print ''

print MMPBSA_Timer.Print('output')
print MMPBSA_Timer.Print('global')

utils.remove(INPUT['keep_files'])

print "\n\nMMPBSA.py Finished! Thank you for using. Please report any " + \
      "bugs to amber@ambermd.org"

sys.exit(0) # quit without error
