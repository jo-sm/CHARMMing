#! /usr/bin/python2.6

# Load system modules. Do NOT load os or sys, since we run "exec", and that can
# potentially allow people to do nasty things if given access to os and sys
import sys, math
from os.path import exists
from optparse import OptionParser

# Load custom modules
from ParmedTools.logos import Logo
from ParmedTools.exceptions import ParmError
from ParmedTools.utils import LineToCmd, PrintHelp, PrintHelpSummary
from chemistry.amber.readparm import AmberParm
from ParmedTools import ParmedActions
from ParmedTools import __version__

# Set up new excepthook to clean up fatal exception printouts
debug = False
prompt = '' # No prompt by default, unless we are reading from sys.stdin
stripped_once = False # We can't strip more than once

def excepthook(exception_type, exception_value, tb):
   """ Replaces excepthook so we can control the printing of 
       tracebacks. Those are helpful for debugging purposes, but 
       may be unsightly to users. debug set above controls this
       behavior, and a CL flag allows users to set this.
   """
   import traceback
   if debug: traceback.print_tb(tb)
   sys.stderr.write('%s: %s\n' % (exception_type.__name__, exception_value))

sys.excepthook = excepthook

# Set up parser
parser = OptionParser("Usage: %prog [Options] <prmtop> <input_file>",
                      version="%%prog: Version %s" % __version__)
parser.add_option('-d', '--debug', dest='debug', 
                  action='store_true', default=False,
                  help='Print how each action is parsed and show which line ' +
                  'of which file an error occurred on (verbose tracebacks). ' +
                  'OFF by default.')
parser.add_option('-e', '--enable-interpreter', dest='interpreter',
                  action='store_true', default=False,
                  help='Allow the use of ! and !! to permit the user to drop ' +
                  'into a limited Python interpreter with access to the ' +
                  'topology file API. This is an advanced option, use with ' +
                  'care.')
parser.add_option('-q', '--quiet', dest='debug',
                  action='store_false', 
                  help='Disable verbose tracebacks. Reverses -d/--debug')
parser.add_option('-p', '--prompt', dest='prompt', default='>',
                  help='Which character/string to use as a command prompt')
parser.add_option('-n', '--no-splash', dest='printlogo', action='store_false',
                  help='Prevent printing the greeting logo', default=True)

(opt, args) = parser.parse_args()

debug = opt.debug

# Load the splash screen
if opt.printlogo:
   splash = Logo()
   print splash

if len(args) > 2:
   print 'Unexpected command-line arguments!'
   parser.print_usage()
   sys.exit(1)

if len(args) > 0:
   amber_prmtop = AmberParm(args[0])
   if not amber_prmtop.valid:
      raise ParmError('Invalid prmtop %s' % args[0])
   print 'Loaded Amber topology file %s\n' % args[0]

else:
   sys.stdout.write('Choose an Amber Topology file to load: ')
   prm_name = sys.stdin.readline().strip()
   amber_prmtop = AmberParm(prm_name)
   if not amber_prmtop.valid:
      raise ParmError('Invalid prmtop %s' % prm_name)
   print 'Loaded Amber topology file %s\n' % prm_name

if len(args) > 1:
   if not exists(args[1]):
      raise ParmError('Missing file %s' % args[1])
   print 'Reading actions from %s\n' % args[1]
   input = open(args[1], 'r')

else:
   sys.stdout.write('Reading input from STDIN...\n')
   input = sys.stdin
   prompt = opt.prompt.strip() + ' '

output_parm = None

sys.stdout.write(prompt)
line = input.readline()

while line and line.strip().lower() != 'go':
   # Skip over comments
   if line.strip().startswith('#') or len(line.strip()) == 0: 
      sys.stdout.write(prompt)
      line = input.readline()
      continue
   # Allow for in-line comments
   if '#' in line: line = line[:line.index('#')].strip()

   # We allow any line that starts with ! to represent additional commands not
   # present in ParmedActions. This allows for more direct access to the 
   # AmberParm class, and effectively opens up a limited python interpreter
   # inside parmed.py. Multiline, formatted code is read after a !! line, and
   # stops reading after the next !! line
   if line.strip() == '!!':
      if not opt.interpreter:
         raise ParmError('ParmEd not in interpreter mode! Use the -e flag ' +
                         'to enable this feature. Use with care.')
      if prompt: sys.stdout.write('py >>> ')
      line = input.readline()
      code = ''
      while line and line.strip() != '!!':
         if 'import' in line.split():
            raise ParmError('import statements not allowed via ! for ' +
                            'security reasons')
         code += line
         if prompt: sys.stdout.write('py >>> ')
         line = input.readline()
      exec(code)
      # If an EOF appeared in the Python section, bail out here
      if not line: break
      sys.stdout.write(prompt)
      line = input.readline()
      continue
   elif line.strip().startswith('!'):
      if not opt.interpreter:
         raise ParmError('ParmEd not in interpreter mode! Use the -e flag ' +
                         'to enable this feature. Use with care.')
      line = line[1:].strip()
      if 'import' in line.split():
         raise ParmError('import statements not allowed via ! for ' +
                         'security reasons')
      exec(line)
      sys.stdout.write(prompt)
      line = input.readline()
      continue

   # Tokenize the string and turn keywords into arguments
   evalstring = LineToCmd(line)

   if debug and not line.strip().startswith('help'):
      print '  Line parsed as [%s]' % evalstring

   # Tack parmout onto the end
   if line.strip().startswith('parmout'):
      output_parm = evalstring

   # Answer calls for help
   elif line.strip().startswith('help'):
      help_parts = line.split()
      try: 
         usage = ParmedActions.usages[help_parts[1].lower()]
         PrintHelp(usage, getattr(ParmedActions, help_parts[1].lower()).__doc__)
      except KeyError:
         print 'No action %s\n' % help_parts[1]
      except IndexError:
         PrintHelpSummary(ParmedActions.usages)

   # Queue the rest
   else:
      # Make sure we only strip masks once
      if line.split()[0].lower() == 'strip':
         if stripped_once: 
            raise ParmError('Can only strip atoms from a topology file once!')
         else: stripped_once = True
      # Catch any non-existent commands (don't make this quit in error)
      if not hasattr(ParmedActions, line.split()[0].lower()):
         print '  Unrecognized action %s' % (line.strip().split()[0])
      else:
         action = eval(evalstring)
         print '  %s' % action
         action.execute()

   sys.stdout.write(prompt)
   line = input.readline()

# end while input

# Evaluate the parmout command if there was one
if output_parm: 
   parmout_action = eval(output_parm)
   parmout_action.execute()

print 'Done!'
