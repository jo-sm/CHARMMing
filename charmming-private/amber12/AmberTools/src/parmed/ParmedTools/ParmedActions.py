""" 
All of the prmtop actions used in PARMED. Each class is a separate action. You
MUST assume that each argument comes into the action's __init__ as a string, except
parm, which is passed as an amberParm object!
"""

import sys
from ParmedTools.exceptions import (WriteOFFError, ParmError, ParmedMoleculeError,
                                    ChangeStateError, CoarseGrainError,
                                    ChangeLJPairError, ParmedChangeError,
                                    SetParamError, DeleteDihedralError)
from chemistry.amber.mask import AmberMask
from chemistry.amber.readparm import (_Bond, _BondType, _Angle, _AngleType, 
                                      _Dihedral, _DihedralType)
import math

# Add a help dictionary entry for each additional Action added to this class!
# Each help entry should be a list with 2 elements: [Usage, description]
usages = {
              'help' : 'help [<action>]',
           'parmout' : 'parmout <prmtop_name>',
      'setoverwrite' : 'setOverwrite [True|False]', 
       'writefrcmod' : 'writeFrcmod <frcmod_name>',
        'loadrestrt' : 'loadRestrt <restrt_filename>',
          'writeoff' : 'writeOFF <OFF_filename>',
       'changeradii' : 'changeRadii <radii_set>',
      'changeljpair' : 'changeLJPair <mask1> <mask2> <Rmin> <epsilon>',
    'changelj14pair' : 'changeLJ14Pair <mask1> <mask2> <Rmin> <epsilon>',
     'checkvalidity' : 'checkValidity',
            'change' : 'change <property> <mask> <new_value>',
         'printinfo' : 'printInfo <flag>',
         'addljtype' : 'addLJType <mask> [<new_radius> <new_epsilon>]',
           'outparm' : 'outparm <prmtop_name>',
      'printljtypes' : 'printLJTypes [<mask>]',
              'scee' : 'scee <scee_value>',
              'scnb' : 'scnb <scnb_value>',
'changeljsingletype' : 'changeLJSingleType <mask> <radius> <depth>',
      'printdetails' : 'printDetails <mask>',
        'printflags' : 'printFlags',
     'printpointers' : 'printPointers',
        'printbonds' : 'printBonds <mask>',
       'printangles' : 'printAngles <mask>',
    'printdihedrals' : 'printDihedrals <mask>',
      'setmolecules' : 'setMolecules [solute_ions=True|False]',
  'combinemolecules' : 'combineMolecules <mol_id1> [<mol_id2>]',
                'go' : 'go',
              'quit' : 'quit',
#   'addcoarsegrain' : 'addCoarseGrain <parameter_file>',
   'changeprotstate' : 'changeProtState <mask> <state #>',
         'netcharge' : 'netCharge [<mask>]',
             'strip' : 'strip <mask>',
     'definesolvent' : 'defineSolvent <residue list>',
     'addexclusions' : 'addExclusions <mask1> <mask2>',
       'adddihedral' : 'addDihedral <mask1> <mask2> <mask3> <mask4> <phi_k> ' +
                                   '<per> <phase> <scee> <scnb> [<type>]',
           'setbond' : 'setBond <mask1> <mask2> <k> <Req>',
          'setangle' : 'setAngle <mask1> <mask2> <mask3> <k> <THETeq>',
   'addatomicnumber' : 'addAtomicNumber',
    'deletedihedral' : 'deleteDihedral <mask1> <mask2> <mask3> <mask4>'
}

# Add help and go as a class here to basically hold docstrings for the help
# function

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class help(object):
   """ Prints usage message and short description for <action>, or prints all
       allowable options with their usage statement (but no description)
   """

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class go(object):
   'Stops waiting for commands and executes the last parmout statement'

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Action(object):
   """ Base action class """
   def __init__(self, parm):
      """ Defines the prmtop """
      self.parm = parm

   def execute(self):
      """ Commands involved in executing the action """
      pass

   def __str__(self):
      return ''

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class parmout(Action):
   """ Final prmtop written after all actions are complete """
   def __init__(self, parm, filename, rst_name=None):
      Action.__init__(self, parm)
      self.filename = filename
      self.rst_name = rst_name

   def __str__(self):
      return 'Outputting Amber topology file %s' % self.filename

   def execute(self):
      self.parm.writeParm(self.filename)
      if self.rst_name:
         try: self.parm.writeRst7(self.rst_name)
         except AttributeError:
            raise ParmError('You must load coordinates to write a restart')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class setoverwrite(Action):
   """ Necessary to overwrite original topology file name. """
   def __init__(self, parm, owrite = True):
      Action.__init__(self, parm)
      self.overwrite = bool(eval(owrite))

   def __str__(self):
      if self.overwrite:
         return 'Prmtop is overwritable'
      else:
         return 'Prmtop is NOT overwritable'
   
   def execute(self):
      self.parm.overwrite = self.overwrite

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class writefrcmod(Action):
   " Writes an frcmod file from all of the parameters in the topology file. "
   def __init__(self, parm, frcmod='frcmod'):
      Action.__init__(self, parm)
      self.frcmod_name = frcmod

   def __str__(self):
      return 'Dumping FRCMOD file %s' % self.frcmod_name

   def execute(self):
      self.parm.frcmod(self.frcmod_name)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class loadrestrt(Action):
   """Loads a restart file so we have coordinates. Necessary for distance-based
      mask criteria and writeOFF
   """
   def __init__(self, parm, rst7):
      Action.__init__(self, parm)
      self.rst_name = rst7
      self.parm.LoadRst7(rst7)

   def __str__(self):
      return 'Loading restart file %s' % self.rst_name

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class writeoff(Action):
   "Writes an Amber OFF Library with all of the residues found in the topology"
   def __init__(self, parm, filename):
      Action.__init__(self, parm)
      self.off_file = filename

   def __str__(self):
      return 'Writing Amber OFF file %s' % self.off_file

   def execute(self):
      try: self.parm.rst7
      except: raise WriteOFFError('You must load a restart for WriteOFF!')

      self.parm.writeOFF(self.off_file)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changeradii(Action):
   """ Changes intrinsic GB radii to the specified set: bondi, mbondi, amber6,
       mbondi2, or mbondi3
   """
   def __init__(self, parm, radii_set='mbondi'):
      Action.__init__(self, parm)
      self.radii = radii_set

   def __str__(self):
      return 'Changing PB/GB radii to %s' % self.radii

   def execute(self):
      from ParmedTools.changeradii import ChRad
      # Add RADIUS_SET to prmtop if it's not there already, and a blank 
      # description, since it's about to be set here
      if not 'RADIUS_SET' in self.parm.flag_list:
         self.parm.addFlag('RADIUS_SET', '1a80', num_items=1)
      # Add RADII prmtop section if it doesn't exist already. Just make it a
      # zeroed array, since it's all about to be set here
      if not 'RADII' in self.parm.flag_list:
         self.parm.addFlag('RADII', '5E16.8', num_items=self.parm.ptr('natom'))
      ChRad(self.parm, self.radii)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changeljpair(Action):
   """ Changes a particular Lennard Jones pair based on a given (pre-combined)
       epsilon/Rmin 
   """
   def __init__(self, parm, atom_type1, atom_type2, rmin, eps):
      Action.__init__(self, parm)
      self.mask1 = AmberMask(self.parm, atom_type1)
      self.mask2 = AmberMask(self.parm, atom_type2)
      self.rmin = float(rmin)
      self.eps = float(eps)

   def __str__(self):
      return ('Setting LJ %s-%s pairwise interaction to have ' +
              'Rmin = %16.5f and Epsilon = %16.5f') % (self.mask1, self.mask2,
              self.rmin, self.eps)

   def execute(self):
      from ParmedTools.changeljpair import ChLJPair
      selection1 = self.mask1.Selection()
      selection2 = self.mask2.Selection()
      if sum(selection1) == 0 or sum(selection2) == 0:
         print >> sys.stderr, 'Skipping empty masks in changeLJPair'
         return 0
      # Make sure we only selected 1 atom type in each mask
      attype1 = None
      attype2 = None
      for i in range(self.parm.ptr('natom')):
         if selection1[i] == 1:
            if not attype1:
               attype1 = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
            else:
               if attype1 != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                  raise ChangeLJPairError(
                                      'First mask matches multiple atom types!')
         if selection2[i] == 1:
            if not attype2:
               attype2 = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
            else:
               if attype2 != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                  raise ChangeLJPairError(
                                     'Second mask matches multiple atom types!')
      ChLJPair(self.parm, attype1, attype2, self.rmin, self.eps)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changelj14pair(Action):
   """ Changes a particular 1-4 Lennard Jones pair based on a given
       (pre-combined) epsilon/Rmin. Only valid for CHAMBER prmtops
   """
   def __init__(self, parm, atom_type1, atom_type2, rmin, eps):
      Action.__init__(self, parm)
      # Make sure this is a chamber prmtop
      if not self.parm.chamber:
         raise ChangeLJPairError('Changing 1-4 NB pairs makes no sense for ' +
                                 'non-chamber created prmtops!')
      # If not, initiate instance data
      self.mask1 = AmberMask(self.parm, atom_type1)
      self.mask2 = AmberMask(self.parm, atom_type2)
      self.rmin = float(rmin)
      self.eps = float(eps)

   def __str__(self):
      if self.parm.chamber:
         return ('Setting LJ 1-4 %s-%s pairwise interaction to have ' +
              '1-4 Rmin = %16.5f and 1-4 Epsilon = %16.5f') % (self.mask1,
              self.mask2, self.rmin, self.eps)

   def execute(self):
      from ParmedTools.changeljpair import ChLJPair
      selection1 = self.mask1.Selection()
      selection2 = self.mask2.Selection()
      if sum(selection1) == 0 or sum(selection2) == 0:
         print >> sys.stderr, 'Skipping empty masks in changeLJ14Pair'
         return None
      # Make sure we only selected 1 atom type, and figure out what it is
      attype1 = None
      attype2 = None
      for i in range(self.parm.ptr('natom')):
         if selection1[i] == 1:
            if not attype1:
               attype1 = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
            else:
               if attype1 != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                  raise ChangeLJPairError(
                                      'First mask matches multiple atom types!')
         if selection2[i] == 1:
            if not attype2:
               attype2 = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
            else:
               if attype2 != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                  raise ChangeLJPairError(
                                     'Second mask matches multiple atom types!')

      # Adjust 1-4 non-bonded terms as well if we're using a chamber-prmtop
      ChLJPair(self.parm, attype1, attype2, self.rmin, self.eps, True)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class checkvalidity(Action):
   """Basic checks for prmtop validity."""
   def __str__(self):
      return 'Determining validity of prmtop'

   def execute(self):
      valid = True
      if not self.parm.valid:
         print >> sys.stderr, 'AmberParm detects invalid prmtop file'
         valid = False

      # Check that ATOMS_PER_MOLECULE == NATOM
      try:
         if (sum(self.parm.parm_data['ATOMS_PER_MOLECULE']) != 
             self.parm.ptr('natom')):
            print >> sys.stderr, 'sum(ATOMS_PER_MOLECULE) != NATOM'
            valid = False
      except KeyError:
         pass

      # Duplicate pmemd's checks
      if self.parm.ptr('ifpert') != 0:
         print >> sys.stderr, 'IFPERT must be 0!'
         valid = False
      if self.parm.ptr('mbona') != self.parm.ptr('nbona') or \
         self.parm.ptr('mtheta') != self.parm.ptr('ntheta') or \
         self.parm.ptr('mphia') != self.parm.ptr('nphia'):
         print >> sys.stderr, 'Constraints can no longer be put in the prmtop!'
         valid = False

      return valid

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class change(Action):
   """Changes the property of given atoms to a new value. <property> can be
      CHARGE, MASS, RADII, SCREEN, ATOM_NAME, or AMBER_ATOM_TYPE
   """
   def __init__(self, parm, prop, atom_mask, new_val):
      Action.__init__(self, parm)
      self.prop = prop.upper()
      self.mask = AmberMask(self.parm, atom_mask)
      if self.prop in ('CHARGE', 'RADII', 'SCREEN', 'MASS'):
         self.new_val = float(new_val)
         self.new_val_str = '%.4f' % self.new_val
      elif self.prop in ('ATOM_NAME', 'AMBER_ATOM_TYPE'):
         self.new_val = new_val[:4]
         if len(new_val) > 4:
            print >> sys.stderr, ('Warning: Only 4 letters allowed for %s ' +
                     'entries! Truncating remaining letters.') % (self.prop)
         self.new_val_str = '%-4s' % self.new_val
      else:
         raise ParmedChangeError('You may only use "change" with CHARGE, ' +
            'MASS, RADII, SCREEN, ATOM_NAME, or AMBER_ATOM_TYPE!')

   def __str__(self):
      atnums = self.mask.Selection()
      if sum(atnums) == 0:
         return "change %s: Nothing to do" % self.prop
      string = '\n'
      for i in range(self.parm.ptr('natom')):
         if atnums[i] == 1:
            string += "Changing %s of atom # %d (%s) from %s to %s\n" % (
               self.prop, i+1, self.parm.parm_data['ATOM_NAME'][i],
               self.parm.parm_data[self.prop][i], self.new_val_str)
      return string

   def execute(self):
      atnums = self.mask.Selection()
      if sum(atnums) == 0:
         print 'change %s: %s matches no atoms' % (self.prop, self.mask)
         return
      for i in range(len(atnums)):
         if atnums[i] == 1:
            self.parm.parm_data[self.prop][i] = self.new_val

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printinfo(Action):
   """ Prints all prmtop data corresponding to the given %FLAG """
   outfile = sys.stdout
   def __init__(self, parm, flag):
      Action.__init__(self, parm)
      self.flag = flag.upper() # FLAGs are always upper-case
      if not self.flag in self.parm.flag_list:
         print '%%FLAG %s not found!' % self.flag
         self.found = False
      else:
         if 'E' in self.parm.formats[self.flag].upper() \
               or 'F' in self.parm.formats[self.flag.upper()]:
            num = 0
            for item in self.parm.parm_data[self.flag]:
               if num % 5 == 0: self.outfile.write('\n')
               self.outfile.write('%16.5f ' % item)
               num += 1
         else:
            num = 0
            for item in self.parm.parm_data[self.flag]:
               if num % 5 == 0: self.outfile.write('\n')
               self.outfile.write('%-16s ' % item)
               num += 1

      self.outfile.write('\n')
      self.found = True

   def __str__(self):
      return 'Print all data for %%FLAG %s' % self.flag

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class addljtype(Action):
   """Turns given mask into a new LJ atom type. It uses the radius and Rmin from
      the first atom type in <mask> if new_radius or new_epsilon aren't provided
   """
   def __init__(self, parm, mask, new_radius=None, new_epsilon=None, 
                new_radius_14=None, new_epsilon_14=None):
      Action.__init__(self, parm)
      self.mask = AmberMask(self.parm, mask)
      self.new_radius = new_radius
      self.new_epsilon = new_epsilon
      self.new_radius_14 = new_radius_14
      self.new_epsilon_14 = new_epsilon_14

   def __str__(self):
      return 'Making atoms %s into a new LJ atom type' % self.mask

   def execute(self):
      from ParmedTools.addljtype import AddLJType
      # Find the first atom that's selected in this selection. We've
      # already made sure that at least one atom was selected
      sel_atms = self.mask.Selection()
      for i in range(self.parm.ptr('natom')):
         if sel_atms[i] == 1: 
            first_atm = i
            break
      # If either the radius or epsilon were not specified, then pull it
      # from the *first* atom in the selection
      if self.new_radius == None:
         self.new_radius = self.parm.LJ_radius[
             self.parm.parm_data['ATOM_TYPE_INDEX'][first_atm]-1]
      else:
         self.new_radius = float(self.new_radius)
      if self.new_epsilon == None:
         self.new_epsilon = self.parm.LJ_depth[
               self.parm.parm_data['ATOM_TYPE_INDEX'][first_atm]-1]
      else:
         self.new_epsilon = float(self.new_epsilon)
      # Now do the same for chamber prmtops
      if self.new_radius_14 == None and self.parm.chamber:
         self.new_radius_14 = self.parm.LJ_14_radius[
               self.parm.parm_data['ATOM_TYPE_INDEX'][first_atm]-1]
      elif self.parm.chamber:
         self.new_radius_14 = self.new_radius_14
      if self.new_epsilon_14 == None and self.parm.chamber:
         self.new_epsilon_14 = self.parm.LJ_14_depth[
               self.parm.parm_data['ATOM_TYPE_INDEX'][first_atm]-1]
      elif self.parm.chamber:
         self.new_epsilon_14 = self.new_epsilon_14
      
      AddLJType(self.parm, sel_atms, self.new_radius, self.new_epsilon, 
                self.new_radius_14, self.new_epsilon_14)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class outparm(Action):
   """ Prints a new prmtop like parmout, but keeps its place in the action stack
       so several can be written out in 1 parmed session
   """
   def __init__(self, parm, filename, rst_name=None):
      Action.__init__(self, parm)
      self.filename = filename
      self.rst_name = rst_name

   def __str__(self):
      return 'Outputting Amber topology file %s' % self.filename

   def execute(self):
      self.parm.writeParm(self.filename)
      if self.rst_name:
         try: self.parm.writeRst7(self.rst_name)
         except IndexError:
            raise ParmError('You must load coordinates to write a restart')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printljtypes(Action):
   """ Prints the Lennard Jones type index for the given atom mask or, if no 
       value is given, the LJ type index for each atom.
   """
   def __init__(self, parm, idx=None):
      Action.__init__(self, parm)
      # Compile the list of indices
      if idx == None:
         print 'Printing all nonbonded indices'
         self.idx = [i for i in range(1,self.parm.ptr('ntypes')+1)]
         print self.idx
      elif '@' in idx or ':' in idx: # this is a mask
         self.mask = AmberMask(self.parm, idx)
         self.type_list = None
      else:
         self.mask = None
         self.type_list = []
         type_fields = idx.strip().split(',')
         for field in type_fields:
            if len(field.strip()) == 0: continue
            if '-' in field:
               begin = int(field.split('-')[0])
               end = min(int(field.split('-')[1]), self.parm.ptr('ntypes'))
               if begin < 0 or end < begin: 
                  raise ParmError('printLJTypes: Bad atom type range')
               self.type_list.extend([i for i in range(begin, end+1)])
            else:
               self.type_list.append(int(field))

   def __str__(self):
      # Construct the atom selections and related atom types lists
      if self.mask:
         selection = self.mask.Selection()
      elif self.type_list:
         selection = [0 for i in range(self.parm.ptr('natom'))]
         for item in self.type_list:
            selection[item-1] = 1
      else: return 'Nothing to do for printLJTypes'

      self.idx = []

      for i in range(self.parm.ptr('natom')):
         if selection[i] == 1:
            if not self.parm.parm_data['ATOM_TYPE_INDEX'][i] in self.idx:
               self.idx.append(self.parm.parm_data['ATOM_TYPE_INDEX'][i])

      string = '\n%15s %4s %4s\n' % ("  ATOM NUMBER  ", 'NAME', 'TYPE')
      string += '---------------------------------------------\n'
      for i in range(self.parm.ptr('natom')):
         if self.parm.parm_data['ATOM_TYPE_INDEX'][i] in self.idx:
            string += 'ATOM %-10d %-4s %-4s: Type index: %d\n' % (
                                      i+1, self.parm.parm_data['ATOM_NAME'][i],
                                      self.parm.parm_data['AMBER_ATOM_TYPE'][i],
                                      self.parm.parm_data['ATOM_TYPE_INDEX'][i])

      return string

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class scee(Action):
   """ Sets the 1-4 EEL scaling factor in the prmtop """
   def __init__(self, parm, value):
      Action.__init__(self, parm)
      self.scee_value = float(value)

   def __str__(self):
      return ("Setting default SCEE electrostatic scaling value to %.4f" % 
              self.scee_value)

   def execute(self):
      if not 'SCEE_SCALE_FACTOR' in self.parm.flag_list:
         self.parm.addFlag('SCEE_SCALE_FACTOR', '5E16.8', data=[self.scee_value
                                        for i in range(self.parm.ptr('nptra'))])
      else:
         self.parm.parm_data['SCEE_SCALE_FACTOR'] = [self.scee_value 
                                         for i in range(self.parm.ptr('nptra'))]

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class scnb(Action):
   """ Sets the 1-4 VDW scaling factor in the prmtop """
   def __init__(self, parm, value):
      Action.__init__(self, parm)
      self.scnb_value = float(value)

   def __str__(self):
      return "Setting default SCNB van der Waals scaling value to %.4f" % \
             self.scnb_value

   def execute(self):
      if not 'SCNB_SCALE_FACTOR' in self.parm.flag_list:
         self.parm.addFlag('SCNB_SCALE_FACTOR','5E16.8', data=[self.scnb_value 
                                        for i in range(self.parm.ptr('nptra'))])
      else:
         self.parm.parm_data['SCNB_SCALE_FACTOR'] = [self.scnb_value 
                                         for i in range(self.parm.ptr('nptra'))]

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changeljsingletype(Action):
   """ Allows you to change the radius/well depth of a single LJ type
       specified by <mask>
   """
   def __init__(self, parm, mask, radius, depth):
      Action.__init__(self, parm)
      self.mask = AmberMask(self.parm, mask)
      self.radius = float(radius)
      self.depth = float(depth)
      self.orig_radius = self.orig_depth = None
      if 1 in self.mask.Selection():
         first_loc = self.mask.Selection().index(1)
         attype = self.parm.parm_data['ATOM_TYPE_INDEX'][first_loc]
         self.orig_radius = self.parm.LJ_radius[attype-1]
         self.orig_depth = self.parm.LJ_depth[attype-1]

   def __str__(self):
      return ("Changing %s Lennard-Jones well depth from %.4f to %.4f " + 
             "(kal/mol) and radius from %.4f to %.4f (Angstroms)") % (self.mask,
             self.orig_depth, self.depth, self.orig_radius, self.radius)

   def execute(self):
      from math import sqrt
      from ParmedTools.exceptions import LJ_TypeError
      # If this is an empty mask do nothing
      if self.orig_radius is None: return
      # Make sure we've only selected a single atom type with our mask
      attype = None
      for i, sel in enumerate(self.mask.Selection()):
         if sel == 1:
            if attype is None:
               attype = self.parm.parm_data['ATOM_TYPE_INDEX'][i]
            else:
               if attype != self.parm.parm_data['ATOM_TYPE_INDEX'][i]:
                  raise LJ_TypeError('changeLJSingleType: ' +
                          'Selection mask has multiple atom types!')
      # Fill the Lennard-Jones radius and depth arrays to make sure they're
      # up-to-date
      self.parm.fill_LJ()
      self.parm.LJ_radius[attype-1] = self.radius
      self.parm.LJ_depth[attype-1] = self.depth

      for i in range(self.parm.ptr('ntypes')):
         lj_index = self.parm.parm_data['NONBONDED_PARM_INDEX'][
            self.parm.ptr('ntypes') * i + attype-1] - 1
         rij = self.parm.LJ_radius[i] + self.radius
         wij = sqrt(self.parm.LJ_depth[i] * self.depth)
         acoef = wij * rij ** 12
         bcoef = 2 * wij * rij ** 6
         self.parm.parm_data['LENNARD_JONES_ACOEF'][lj_index] = acoef
         self.parm.parm_data['LENNARD_JONES_BCOEF'][lj_index] = bcoef

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printdetails(Action):
   """ Returns information about all atoms in a given mask """
   def __init__(self, parm, mask):
      Action.__init__(self, parm)
      self.mask = AmberMask(self.parm, mask)

   def __str__(self):
      selection = self.mask.Selection()
      retstr = "\nThe mask %s matches %d atoms:\n\n"%(self.mask,sum(selection))
      retstr += "%7s%7s%9s%6s%6s%12s%12s%10s%10s%10s%10s\n" % ('ATOM',
                'RES','RESNAME','NAME','TYPE','LJ Radius','LJ Depth','Mass',
                'Charge','GB Radius','GB Screen')
      for i in range(self.parm.ptr('natom')):
         if selection[i] == 1:
            retstr += "%7d%7d%9s%6s%6s%12.4f%12.4f%10.4f%10.4f%10.4f%10.4f\n"%(
              i+1, self.parm.residue_container[i],
              self.parm.parm_data['RESIDUE_LABEL'][
                                             self.parm.residue_container[i]-1],
              self.parm.parm_data['ATOM_NAME'][i], 
              self.parm.parm_data['AMBER_ATOM_TYPE'][i],
              self.parm.LJ_radius[self.parm.parm_data['ATOM_TYPE_INDEX'][i]-1],
              self.parm.LJ_depth[self.parm.parm_data['ATOM_TYPE_INDEX'][i]-1],
              self.parm.parm_data['MASS'][i], self.parm.parm_data['CHARGE'][i],
              self.parm.parm_data['RADII'][i], self.parm.parm_data['SCREEN'][i])
      return retstr

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printflags(Action):
   """ Prints all %FLAGs found in the topology file """
   def __str__(self):
      string = '\n'
      for flag in self.parm.flag_list:
         string += '%%FLAG %s\n' % flag
      return string

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printpointers(Action):
   """ Prints a list of all the POINTERS and their values """
   def __str__(self):
      ptrs = self.parm.parm_data['POINTERS']
      ret_str = """\nNATOM (number of atoms in system)................= %d
NTYPES (number of atom type names)...............= %d
NBONH (number of bonds containing H).............= %d
MBONA (number of bonds without H)................= %d
NTHETH (number of angles containing H)...........= %d
MTHETA (number of angles without H)..............= %d
NPHIH (number of dihedrals containing H).........= %d
MPHIA (number of dihedrals without H)............= %d
NHPARM (currently unused)........................= %d
NPARM (1 if made with addles, 0 if not)..........= %d
NNB (number of excluded atoms)...................= %d
NRES (number of residues in system)..............= %d
NBONA (MBONA + constraint bonds).................= %d
NTHETA (MTHETA + constraint angles)..............= %d
NPHIA (MPHIA + constraint dihedrals).............= %d
NUMBND (number of unique bond types).............= %d
NUMANG (number of unique angle types)............= %d
NPTRA (number of unique dihedral types)..........= %d
NATYP (number of nonbonded atom types)...........= %d
NPHB (number of distinct 10-12 H-bond pairs).....= %d
IFPERT (1 if prmtop is perturbed; not used)......= %d
NBPER (perturbed bonds; not used)................= %d
NGPER (perturbed angles; not used)...............= %d
NDPER (perturbed dihedrals; not used)............= %d
MBPER (bonds in perturbed group; not used).......= %d
MGPER (angles in perturbed group; not used)......= %d
MDPER (diheds in perturbed group; not used)......= %d
IFBOX (Type of box: 1=orthogonal, 2=not, 0=none).= %d
NMXRS (number of atoms in largest residue).......= %d
IFCAP (1 if solvent cap exists)..................= %d
NUMEXTRA (number of extra points in topology)....= %d
""" % (
                       ptrs[0], ptrs[1], ptrs[2], ptrs[3], ptrs[4], ptrs[5], 
                       ptrs[6], ptrs[7], ptrs[8], ptrs[9], ptrs[10], ptrs[11], 
                       ptrs[12], ptrs[13], ptrs[14], ptrs[15], ptrs[16], 
                       ptrs[17], ptrs[18], ptrs[19], ptrs[20], ptrs[21], 
                       ptrs[22], ptrs[23], ptrs[24], ptrs[25], ptrs[26], 
                       ptrs[27], ptrs[28], ptrs[29], ptrs[30] )
      if len(ptrs) == 32: ret_str += \
         "NCOPY (number of PIMD slices/number of beads)....= %d\n" % ptrs[31]
      if self.parm.ptr('IFBOX'):
         ret_str += "\nSOLVENT POINTERS\n" + """
IPTRES (Final solute residue)....................= %d
NSPM (Total number of molecules).................= %d
NSPSOL (The first solvent "molecule")............= %d
""" % (self.parm.parm_data['SOLVENT_POINTERS'][0],
       self.parm.parm_data['SOLVENT_POINTERS'][1],
       self.parm.parm_data['SOLVENT_POINTERS'][2])
         
      return ret_str

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class setmolecules(Action):
   """Determines the molecularity of the system based on the bonding network and
      correctly determines the SOLVENT_POINTERS and ATOMS_PER_MOLECULE sections
      of the topology file. It will consider the ions to be part of the solute
      if True is passed or not if False is passed. Defaults to False.
   """
   def __init__(self, parm, solute_ions=False):
      Action.__init__(self, parm)
      if solute_ions == 'True' or solute_ions == 'solute_ions=True':
         self.solute_ions = True
      elif solute_ions == 'False' or solute_ions == 'solute_ions=False':
         self.solute_ions = False
      else:
         raise ParmedMoleculeError('Illegal value for solute_ions')

   def __str__(self):
      return "Setting MOLECULE properties of the prmtop (SOLVENT_POINTERS "+ \
             "and ATOMS_PER_MOLECULE)"

   def execute(self):
      self.parm.rediscover_molecules(self.solute_ions)
   
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class combinemolecules(Action):
   """Combines the molecule <mol_id1> with the adjacent molecule <mol_id2> to
      ensure that those two molecules are imaged together. Most commonly used
      with protein/small ligand systems to keep them together during wrapping. 
      This will slightly affect the pressure calculation, but not by very much
      for small ligands. Note that <mol_id1> and <mol_id2> must be sequential
      if <mol_id2> is supplied
   """
   def __init__(self, parm, mol_id, mol_id2=None):
      Action.__init__(self, parm)
      self.mol_id = int(mol_id)
      if mol_id2:
         if self.mol_id + 1 != int(mol_id2):
            raise ParmedMoleculeError('Can only combine adjacent molecules!')

   def __str__(self):
      return "Combining molecules %d and %d into the same molecule" % (
              self.mol_id, self.mol_id + 1)

   def execute(self):
      from ParmedTools.mod_molecules import combineMolecules as cm
      cm(self.parm, self.mol_id)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class quit(Action):
   """ Stops where we are without executing any parmout command """
   def __init__(self,parm):
      import sys
      print 'Quitting.'
      sys.exit(0)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#class addcoarsegrain(Action):
#   """ Adds coarse graining information to the Amber topology file according to
#       a given coarse graining parameter file
#   """
#   def __init__(self, parm, cg_param_file):
#      from os.path import exists
#      Action.__init__(self, parm)
#      if not exists(cg_param_file):
#         raise CoarseGrainError('Cannot find parameter file %s' % cg_param_file)
#      self.cg_param_file = cg_param_file
#      # Check to see if we've already made this a coarsegrained file...
#      if 'ANGLE_COEF_A' in self.parm.flag_list:
#         raise CoarseGrainError('Prmtop already has coarse grained sections')
#
#   def __str__(self):
#      return ("Setting up coarse graining for topology file using parameter " +
#              "file " + self.cg_param_file)
#
#   def execute(self):
#      from ParmedTools.coarsegrain import addCoarseGrain as addCG
#      addCG(self.parm, self.cg_param_file)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class changeprotstate(Action):
   """Changes the protonation state of a given titratable residue that can be 
      treated via constant pH MD in Amber. This corresponds to all residues 
      found in $AMBERHOME/AmberTools/src/etc/cpin_data.py
   """
   def __init__(self, parm, mask, state):
      Action.__init__(self, parm)
      self.state = int(state)
      self.mask = AmberMask(self.parm, mask)

   def __str__(self):
      sel = self.mask.Selection()
      if sum(sel) == 0:
         return "No residues selected for state change"
      res = self.parm.residue_container[sel.index(1)]
      return 'Changing protonation state of residue %d (%s) to %d' % (res,
         self.parm.parm_data['RESIDUE_LABEL'][res-1], self.state)
   
   def execute(self):
      from cpinutils import cpin_data
      sel = self.mask.Selection()
      # If we didn't select any residues, just return
      if sum(sel) == 0: return
      resnum = self.parm.residue_container[sel.index(1)]
      resname = self.parm.parm_data['RESIDUE_LABEL'][resnum-1]
      # Get the charges from cpin_data. The first 2 elements are energy and 
      # proton count so the charges are chgs[2:]
      chgs = cpin_data.getData(resname, 2)
      # First catch any errors
      if chgs == -1: 
         raise ChangeStateError("Residue %s isn't defined as a titratable " +
                                "residue in cpin_data.py")

      if self.state > len(chgs):
         raise ChangeStateError(('Residue %s only has %d titratable states. ' +
                       'You chose state %d') % (resname, len(chgs), self.state))

      if sum(sel) != len(chgs[self.state])-2:
         raise ChangeStateError('You must select one and only one entire ' +
                                'titratable residue')
      
      chgnum = 2
      for i in range(self.parm.parm_data['RESIDUE_POINTER'][resnum-1]-1,
                     self.parm.parm_data['RESIDUE_POINTER'][resnum]-1):
         if sel[i] != 1:
            raise ChangeStateError('You must select 1 and only 1 entire ' +
                                   'residue to change the protonation state of')
         # Print some info about the change we're making
         print ('  Changing state of atom %-4s (residue %-4s) from %8.4f to ' +
                '%8.4f') % ( self.parm.parm_data['ATOM_NAME'][i], resname, 
            self.parm.parm_data['CHARGE'][i],  chgs[self.state][chgnum])
         # Actually make the change
         self.parm.parm_data['CHARGE'][i] = chgs[self.state][chgnum]
         chgnum += 1

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class netcharge(Action):
   """ Prints the total charge of all of the atoms given by the mask.
       Defaults to all atoms
   """
   def __init__(self, parm, mask=':*'):
      Action.__init__(self, parm)
      self.mask = AmberMask(self.parm, mask)

   def execute(self):
      """ Calculates the charge of all atoms selected in mask """
      sel = self.mask.Selection()

      netchg = 0.0
      for i in range(len(sel)):
         if sel[i]: netchg += self.parm.parm_data['CHARGE'][i]

      print '  The net charge of %s is %.4f' % (self.mask, netchg)
      return netchg

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class strip(Action):
   """ Deletes the atoms specified by <mask> from the topology file and rebuilds
       the topology file according to the parameters that remain.
   """
   def __init__(self, parm, mask):
      Action.__init__(self, parm)
      self.mask = AmberMask(self.parm, mask)
      self.num_atms = sum(self.mask.Selection())
   
   def __str__(self):
      return "Removing mask '%s' (%d atoms) from the topology file." % (
         self.mask, self.num_atms)

   def execute(self):
      self.parm.delete_mask(self.mask)

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class definesolvent(Action):
   """ Allows you to change what parmed.py will consider to be "solvent".
       <residue list> must be a comma-separated set of residue names with
       no spaces between them.
   """
   def __init__(self, parm, res_list):
      res_list.replace(' ', '')
      Action.__init__(self, parm)
      if res_list.endswith(','): self.res_list = res_list[:len(res_list)-1]
      else: self.res_list = res_list
      self.parm.solvent_residues = res_list.split(',')

   def __str__(self):
      return "Residues %s are now considered to be solvent" % self.res_list

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class addexclusions(Action):
   """ 
   Allows you to add arbitrary exclusions to the exclusion list. Every atom in
   <mask2> is added to the exclusion list for each atom in <mask1> so that
   non-bonded interactions between those atom pairs will not be computed. NOTE
   that this ONLY applies to direct-space (short-range) non-bonded potentials.
   For PME simulations, long-range electrostatics between these atom pairs are
   still computed (in different unit cells).
   """
   def __init__(self, parm, mask1, mask2):
      Action.__init__(self, parm)
      self.mask1 = AmberMask(self.parm, mask1)
      self.mask2 = AmberMask(self.parm, mask2)

   def __str__(self):
      return 'Adding atoms from %s to exclusion lists of atoms in %s' % (
         self.mask2, self.mask1)
   
   def execute(self):
      sel1 = self.mask1.Selection()
      sel2 = self.mask2.Selection()
      # Loop through both selections and add each selected atom in sel2 to the
      # exclusion list for selected atoms in sel1 (and vice-versa).
      for i in range(len(sel1)):
         if sel1[i]:
            for j in range(len(sel2)):
               if sel2[j]:
                  # Make sure that this atom isn't already in the exclusion list
                  # by virtue of being a bonded partner.
                  atm1 = self.parm.atom_list[i]
                  atm2 = self.parm.atom_list[j]
                  # Skip over atm1 == atm2
                  if atm1 is atm2: continue
                  # Add each other to each other's exclusion lists.
                  atm1.exclude(atm2)
                  self.parm.atom_list.changed = True

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printbonds(Action):
   """
   Prints all of the bonds (with their details) for the given atoms in the
   mask
   """
   def __init__(self, parm, mask):
      Action.__init__(self, parm)
      self.mask = AmberMask(self.parm, mask)

   def __str__(self):
      retstr = '%-20s %-20s %-10s %-10s\n' % (
               'Atom 1', 'Atom 2', 'R eq', 'Frc Cnst')
      # Loop through all of the bonds without and inc hydrogen
      atomsel = self.mask.Selection()
      for bond in self.parm.bonds_without_h:
         idx1 = bond.atom1.starting_index
         idx2 = bond.atom2.starting_index
         if not (atomsel[idx1] or atomsel[idx2]): continue

         atm1 = self.parm.atom_list[idx1]
         atm2 = self.parm.atom_list[idx2]
         retstr += '%7d %4s (%4s) %7d %4s (%4s) %10.4f %10.4f\n' % (
            idx1+1, atm1.atname, atm1.attype, idx2+1, atm2.atname,
            atm2.attype, bond.bond_type.req, bond.bond_type.k)

      for bond in self.parm.bonds_inc_h:
         idx1 = bond.atom1.starting_index
         idx2 = bond.atom2.starting_index
         if not (atomsel[idx1] or atomsel[idx2]): continue

         atm1 = self.parm.atom_list[idx1]
         atm2 = self.parm.atom_list[idx2]
         retstr += '%7d %4s (%4s)  %7d %4s (%4s) %10.4f %10.4f\n' % (
            idx1+1, atm1.atname, atm1.attype, idx2+1, atm2.atname,
            atm2.attype, bond.bond_type.req, bond.bond_type.k)
      
      return retstr

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printangles(Action):
   """
   Prints all of the angles (with their details) for the given atoms in the
   mask
   """
   def __init__(self, parm, mask):
      Action.__init__(self, parm)
      self.mask = AmberMask(self.parm, mask)

   def __str__(self):
      retstr = '%-20s %-20s %-20s %-10s %-10s\n' % (
               'Atom 1', 'Atom 2', 'Atom 3', 'Frc Cnst', 'Theta eq')
      # Loop through all of the bonds without and inc hydrogen
      atomsel = self.mask.Selection()
      for angle in self.parm.angles_without_h:
         idx1 = angle.atom1.starting_index
         idx2 = angle.atom2.starting_index
         idx3 = angle.atom3.starting_index
         if not (atomsel[idx1] or atomsel[idx2] or atomsel[idx3]): continue

         atm1 = self.parm.atom_list[idx1]
         atm2 = self.parm.atom_list[idx2]
         atm3 = self.parm.atom_list[idx3]
         retstr += \
            '%7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s) %10.4f %10.4f\n' % (
            idx1+1, atm1.atname, atm1.attype, idx2+1, atm2.atname,
            atm2.attype, idx3+1, atm3.atname, atm3.attype,
            angle.angle_type.k, angle.angle_type.theteq*180/math.pi)

      for angle in self.parm.angles_inc_h:
         idx1 = angle.atom1.starting_index
         idx2 = angle.atom2.starting_index
         idx3 = angle.atom3.starting_index
         if not (atomsel[idx1] or atomsel[idx2] or atomsel[idx3]): continue

         atm1 = self.parm.atom_list[idx1]
         atm2 = self.parm.atom_list[idx2]
         atm3 = self.parm.atom_list[idx3]
         retstr += \
            '%7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s) %10.4f %10.4f\n' % (
            idx1+1, atm1.atname, atm1.attype, idx2+1, atm2.atname,
            atm2.attype, idx3+1, atm3.atname, atm3.attype,
            angle.angle_type.k, angle.angle_type.theteq*180/math.pi)
      
      return retstr

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class printdihedrals(Action):
   """
   Prints all of the dihedrals (with their details) for the given atoms in the
   mask
   """
   def __init__(self, parm, mask):
      Action.__init__(self, parm)
      self.mask = AmberMask(self.parm, mask)

   def __str__(self):
      retstr = '%-20s %-20s %-20s %-20s  %-10s %-10s %-10s %-10s %-10s\n' % (
               'Atom 1', 'Atom 2', 'Atom 3', 'Atom 4', 'Height', 'Periodic.',
               'Phase', 'EEL Scale', 'VDW Scale')
      # Loop through all of the bonds without and inc hydrogen
      atomsel = self.mask.Selection()
      for dihedral in self.parm.dihedrals_without_h:
         idx1 = dihedral.atom1.starting_index
         idx2 = dihedral.atom2.starting_index
         idx3 = dihedral.atom3.starting_index
         idx4 = dihedral.atom4.starting_index
         if not (atomsel[idx1] or atomsel[idx2] or atomsel[idx3] or 
                 atomsel[idx4]): continue

         atm1 = self.parm.atom_list[idx1]
         atm2 = self.parm.atom_list[idx2]
         atm3 = self.parm.atom_list[idx3]
         atm4 = self.parm.atom_list[idx4]
         # Determine if it's an Improper, Multiterm, or neither
         if dihedral.signs[1] < 0: char = 'I'
         elif dihedral.signs[0] < 0: char = 'M'
         else: char = ' '
         retstr += \
            ('%1s %7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s) ' +
            '%10.4f %10.4f %10.4f %10.4f %10.4f\n') % (char,
            idx1+1, atm1.atname, atm1.attype, idx2+1, atm2.atname, atm2.attype,
            idx3+1, atm3.atname, atm3.attype, idx4+1, atm4.atname, atm4.attype,
            dihedral.dihed_type.phi_k, dihedral.dihed_type.per,
            dihedral.dihed_type.phase*180/math.pi, dihedral.dihed_type.scee,
            dihedral.dihed_type.scnb)

      for dihedral in self.parm.dihedrals_inc_h:
         idx1 = dihedral.atom1.starting_index
         idx2 = dihedral.atom2.starting_index
         idx3 = dihedral.atom3.starting_index
         idx4 = dihedral.atom4.starting_index
         if not (atomsel[idx1] or atomsel[idx2] or atomsel[idx3] or 
                 atomsel[idx4]): continue

         atm1 = self.parm.atom_list[idx1]
         atm2 = self.parm.atom_list[idx2]
         atm3 = self.parm.atom_list[idx3]
         atm4 = self.parm.atom_list[idx4]
         if dihedral.signs[1] < 0: char = 'I'
         elif dihedral.signs[0] < 0: char = 'M'
         else: char = ' '
         retstr += \
            ('%1s %7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s)  %7d %4s (%4s) ' +
            '%10.4f %10.4f %10.4f %10.4f %10.4f\n') % (char,
            idx1+1, atm1.atname, atm1.attype, idx2+1, atm2.atname, atm2.attype,
            idx3+1, atm3.atname, atm3.attype, idx4+1, atm4.atname, atm4.attype,
            dihedral.dihed_type.phi_k, dihedral.dihed_type.per,
            dihedral.dihed_type.phase*180/math.pi, dihedral.dihed_type.scee,
            dihedral.dihed_type.scnb)

      return retstr

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class setbond(Action):
   """
   Changes (or adds a non-existent) bond in the topology file. Each mask must
   select the same number of atoms, and a bond will be placed between the
   atoms in mask1 and mask2 (one bond between atom1 from mask1 and atom1
   from mask2 and another bond between atom2 from mask1 and atom2 from mask2,
   etc.)
   """

   def __init__(self, parm, mask1, mask2, k, req):
      Action.__init__(self, parm)
      self.mask1 = AmberMask(self.parm, mask1)
      self.mask2 = AmberMask(self.parm, mask2)
      self.k = float(k)
      self.req = float(req)
   
   def __str__(self):
      return ('Set a bond between %s and %s with k = %f kcal/(mol Angstrom' + 
              '**2) and Req = %f Angstroms') % (self.mask1, self.mask2,
              self.k, self.req)

   def execute(self):
      sel1 = self.mask1.Selection()
      sel2 = self.mask2.Selection()

      if sum(sel1) != sum(sel2):
         raise SetParamError('setBond: Each mask must select the same number ' +
                             'of atoms!')

      # If no atoms, nothing to do
      if sum(sel1) == 0: return
   
      # Create the new bond type
      new_bnd_typ = _BondType(self.k, self.req, -1)
      # Does that bond type exist in the list already? If it does, re-bind
      # new_bnd to that bond type reference
      exists = False
      for bnd_typ in self.parm.bond_type_list:
         if new_bnd_typ == bnd_typ:
            new_bnd_typ = bnd_typ
            exists = True
            break
      # If the bond is new, add it to the type list
      if not exists:
         self.parm.bond_type_list.append(new_bnd_typ)

      atnum1, atnum2 = -1, -1
      # Loop through all of the selected atoms
      for it in range(sum(sel1)):
         # Collect the atoms involved
         atnum1 = sel1.index(1, atnum1+1)
         atnum2 = sel2.index(1, atnum2+1)
         atm1 = self.parm.atom_list[atnum1]
         atm2 = self.parm.atom_list[atnum2]

         # See if any atom is Hydrogen (allow for deuteriums)
         has_h = False
         if ((atm1.atname[0] in ['H', 'h'] and abs(atm1.mass - 1.01) < 1.1) or
             (atm2.atname[0] in ['H', 'h'] and abs(atm2.mass - 1.01) < 1.1)):
            has_h = True
   
         # See if the bond exists in the first place, and if so, replace its bond
         # type with our new bond type (new_bnd)
         if atm2 in atm1.bond_partners or atm1 in atm2.bond_partners:
            if has_h:
               bond_list = self.parm.bonds_inc_h
            else:
               bond_list = self.parm.bonds_without_h
   
            for bnd in bond_list:
               if atm1 in bnd and atm2 in bnd:
                  bnd.bond_type = new_bnd_typ
                  bond_list.changed = True
                  break
   
         # Otherwise, it doesn't exist, so we just create a new one
         else:
            if has_h:
               self.parm.bonds_inc_h.append(_Bond(atm1, atm2, new_bnd_typ))
            else:
               self.parm.bonds_without_h.append(_Bond(atm1, atm2, new_bnd_typ))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class setangle(Action):
   """
   Changes (or adds a non-existent) angle in the topology file. Each mask must
   select the same number of atoms, and an angle will be placed between the
   atoms in mask1, mask2, and mask3 (one angle between atom1 from mask1, atom1
   from mask2, and atom1 from mask3, another angle between atom2 from mask1,
   atom2 from mask2, and atom2 from mask3, etc.)
   """

   def __init__(self, parm, mask1, mask2, mask3, k, theteq):
      Action.__init__(self, parm)
      self.mask1 = AmberMask(self.parm, mask1)
      self.mask2 = AmberMask(self.parm, mask2)
      self.mask3 = AmberMask(self.parm, mask3)
      self.k = float(k)
      self.theteq = float(theteq) * math.pi / 180.0
   
   def __str__(self):
      return ('Set an angle between %s, %s and %s with k = %f kcal/(mol ' +
              'rad**2) and THETeq = %f degrees') % (self.mask1, self.mask2,
              self.mask3, self.k, self.theteq * 180/math.pi)

   def execute(self):
      sel1 = self.mask1.Selection()
      sel2 = self.mask2.Selection()
      sel3 = self.mask3.Selection()

      if sum(sel1) != sum(sel2) or sum(sel1) != sum(sel3):
         raise SetParamError('Each mask in setAngle must select the same ' +
                             'number of atoms!')

      if sum(sel1) == 0: return

      # Create the new angle type
      new_ang_typ = _AngleType(self.k, self.theteq, -1)
      # Does that angle type exist in the list already? If it does, re-bind
      # new_ang to that angle type reference
      exists = False
      for ang_typ in self.parm.angle_type_list:
         if new_ang_typ == ang_typ:
            new_ang_typ = ang_typ
            exists = True
            break
      # If the angle is new, add it to the type list
      if not exists:
         self.parm.angle_type_list.append(new_ang_typ)

      atnum1, atnum2, atnum3 = -1, -1, -1

      # Loop through all of the selections
      for it in range(sum(sel1)):
         # Collect the atoms involved
         atnum1 = sel1.index(1, atnum1+1)
         atnum2 = sel2.index(1, atnum2+1)
         atnum3 = sel3.index(1, atnum3+1)
         atm1 = self.parm.atom_list[atnum1]
         atm2 = self.parm.atom_list[atnum2]
         atm3 = self.parm.atom_list[atnum3]
         # See if any atom is Hydrogen (allow for deuteriums)
         has_h = False
         if ((atm1.atname[0] in ['H', 'h'] and abs(atm1.mass - 1.01) < 1.1) or
             (atm2.atname[0] in ['H', 'h'] and abs(atm2.mass - 1.01) < 1.1) or
             (atm3.atname[0] in ['H', 'h'] and abs(atm3.mass - 1.01) < 1.1)):
            has_h = True
   
         # See if the angle exists in the first place, and if so, replace its angle
         # type with our new angle type (new_ang)
         if ((atm1 in atm2.bond_partners and atm1 in atm3.angle_partners) and
             (atm2 in atm3.bond_partners)):
            if has_h:
               angle_list = self.parm.angles_inc_h
            else:
               angle_list = self.parm.angles_without_h
   
            for ang in angle_list:
               if atm1 in ang and atm2 in ang and atm3 in ang:
                  ang.angle_type = new_ang_typ
                  angle_list.changed = True
                  break
   
         # Otherwise, it doesn't exist, so we just create a new one
         else:
            if has_h:
               self.parm.angles_inc_h.append(_Angle(atm1, atm2, atm3, new_ang_typ))
            else:
               self.parm.angles_without_h.append(_Angle(atm1, atm2, atm3,
                                                        new_ang_typ))
   
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class adddihedral(Action):
   """ Adds a dihedral between mask1, mask2, mask3, and mask4. Each mask must
       specify the same number of atoms, and the dihedral is defined around the 
       bond between atoms in mask 2 and 3. If each mask selects 2 atoms, for
       instance, a dihedral will be placed around atom1 in mask 1, atom1 in 
       mask 2, atom1 in mask 3, and atom1 in mask 4.  A second dihedral will
       be placed around atom2 in mask 1, atom2 in mask 2, atom2 in mask 3, and
       atom2 in mask4. dihed_type can either be "normal", "multiterm", or
       "improper". Note for ring systems of 6 or fewer atoms, you'll need to
       use "multiterm" to avoid double-counting 1-4 interactions for some of
       the dihedrals. <phi_k> is the barrier height of the dihedral term,
       <per> is the periodicity, <phase> is the phase offset, <scee> is the
       1-4 EEL scaling factor for this dihedral (default for AMBER is 1.2,
       default for GLYCAM is 1.0), and <scnb> is the 1-4 VDW scaling factor
       for this dihedral (default for AMBER is 2.0, default for GLYCAM is 1.0)
   """

   def __init__(self, parm, mask1, mask2, mask3, mask4, phi_k, per, phase,
                scee, scnb, dihed_type='normal'):
      Action.__init__(self, parm)
      self.mask1 = AmberMask(self.parm, mask1)
      self.mask2 = AmberMask(self.parm, mask2)
      self.mask3 = AmberMask(self.parm, mask3)
      self.mask4 = AmberMask(self.parm, mask4)
      self.phi_k = float(phi_k)
      self.per = float(per)
      self.phase = float(phase) * math.pi/180.0
      self.scee = float(scee)
      self.scnb = float(scnb)
      if dihed_type.lower() == 'normal'[:len(dihed_type)]:
         self.improper = False
         self.multiterm = False
         self.type = 'a normal'
      elif dihed_type.lower() == 'improper'[:len(dihed_type)]:
         self.improper = True
         self.multiterm = True
         self.type = 'an improper'
      elif dihed_type.lower() == 'multiterm'[:len(dihed_type)]:
         self.improper = False
         self.multiterm = True
         self.type = 'a multiterm'
   
   def __str__(self):
      return ('Set %s dihedral between %s, %s, %s, and %s with phi_k = %f ' +
              'kcal/mol periodicity = %f phase = %f degrees scee = %f ' +
              'scnb = %f') % (self.type,
            self.mask1, self.mask2, self.mask3, self.mask4, self.phi_k,
            self.per, self.phase * 180/math.pi, self.scee, self.scnb)

   def execute(self):
      sel1 = self.mask1.Selection()
      sel2 = self.mask2.Selection()
      sel3 = self.mask3.Selection()
      sel4 = self.mask4.Selection()

      if sum(sel1) != sum(sel2) or sum(sel1) != sum(sel3) or \
                                   sum(sel1) != sum(sel4):
         raise SetParamError('addDihedral: Each mask must select the same ' +
                             'number of atoms!')
      
      # If we do nothing, just return
      if sum(sel1) == 0: return
   
      # Create the new dihedral type
      new_dih_typ = _DihedralType(self.phi_k, self.per, self.phase, self.scee,
                                  self.scnb, -1)
      self.parm.dihedral_type_list.append(new_dih_typ)

      # Loop through all of the atoms
      atnum1, atnum2, atnum3, atnum4 = -1, -1, -1, -1

      for it in range(sum(sel1)):
         # Collect the atoms involved
         atnum1 = sel1.index(1, atnum1+1)
         atnum2 = sel2.index(1, atnum2+1)
         atnum3 = sel3.index(1, atnum3+1)
         atnum4 = sel4.index(1, atnum4+1)
         atm1 = self.parm.atom_list[atnum1]
         atm2 = self.parm.atom_list[atnum2]
         atm3 = self.parm.atom_list[atnum3]
         atm4 = self.parm.atom_list[atnum4]
   
         # Make sure atom 1 doesn't occur in the 3rd or 4th spot since we'll have
         # the -0 effect biting us...
         if atm3.starting_index == 0 or atm4.starting_index == 0:
            atm1, atm2, atm3, atm4 = atm4, atm3, atm2, atm1
   
         # See if any atom is Hydrogen (allow for deuterium)
         has_h = False
         if ((atm1.atname[0] in ['H', 'h'] and abs(atm1.mass - 1.01) < 1.1) or
             (atm2.atname[0] in ['H', 'h'] and abs(atm2.mass - 1.01) < 1.1) or
             (atm3.atname[0] in ['H', 'h'] and abs(atm3.mass - 1.01) < 1.1) or
             (atm4.atname[0] in ['H', 'h'] and abs(atm4.mass - 1.01) < 1.1)):
            has_h = True
   
         # See what signs has to be.  signs is a 2-element list with a 1 or -1
         # for the 3rd or 4th atom index. A -1 for the 3rd means the end groups are
         # not calculated (multiterm, 6-or-lower-membered rings, impropers, etc.)
         # and a -1 for the 4th means improper (-1 for 4th always means -1 for 3rd)
         if self.improper:
            signs = [-1,-1]
         elif self.multiterm:
            signs = [-1,1]
         else:
            signs = [1,1]
         # Create our new dihedral!
         if has_h:
            self.parm.dihedrals_inc_h.append(_Dihedral(atm1, atm2, atm3, atm4,
                                                    new_dih_typ, signs))
         else:
            self.parm.dihedrals_without_h.append(_Dihedral(atm1, atm2, atm3, atm4,
                                                        new_dih_typ, signs))

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class addatomicnumber(Action):
   """
   Adds the atomic number of each atom to a new section titled "ATOMIC_NUMBER"
   in the topology file. Elements are identified by the atomic mass found in
   the MASS section of the topology files.  Elements are matched by picking the
   element whose average atomic mass in the periodic table is closest to each
   atom, which should work appropriately for all isotopes of all atoms, except
   possibly Tritium
   """
   def __init__(self, parm):
      Action.__init__(self, parm)
      if 'ATOMIC_NUMBER' in self.parm.flag_list: self.present = True
      else: self.present = False

   def __str__(self):
      if self.present:
         return 'ATOMIC_NUMBER already in [%s] -- Doing nothing.' % self.parm
      return "Adding ATOMIC_NUMBER to [%s]" % self.parm

   def execute(self):
      from chemistry.periodic_table import AtomicNum
      from chemistry.amber.readparm import Element
      if self.present: return
      self.parm.addFlag('ATOMIC_NUMBER','10I8',num_items=self.parm.ptr('natom'))
      for i, mass in enumerate(self.parm.parm_data['MASS']):
         self.parm.parm_data['ATOMIC_NUMBER'][i] = AtomicNum[Element(mass)]

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class deletedihedral(Action):
   """
   Deletes the dihedral around <mask2> and <mask3> in which the end-groups are
   <mask1> and <mask4>. For multi-term dihedrals, it removes each term.
   """
   def __init__(self, parm, mask1, mask2, mask3, mask4):
      Action.__init__(self, parm)
      self.mask1 = AmberMask(self.parm, mask1)
      self.mask2 = AmberMask(self.parm, mask2)
      self.mask3 = AmberMask(self.parm, mask3)
      self.mask4 = AmberMask(self.parm, mask4)
      if sum(self.mask1.Selection()) != sum(self.mask2.Selection()) or \
         sum(self.mask1.Selection()) != sum(self.mask3.Selection()) or \
         sum(self.mask1.Selection()) != sum(self.mask4.Selection()):
         raise DeleteDihedralError('All masks must select the same number of ' +
               'atoms!. They selected %d, %d, %d, and %d, respectively' % (
               sum(self.mask1.Selection()), sum(self.mask2.Selection()),
               sum(self.mask3.Selection()), sum(self.mask4.Selection())))

   def __str__(self):
      if sum(self.mask1.Selection()) == 0:
         return 'No specified dihedrals to delete'
      return 'Deleting dihedral terms involving [%s]-[%s]-[%s]-[%s]' % (
             self.mask1, self.mask2, self.mask3, self.mask4) + \
             ' (At most %d total, distinct, dihedrals)' % (
             sum(self.mask1.Selection()))

   def execute(self):
      from chemistry.periodic_table import AtomicNum
      from chemistry.amber.readparm import Element
      sel1, sel2 = self.mask1.Selection(), self.mask2.Selection()
      sel3, sel4 = self.mask3.Selection(), self.mask4.Selection()
      # Bail out if we're deleting nothing
      if sum(sel1) == 0: return

      # Keep track of the dihedrals we want to delete from each
      # dihedral list (dihedrals_inc_h, dihedrals_without_h)
      deleting_dihedrals = [[],[]]
      # We have already checked that they are the same number of atoms
      # Now, loop through the atoms and see if any dihedrals match that spec
      idx1 = idx2 = idx3 = idx4 = 0
      total_diheds = 0
      for i in range(sum(sel1)):
         # This helps us keep track of multi-term dihedrals so we don't confuse
         # users
         found_this_dihedral = False
         # Get the first selected atom from each selection then update our index
         idx1, idx2 = sel1.index(1), sel2.index(1)
         idx3, idx4 = sel3.index(1), sel4.index(1)
         # Make sure none of the indices are the same
         if idx1 == idx2 or idx1 == idx3 or idx1 == idx4 or idx2 == idx3 or \
            idx2 == idx4 or idx3 == idx4:
            print 'Skipping %d-%d-%d-%d dihedral deletion -- duplicate atoms!'%(
               idx1, idx2, idx3, idx4)
            continue
         # Figure out if our dihedral would have hydrogen or not (limits what
         # dihedral list we have to search...)
         if AtomicNum[Element(self.parm.parm_data['MASS'][idx1])] == 1 or \
            AtomicNum[Element(self.parm.parm_data['MASS'][idx2])] == 1 or \
            AtomicNum[Element(self.parm.parm_data['MASS'][idx3])] == 1 or \
            AtomicNum[Element(self.parm.parm_data['MASS'][idx4])] == 1:
            dihed_list = self.parm.dihedrals_inc_h
            dihed_list_idx = 0
         else:
            dihed_list = self.parm.dihedrals_without_h
            dihed_list_idx = 1
         # Now search through our dihedral list to see which indexes (if any) we
         # have to remove. Keep tabs of them so we can pop them in reverse order
         # (so we don't have to re-figure indices) afterwards
         proposed_dihedral = (idx1, idx2, idx3, idx4)
         for j, dihed in enumerate(dihed_list):
            if dihed == proposed_dihedral:
               if not found_this_dihedral:
                  print 'Matched dihedral number %d' % j
                  found_this_dihedral = True
                  total_diheds += 1
               else:
                  print '  Matched multi-term dihedral number %d' % j
               deleting_dihedrals[dihed_list_idx].append(j)
      if not deleting_dihedrals[0] and not deleting_dihedrals[1]:
         print 'No dihedrals matched -- not deleting any dihedrals'
         return

      print 'Deleting %d dihedrals' % (len(deleting_dihedrals[0]) + 
                                       len(deleting_dihedrals[1])),
      print ' (%d distinct dihedrals)' % total_diheds

      # At this point, we've collected all of our dihedrals, now sort them
      if deleting_dihedrals[0]: deleting_dihedrals[0].sort()
      if deleting_dihedrals[1]: deleting_dihedrals[1].sort()
      # deleting_dihedrals now contains all of our dihedral indexes
      if deleting_dihedrals[0]:
         while deleting_dihedrals[0]:
            del self.parm.dihedrals_inc_h[deleting_dihedrals[0].pop()]
      if deleting_dihedrals[1]:
         while deleting_dihedrals[1]:
            del self.parm.dihedrals_without_h[deleting_dihedrals[1].pop()]
