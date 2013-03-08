"""
This is a module that contains functions responsible for parsing the
input file for MMPBSA.py. It must be included with MMPBSA.py to
ensure proper functioning.

Last updated: 04/17/2011

                           GPL LICENSE INFO

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
"""
from MMPBSA_mods.exceptions import InputError

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Variable(object):
   """ Base variable class. It has a name and a single value """

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self, varname, dat_type='int', default=None, chars_to_match=4,
                case_sensitive=False):
      """ Initializes the variable type. Sets the default value as well as
          specifying how many characters are required by the parser to trigger
          recognition
      """
      # Catch illegalities
      if not dat_type in ('int', 'str', 'float'):
         raise InputError('Variable has unknown data type %s' % dat_type)
      
      # You can't match more characters than you have characters!
      chars_to_match = min(chars_to_match, len(varname))

      self.name = varname
      self.datatype = dat_type
      if default != None:
         if self.datatype == 'str': 
            self.value = default.replace("'",'').replace('"','')
         else:
            self.value = eval('%s("%s")' % (self.datatype, default))
      else: self.value = None
      self.tomatch = chars_to_match
      self.case_sensitive = case_sensitive

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __str__(self):
      """ Prints statistics for the variable """
      string  = 'Variable name:  %s\n' % self.name
      string += 'Variable type:  %s\n' % self.datatype
      string += 'Variable value: %s\n' % self.value
      string += 'Matching chars: %s\n' % self.tomatch
      return string

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __eq__(self, teststring):
      """ Determines if a variable string matches this variable """

      if len(teststring) > len(self.name):
         return False

      if len(teststring) < self.tomatch:
         return False

      myvar = self.name
      if not self.case_sensitive:
         myvar = self.name.lower()
         teststring = teststring.lower()

      return myvar[:len(teststring)] == teststring

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __ne__(self, teststring):
      """ Not equal """
      return not self.__eq__(teststring)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def SetValue(self, value):
      """ Sets the value of the variable """
      if self.datatype != 'str':
         self.value = eval('%s("%s")' % (self.datatype, value))
      else:
         self.value = value.replace('"','').replace("'",'')

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class Namelist(object):
   """ Sets up a namelist. This holds many different Variables, and these
       variables will only be recognized when parsing this particular namelist.
       Set up to mimic the behavior of a Fortran namelist (so that the input is
       similar to the rest of Amber). Some of the known deficiencies: 
       
         o the body of the namelist cannot start on the same line as the start 
           or end of the namelist

         o the end of the namelist must be &end or / and must be on its own line

         o It will not (yet) recognize array lengths -- those are simply parsed
           as strings
   """

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self, trigger, full_name, to_match=3):
      """ Sets up the list of variables, which is just the trigger for now. The
          trigger is a logical variable that gets set to true if this namelist
          is parsed. Open() trips the trigger if it exists. It can be passed in
          as anything that evaluates to False if it doesn't exist.
      """
      self.trigger = trigger
      self.variables = {}
      if self.trigger: self.variables = {self.trigger : False}
      self.open = False
      self.full_name = full_name
      self.to_match = to_match

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __eq__(self, nml):
      """ A namelist is equal if the name matches properly """
      return nml == self.full_name[:len(nml)] and len(nml) >= \
                    min(self.to_match, len(self.full_name))

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __ne__(self, nml):
      """ Not equal """
      return not self.__eq__(nml)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def addVariable(self, varname, datatype, default=None, to_match=4):
      """ Adds a variable to this namelist. It checks to make sure that it's
          going to create a conflict with an existing variable.
      """

      if varname in self.variables.keys():
         raise InputError('Duplicated variable %s in Namelist' % varname)

      for var in self.variables.keys():

         if var == varname:
            raise InputError(('New variable %s triggered a match with existing '
                            + 'variable %s') % (varname, var))

      self.variables[varname] = Variable(varname, datatype, default, to_match)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Open(self):
      """ Signifies that the namelist is open """
      if self.open:
         raise InputError('Namelist already open. Cannot open before closing')

      if self.trigger: self.variables[self.trigger] = True
      self.open = True

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

class InputFile(object):
   """ Defines the Input File and parses it. You have to add stuff to the parser
       via addNamelist. Use it as follows:

       input = InputFile()

       input.addNamelist('gb', 'gb', [['saltcon', 'float', 0, 4], ...], 
                         trigger='gbrun')
       input.addNamelist('ala', 'alanine_scanning', [['mutant', 'int', 0, 4]],
                         trigger='alarun')

       INPUT = input.Parse('mmpbsa.in')
   """

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __init__(self):
      """ Initializes the input file, sets up empty arrays/dictionaries """
      self.ordered_namelist_keys = []
      self.namelists = {}
      self.text = '' # text of the given input file

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def __str__(self):
      """ Prints out the input file """
      if not self.text:
         return

      self.text = self.text.replace('\n', '\n|') # Add | to start of each line

      return ('|Input file:\n|---------------------------------------' + 
              '-----------------------\n|' + self.text + 
              '-----------------------------------------------------' + 
              '---------\n')

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def addNamelist(self, name, full_name, variable_list, trigger=None):
      """ Adds a namelist to the input file that will be parsed. Variable list
          must be an array of arrays. Each array element must be an array that
          has the information [varname, datatype, default, chars to match]. The
          'trigger' is the logical variable that gets set to true if this
          namelist is specified in the input file.
      """

      if name in self.ordered_namelist_keys:
         raise InputError('Namelist %s defined multiple times' % name)

      self.ordered_namelist_keys.append(name)
      self.namelists[name] = Namelist(trigger, full_name)

      for var in variable_list:

         if type(var).__name__ != 'list' and len(var) != 4:
            raise InputError('variables in variable_list must be lists of ' +
                             'length 4. [varname, datatype, chars to match, ' +
                             'default]')

         self.namelists[name].addVariable(var[0], var[1], var[2], var[3])

      # end for var

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def _full_namelist_name(self, nml):
      """ Determines what the full namelist name is. We try to make as many
          allowances as possible. We will match the first 3 characters and
          replace all _'s with 
      """
      nml = nml.replace(' ', '_') # replaces spaces with _'s
      for key in self.ordered_namelist_keys:
         if self.namelists[key] == nml: return key

      raise InputError('Unrecognized namelist %s' % nml)

   #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

   def Parse(self, filename):
      """ This subroutine parses the input file. Only data in namelists are 
          parsed, and all namelists must be set prior to calling this routine.

          It will create a dictionary of Input variables for all variables in
          all namelists. They all flood the same namespace. If there are any
          conflicts between variables in namelists, an error will be raised.
          Make sure all input variables are unique!
      """
      from os.path import exists

      # Make sure our file exists

      if not exists(filename):
         raise InputError("Can't find input file (%s)" % filename)

      # Load the whole thing into memory. This should be plenty short enough.
      infile = open(filename, 'r')
      lines = infile.readlines()
      infile.close()

      for line in lines: # save the input file so we can echo it back later
         self.text += line

      # We will loop through the input file three times: 
      # 
      # 1st: Load all of the data into an array (namelist_fields)
      # 2nd: Combine multi-element values (to allow commas in input variables)
      # 3rd: Loop through the fields to change the values of the variables.

      declared_namelists = [] # keep track of the namelists we've found so far
      namelist_fields = []    # entries in a given namelist
      innml = False           # are we in a namelist now? Don't enter multiple
       
      # split up the input file into separate fields by comma
      
      for line in lines:
         # Skip title lines (we are flexible here) and comments
         if not innml and not line.strip().startswith('&'): 
            continue
         if line.strip().startswith('#') or line.strip().startswith('!'): 
            continue

         # Catch some errors
         if innml and line.strip().startswith('&'):
            raise InputError('Invalid input. Terminate each namelist prior ' +
                             'to starting another one.')
         
         # End of a namelist
         elif innml and line.strip() in ['/', '&end']:
            innml = False

         # Now if we finally find a namelist
         elif not innml and line.strip().startswith('&'):
            innml = True
            namelist = line.strip()[1:].lower()
            namelist = self._full_namelist_name(namelist)

            if namelist in declared_namelists:
               raise InputError('Namelist %s specified multiple times' % 
                                namelist)
            
            self.namelists[namelist].Open()
            declared_namelists.append(namelist)
            namelist_fields.append([])

         # We are in a namelist here, now fill in the fields
         elif innml:
            items = line.strip().split(',')

            # Screen any blank fields
            j = 0
            while j < len(items):
               items[j] = items[j].strip()
               if len(items[j]) == 0:
                  items.pop(j)
               else:
                  j += 1

            namelist_fields[len(namelist_fields)-1].extend(items)

         # end if [elif innml]

      # end for line in lines

      # Combine any multi-element fields into the last field that has a = in it

      begin_field = -1
      for i in range(len(namelist_fields)):
         for j in range(len(namelist_fields[i])):
            if not '=' in namelist_fields[i][j]:
               if begin_field == -1:
                  raise InputError('Invalid input file! Error reading ' +
                                   'namelist %s' % declared_namelists[i])
               else:
                  namelist_fields[i][begin_field] += \
                                       ',%s' % namelist_fields[i][j]
            else:
               begin_field = j

      # Now parse through the items to add them to the master dictionary. Note
      # that thanks to the last step, all data in namelist_fields will be 
      # contained within fields that have a '='. All others can be ignored

      for i in range(len(namelist_fields)):
         for j in range(len(namelist_fields[i])):

            if not '=' in namelist_fields[i][j]:
               continue
            else:
               var = namelist_fields[i][j].split('=')
               var[0] = var[0].strip()
               var[1] = var[1].strip()

               # Now we have to loop through all variables in that namelist to
               # see if this is the variable we want.

               found = False
               for key in self.namelists[declared_namelists[i]].\
                           variables.keys():
                  if self.namelists[declared_namelists[i]].variables[key] == \
                            var[0]:
                     self.namelists[declared_namelists[i]].variables[key].\
                                                  SetValue(var[1])
                     found = True
                     break

               if not found: 
                  raise InputError('Unknown variable %s in &%s' % (var[0],
                                   declared_namelists[i]))

      # Now it's time to fill the INPUT dictionary

      INPUT = {}

      for nml in self.ordered_namelist_keys:
         for var in self.namelists[nml].variables.keys():

            # Here, the triggers are just bool types, so protect from accessing
            # an attribute that doesn't exist! We only allow Variable types and
            # bool types

            var_object = self.namelists[nml].variables[var]
            try:
               INPUT[var] = self.namelists[nml].variables[var].value
            except AttributeError:
               if type(var_object).__name__ == 'bool':
                  INPUT[var] = var_object
               else:
                  raise InputError('Disallowed namelist variable type')

      return INPUT

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
