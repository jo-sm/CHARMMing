""" Useful utilities for Parmed """

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def LineToCmd(line):
   """ Translates a parmed user command into an action """
   tokens = line.strip().split()
   cmd = 'ParmedActions.%s(amber_prmtop' % tokens[0].lower()
   for i in range(len(tokens)-1):
      cmd += ",'%s'" % tokens[1+i]
   return cmd + ')'

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def PrintHelp(usage, myhelp):
   """ Prints the help message in a nicely formatted way, for 80-char terms """
   def addOn(line, addition):
      if len(line) + len(addition) > 80:
         print line
         line = '    %s' % addition
      elif len(line) > 0:
         line += ' %s' % addition
      elif len(line) == 0:
         line = '   %s' % addition
      return line

   print ''
   print usage.strip()
   words = myhelp.strip().split()
   line = '   '
   for word in words:
      line = addOn(line, word)
   print line + '\n'

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def PrintHelpSummary(help_dict):
   """ Prints a generic help message with all functions in ParmedActions """
   import sys
   print ''
   key_list = help_dict.keys()
   key_list.sort()
   longest_option = 0
   for key in key_list: longest_option = max(len(key), longest_option)
   print '%%-%ds | Usage' % longest_option % 'Command'
   for i in range(80): sys.stdout.write('-')
   print ''
   for key in key_list:
      print '%%-%ds' % longest_option % help_dict[key].split()[0] + ' | ' + \
            help_dict[key]
   print ''

#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
