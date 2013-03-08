#!/usr/bin/env python

# This is an attempt at a python preprocessor that should work like cpp. It will
# only recognize #ifdefs and #ifndefs right now.

import sys, os

sys.stdout = os.fdopen(sys.stdout.fileno(),'w',0)
sys.stderr = os.fdopen(sys.stderr.fileno(),'w',0)

defines = []
defined = []
elsed = False

def usage():
   print >> sys.stderr, 'ppp.py -DDEFINE1 -DDEFINE2 ... file.py'
   sys.exit()

if len(sys.argv) < 2 or sys.argv[1] == '-h' or sys.argv[1] == '--help':
   usage()

for i in range(len(sys.argv)-1):
   if sys.argv[i].startswith('-D'):
      defines.append(sys.argv[i][2:])

source_file = sys.argv[len(sys.argv)-1]

file = open(source_file,'r')
lines = file.readlines()
file.close()

i = 0
while i < len(lines):
   if lines[i].lower().strip().startswith('#ifdef'):
      defined.append(lines[i].strip()[7:].strip())
      i += 1
      continue
   elif lines[i].lower().strip().startswith('#ifndef'):
      defined.append('!%s' % lines[i].strip()[8:].strip())
      i += 1
      continue
   elif lines[i].lower().strip().startswith("#else"): # fix this so it will do nested elses
      if elsed:
         print >> sys.stderr, 'Error: Two #else constructs found without terminating the first!'
         sys.exit()
      elsed = True
      if defined[len(defined)-1].startswith('!'):
         defined[len(defined)-1] = defined[len(defined)-1][1:]
      else:
         defined[len(defined)-1] = '!' + defined[len(defined)-1]
      i += 1
      continue
   if lines[i].lower().strip().startswith('#endif'):
      elsed = False
      defined.pop()
      i += 1
      continue

   toPrint = True

   for j in range(len(defined)):
      if '!' in defined[j] and defined[j][1:] in defines:
         toPrint = False
         break
      if not '!' in defined[j] and not defined[j] in defines:
         toPrint = False
         break

   if toPrint:
      sys.stdout.write(lines[i])
   i += 1
