#!/usr/bin/python

import sys

if len(sys.argv) < 5:
   print "Insufficient arguments.\n"
   sys.exit(0)

ofp = open(sys.argv[4], 'w+')
ofp.write("* sets greatervalue\n")
ofp.write("*\n\n")
ofp.write("set greatervalue = %s\n" % sys.argv[1])
ofp.write("set secondvalue = %s\n" % sys.argv[2])
ofp.write("set xtaltype = %s\nreturn\n" % sys.argv[3])
ofp.close()
