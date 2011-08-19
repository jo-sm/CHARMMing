#!/usr/bin/python

import sys

if len(sys.argv) < 3:
   print "Insufficient arguments."
   sys.exit(0)

fp = open(sys.argv[2], "w+")
mycharge = round(float(sys.argv[1]))
fp.write("CHARGE=%d\n" % mycharge)
fp.close()
