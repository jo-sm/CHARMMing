#!/usr/bin/env python

import glob,sys,os,time
import lib.Etc as Etc
import lib.future.PDBFile as pdb

inputPath = Etc.expandPath(sys.argv[1])
outputPath = Etc.expandPath(sys.argv[2])
os.chdir(outputPath)
n = 1
for pdbFileName in glob.glob('%s/*.pdb' % inputPath):
    start = time.time()
    n += 1
    print 'parsing %s' % pdbFileName
    try:
        taco = pdb.PDBFile(pdbFileName)
        taco.WritePDB()
    except:
        print 'Fatal Exception Caught >>>', str(sys.exc_info())
    print 'elapsed time: %10.2f s' % (time.time() - start)
#   if n > 1: break

