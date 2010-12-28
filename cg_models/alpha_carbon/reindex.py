#!/usr/bin/env python

from lib.Atom import Atom
import sys

atomList = [ Atom(line,'charmm') for line in open(sys.argv[1]) if line.startswith('ATOM') ]

old_resid = atomList[0].resid
new_resid = 1
for i,atom in enumerate(atomList):
    atom.tag   = 'ATOM'
    atom.segid = 'A'
    atom.atomNumber = i+1
    if atom.resid == old_resid:
        atom.resid = new_resid
    else:
        old_resid = atom.resid
        new_resid += 1
        atom.resid = new_resid

write_to = open(sys.argv[1],'w')
for atom in atomList:
    write_to.write(atom.Print('charmm'))
write_to.write('TER\n')
write_to.close()


