#!/usr/bin/env python

import cPickle,sys
import lib.Etc as Etc
from lib.Atom import Atom

def residMap(pdb,outputName,chk):
    """
    usage: residMap(pdb,outputName,chk)
    standalone usage: $ residMap.py <pdb> <outputName> <chk>

    Takes a pdb file, and based upon the infomation from the chk file, builds a map from
    the charmming resid to the canonical resid.  Then prints a pdb with the same coordinates
    as the input, only with canonical resid indexing.
    """
#   (1) Build Map
    pickledPdb = cPickle.load(open(chk))
    pickledText = [ atom.line for atom in pickledPdb ]
    chk_ver = Etc.get_ver_iter(pickledText)
    pickledOldResid = [ Atom(line,chk_ver).resid for line in pickledText ]
    pickledNewResid = [ (atom.type,atom.segid,atom.resid) for atom in pickledPdb ]
    residMap = dict( zip(pickledNewResid,pickledOldResid) )
#   (2) Remap
    pdb_ver = Etc.get_ver(pdb)
    newPdb = [ Atom(line,pdb_ver) for line in open(pdb) if line.startswith('ATOM') or line.startswith('HETATM') ]
    for atom in newPdb:
        atom.resid = residMap[(atom.type,atom.segid,atom.resid)]
#   (3) Print
    write_to = open(outputName,'w')
    for atom in newPdb:
        write_to.write(atom.Print(pdb_ver))
    write_to.write('TER\n')
    write_to.close()

if __name__ == '__main__':
    residMap(sys.argv[1],sys.argv[2],sys.argv[3])




