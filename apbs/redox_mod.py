#
#                            PUBLIC DOMAIN NOTICE
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the authors' official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software is freely available
#  to the public for use.  There is no restriction on its use or
#  reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, NIH and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. NIH, NHLBI, and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any
#  particular purpose.

"""
PDB parser for [4Fe-4S] proteins

:Author: bsp
:Date: 07/07/2011

:Usage:
        ``setup_4fsr.py --help`` will give you a help message explaining the various
        options.  "Required options" are marked as such, and defaults
        appear in [brackets].

:Options:
    | ``-h, --help``  show this help message and exit
    | ``-I PATH, --input=PATH``
        *Required* PATH of input .pdb file
    | ``-o DIR``, ``--redox_site_o=4FSO``
        Residue name of oxidized [4Fe-4S] redox site [4fso]
    | ``-r DIR``, ``--redox_site_r=4FSR``
        Residue name of reduced [4Fe-4S] redox site [4fsr]
"""


from pychm.io.pdb import PDBFile
from pychm.io.rtf import RTFFile
import pychm

def fesSetup(pdbFilename, clusnameo, clusnamer, **kwargs):
    """
    Parse a *.pdb* plain text file into its constituent chains and segments with consideration
    of iron-sulfur redox sites. Print one CHARMM formatted *.pdb* file per chain/segment 
    combination.

    *kwarg defaults are listed first*

    **kwargs:**
    """
    # Open PDBFile, parse, inte_seg, etc.
    pdbFileArgs = {
                }
    pdb = PDBFile(pdbFilename, **pdbFileArgs)
    thisMol = pdb.iter_models().next()
    thisMol.parse()

    rtf = RTFFile('top_all22_4fe4s_esp_090209.dat')
    #######################################
    # Files split into mol objects        #
    #######################################

    # Select redox site
    fes4 = ['4fso','4fsr','4fss']
    if clusnameo in fes4:
        numFe = 4 # Number of Irons
        numS = 4  # Number of Sulfurs
    else:
        print 'Redox site is not recognized'

    # Collect redox site resnames
    het = list(thisMol.iter_res(segtypes = ['bad'],resName = ['fs4','sf4']))
    #het = list(thisMol.find(segtypes=['bad']))

    renameFeS(thisMol,het,clusnameo,numS,numFe)
    renameLigands(pdb,thisMol,het,clusnameo,numFe)
    chargeFeS(rtf,thisMol,clusnameo) # NEW #

    segDict = {'nuc':'nuc', 'pro':'pro', 'good':'goodhet', 'bad':'het',
                 'dna':'dna', 'rna':'rna'}
    stdoutList = []
    # Write out structure
    for seg in thisMol.iter_seg():
        stdoutList.append('%s-%s' % (seg.chainid, segDict[seg.segType]))
        name = 'new_%s-%s-%s_o.pdb' % (thisMol.code, seg.chainid,segDict[seg.segType])
        seg.write(filename = name)

    chargeFeS(rtf,thisMol,clusnamer,) # NEW #

    stdoutList = []
    # Write out structure
    for seg in thisMol.iter_seg():
        stdoutList.append('%s-%s' % (seg.chainid, segDict[seg.segType]))
        name = 'new_%s-%s-%s_r.pdb' % (thisMol.code, seg.chainid,segDict[seg.segType])
        seg.write(filename = name)


def renameFeS(thisMol,het,clusname,numS,numFe):
    atype = ['fe1 ','fe2 ','fe3 ','fe4 ',' s1 ',' s2 ',' s3 ',' s4 ']
    totinorg = numFe+numS
    sufdist0 = [0]*numFe
    # Rename atoms, Collect distances
    for res in het:
        fe1dist = [0]*numS
        new = [0]*totinorg
        old = [ atom for atom in res ]
        iron = [ atom for atom in thisMol.find(chainid = res.chainid,segtype = 'bad',resid = res.resid) if atom.element == 'fe' ]
        sulfur = [ atom for atom in thisMol.find(chainid = res.chainid,segtype = 'bad',resid = res.resid) if atom.element == 's' ]

        # Rename [4Fe-4S] cubane atoms - SLOPPY
        for s in range(numS):
            fe1dist[s] = iron[0].calc_length(sulfur[s])
        # FE1 is always FE1
        new[0] = old[0]                                                             # SET FE1
        # S4 is sulfur furthest from FE1
        iS4 = fe1dist.index(max(fe1dist))
        new[7] = sulfur[iS4]                                                        # SET S4
        # S3 is sulfur closest to FE1
        iS3 = fe1dist.index(min(fe1dist))
        new[6] = sulfur[iS3]                                                        # SET S3
        # FE2 is iron furthest from S3

        for f in range(numFe):
            sufdist0[f] = sulfur[iS3].calc_length(iron[f])
        iFe2 = sufdist0.index(max(sufdist0))
        new[1] = iron[iFe2]                                                         # SET FE2
        # FE3 & FE4 are just other irons
        if iFe2 == 1: iFe3 = 2; iFe4 = 3
        if iFe2 == 2: iFe3 = 1; iFe4 = 3
        if iFe2 == 3: iFe3 = 1; iFe4 = 2
        new[2] = iron[iFe3]                                                         # SET FE3
        new[3] = iron[iFe4]                                                         # SET FE4
        # S1 & S2 are just other sulfurs
        if iFe2 == 1: iFe3 = 2; iFe4 = 3
        iS1 = 5
        iS2 = 5
        for s in range(numS):
            if s !=  iS1 and s !=  iS2 and s !=  iS3 and s !=  iS4:
                if iS1 == 5: iS1 = s
                else: iS2 = s
        new[4] = sulfur[iS1]                                                        # SET S1
        new[5] = sulfur[iS2]                                                        # SET S2

        # Rename inorgaic atoms
        for j in range(0,8):
            new[j].atomType = atype[j]
        res.resName = clusname


from copy import deepcopy

def renameLigands(pdb,thisMol,het,clusname,numFe):
    # Make Metal/Ligand table from LINK statements, if present
    # ar[l][0/4] = metal/ligand atomtype
    # ar[l][1/5] = metal/ligand resname
    # ar[l][2/6] = metal/ligand chain id
    # ar[l][3/7] = metal/ligand resid
    try:
        ar = [ line.split() for line in pdb.get_metaData()['link'] ]
    except KeyError:
        print 'No LINK statements in PDB. Tell Scott to make a fix.'
    # For each fe in res in het, select bound residue
    # Select side chain of residue
    # Deep copy side chain to het
    # Rename side chain atoms
    for res in het:
        cysid=[]
        # Get resid for ligated cys
        for fe in res:
            for line in ar:
                if line[0].startswith('fe'):
                    t = 0
                    u = 4
                else:
                    t = 4
                    u = 0
                if fe.atomType.strip() == line[t].strip() and fe.resid0 == int(line[t+3]) and res.chainid == line[t+2]:
                    cysid.append(line[u+3])
        # Deepcopy side chains of listed residues
        s=1
        c=1
        for cys in cysid:
            # TO USE MULTIPLE CHAINS/CLUSTER,NEED TO SET CHAIN ID TO THAT OF THE ATOM
            for atom in thisMol.find(resid=cys,chainid=res.chainid):
                if atom.atomType==' cb ':
                    temp = deepcopy(atom)
                    temp.atomType = ' cb'+str(c)
                    temp.resName = clusname
                    temp.segtype = 'bad'
                    temp.resid = res.resid
                    thisMol.append(temp)
                    c=c+1
                if atom.atomType==' sg ':
                    temp = deepcopy(atom)
                    temp.atomType = ' sg'+str(s)
                    temp.resName = clusname
                    temp.segtype = 'bad'
                    temp.resid = res.resid
                    thisMol.append(temp)
                    s=s+1

from pychm.io.rtf import RTFFile

def chargeFeS(rtf,thisMol,clusname):
    het = list(thisMol.iter_res(segtypes = ['bad'],resName = ['fs4','sf4','4fsr','4fso','4fss']))
    for res in het:
        for atom in res:
            atom.resName = clusname
    # Populate charges in protein
    thisMol.populate_charges(rtf)


if __name__ == '__main__':


    import sys
    from pychm.tools import OptionParser

    # Option Parsing
    useText  = \
    """
    '%prog --help' will give you a help message explaining the various
    options.  "Required options" are marked as such, and defaults
    appear in [brackets].
    """
    optparser = OptionParser(usage = useText, version = '%prog 0.1')
    # Required
    optparser.add_option('-I', '--input', required = True, metavar = 'PATH',
                    help = 'PATH of input .pdb file')
    optparser.add_option('-o', '--redox_site_o', default = '4fso',
                    help = 'Residue name of oxidized redox site')
    optparser.add_option('-r', '--redox_site_r', default = '4fsr',
                    help = 'Residue name of reduced redox site')

    # Parse
    (options, args) = optparser.parse_args(sys.argv)
    # Repackage options into kwargs
    kwargs = {
        }
    # Do Work
    fesSetup(options.input, options.redox_site_o, options.redox_site_r, **kwargs)