#! /opt/epd/epd-6.2-2-rh3-x86_64/bin/python

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
:Date: 11/09/2011

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
    | ``-N Number``, ``--redox_site_id=Number``
        Residue IF of redox site [1]
"""


from pychm.io.pdb import PDBFile
from pychm.io.rtf import RTFFile
import pychm

def fesSetup(thisMol, clusnameo, rtf, clusn, mutid, location, identifier, pdb_metadata, **kwargs):
    """
    Parse a *.pdb* plain text file into its constituent chains and segments with consideration
    of iron-sulfur redox sites. Print one CHARMM formatted *.pdb* file per chain/segment 
    combination.

    *kwarg defaults are listed first*

    **kwargs:**
    """
    # Select redox site
    fes4 = ['4fso','4fsr','4fss']
    if clusnameo in fes4:
        numFe = 4 # Number of Irons
        numS = 4  # Number of Sulfurs
    else:
        raise NotImplementedError('Redox site is not recognized')
    # Collect redox site resnames
    het = list(thisMol.iter_res(segtypes = ['bad'],resName = ['fs4','sf4']))
    #het = list(thisMol.find(segtypes=['bad']))
    renameFeS(thisMol,het,clusnameo,numS,numFe)
    renameLigands(pdb_metadata,thisMol,het,clusnameo,numFe)
    chargeFeS(rtf,thisMol,clusnameo) # NEW #
    segDict = {'nuc':'nuc', 'pro':'pro', 'good':'goodhet', 'bad':'het',
                 'dna':'dna', 'rna':'rna'}
    stdoutList = []
    # Mutate Residue if Exist
    if mutid != 0: turnOffChgs(rtf,thisMol,mutid)
    # Write out structure
    for seg in thisMol.iter_seg():
        stdoutList.append('%s-%s' % (seg.chainid, segDict[seg.segType]))
        name = '%s/redox-%s-%s-%s_o.pdb' % (location, identifier, seg.chainid, segDict[seg.segType])
        seg.write(filename = name)
        #RTFFile.write(filename = 'test.rtf')
    # Select Reduced Residue Name
    if clusnameo == '4fso': clusnamer = '4fsr'
    if clusnameo == '4fsr': clusnamer = '4fss'
    reduceFeS(rtf,thisMol,clusnamer,clusn) # NEW #
    stdoutList = []
    # Mutate Residue if Exist
    if mutid != 0: turnOffChgs(rtf,thisMol,mutid)
    # Write out structure
    for seg in thisMol.iter_seg():
        stdoutList.append('%s-%s' % (seg.chainid, segDict[seg.segType]))
        name = '%s/redox-%s-%s-%s_r.pdb' % (location, identifier, seg.chainid, segDict[seg.segType])
        seg.write(filename = name)


def renameFeS(thisMol,het,clusname,numS,numFe):
    """
    Renames Fe & S atoms of [4Fe-4S] cubane to account for          FE1-------S2
    mixed-valence layers. 'Top' layer consists of atoms FE1,       /|         /|
    S1, FE2, and S2, while the bottom layer consist of atoms      / |        / | 
    FE3, S3, FE4, and S4. First, the sulfur furthest from FE1    S1-+-----FE2  |
    is renamed S4 and the sulfur closets to FE1 is renamed S3    |  S3------+FE4
    The iron furthest from S3 is renamed FE2. Next the two       | /        | /
    remaining irons are renames FE3 and FE4. Similarly, two      |/         |/
    remaining two sulfurs are renamed S1 and S2.                 FE3-------S4
    """
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
        new[0] = old[0]                                      # SET FE1
        # S4 is sulfur furthest from FE1
        iS4 = fe1dist.index(max(fe1dist))
        new[7] = sulfur[iS4]                                 # SET S4
        # S3 is sulfur closest to FE1
        iS3 = fe1dist.index(min(fe1dist))
        new[6] = sulfur[iS3]                                 # SET S3
        # FE2 is iron furthest from S3
        for f in range(numFe):
            sufdist0[f] = sulfur[iS3].calc_length(iron[f])
        iFe2 = sufdist0.index(max(sufdist0))
        new[1] = iron[iFe2]                                  # SET FE2
        # FE3 & FE4 are just other irons
        if iFe2 == 1: iFe3 = 2; iFe4 = 3
        if iFe2 == 2: iFe3 = 1; iFe4 = 3
        if iFe2 == 3: iFe3 = 1; iFe4 = 2
        new[2] = iron[iFe3]                                  # SET FE3
        new[3] = iron[iFe4]                                  # SET FE4
        # S1 & S2 are just other sulfurs
        if iFe2 == 1: iFe3 = 2; iFe4 = 3
        iS1 = 5
        iS2 = 5
        for s in range(numS):
            if s !=  iS1 and s !=  iS2 and s !=  iS3 and s !=  iS4:
                if iS1 == 5: iS1 = s
                else: iS2 = s
        new[4] = sulfur[iS1]                                 # SET S1
        new[5] = sulfur[iS2]                                 # SET S2
        # Rename inorgaic atoms
        for j in range(0,8):
            new[j].atomType = atype[j]
        res.resName = clusname


from copy import deepcopy

def renameLigands(pdb_metadata,thisMol,het,clusname,numFe):
    """
    Cysteine residues ligated to Fe atoms of the redox site are selected
    based on the LINK statement in the PDB. The CB, SG, and their bound
    hydrogens are copied to the redox site and atom names are renamed to
    match the redox site atoms in the CHARMM topology file.

    For [4Fe-4S] redox site:

        HB11 HB12
          \  / 
           CB1-SG1            HB21  HB22
          /      \               \  /
         HB13    FE1--------S2    CB2-HB23
                / |        /|    / Plane A - "mixed" oxd pair
              S1--------FE2----SG2
              |   |       | |
              |  S3-------|-FE4
              | /         |/  \    Plane B - "mixed" oxd pair
              FE3---------S4    SG4-CB2-HB21
              /                     /   \
             SG3                 HB22   HB23
             |
             CB3-HB31      (both planes oxd states delocalized)
            /   \
           HB32 HB33

    """
    # Make Metal/Ligand table from LINK statements, if present
    # ar[l][0/4] = metal/ligand atomtype
    # ar[l][1/5] = metal/ligand resname
    # ar[l][2/6] = metal/ligand chain id
    # ar[l][3/7] = metal/ligand resid
    try:
        ar = [ line.split() for line in pdb_metadata['link'] ]
    except KeyError:
        raise NotImplementedError('No LINK statements in PDB. Tell Scott to make a fix.')
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
    """
    All atoms in thisMol are assigned parial charges. Partial charges for the
    redox site are determined by clusname. 
    """
    het = list(thisMol.iter_res(segtypes = ['bad'],resName = ['fs4','sf4','4fsr','4fso','4fss']))
    for res in het:
        for atom in res:
            atom.resName = clusname
    # Populate charges in protein
    thisMol.populate_charges(rtf)


def reduceFeS(rtf,thisMol,clusname,clusn):
    """
    All atoms in thisMol are assigned partial charges. Partial charges for the
    redox site with resid == clusn are determined by clusname.
    """
    het = list(thisMol.iter_res(segtypes = ['bad'],resName = ['fs4','sf4','4fsr','4fso','4fss']))
    for res in het:
        if res.resid == clusn:
            for atom in res:
                atom.resName = clusname
    # Populate charges in protein
    thisMol.populate_charges(rtf)

def turnOffChgs(rtf,thisMol,mutid):
    """
    All atoms of the residue with resid == mutid are given a partial charge of zero.
    """
    pro = list(thisMol.iter_res(segtypes = ['pro']))
    for res in pro:
        if res.resid == mutid:
            for atom in res:
                atom.charge = 0.000

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
    optparser.add_option('-N', '--redox_site_id', default = '1',
                    help = 'Residue ID of redox site')
    optparser.add_option('-M', '--mutID', default = '0',
                    help = 'ID for residue ID')
    # Parse arguments
    (options, args) = optparser.parse_args(sys.argv)
    # Repackage options into kwargs
    kwargs = {
        }
    # Convert PDB to mol object
    pdbFileArgs = {
                }
    pdbFilename=options.input
    pdb = PDBFile(pdbFilename, **pdbFileArgs)
    thisMol = pdb.iter_models().next()
    thisMol.parse()
    # RTF File
    rtf = RTFFile('top_all22_4fe4s_esp_090209.dat')
    options.redox_site_id=1
    mutid=0
    # Do Work
    fesSetup(thisMol, options.redox_site_o, rtf, options.redox_site_id, mutid, **kwargs)
