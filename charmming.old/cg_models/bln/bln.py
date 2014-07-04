#!/usr/bin/env python
# Filename bln.py
"""
# Outline
# ( 0a) Option Parsing
# ( 0b) Option Processing
# ( 1 ) Initialize All-Atom
# ( 2a) Run stride to get H-bonds
# ( 2b) Load contact potential
# ( 2c) Load mapping of residues -> BLN types
# ( 3 ) Define more efficient Rij machinery
# ( 4 ) Reindex CG.atomNumber
# ( 5 ) Write .crd file
# ( 6 ) Write a .seq file
# ( 7 ) Write .pdb file
# ( 8 ) Write a .prm file
"""

# BLN protein model. This code is heavily based on the code implementing
# the Klimov-Thirumulai Go-like model in CHARMMing. The BLN model interface
# was implemented by Frank Pickard IV and Benjamin T. Miller, NHLBI,
# NIH, March, 2010.

# when using this code, please cite the following references:
#
# SJ Marrink, AH de Vries, and AE Mark. J. Phys. Chem B., vol. 108, pp. 750-750. 2004.
# PJ Bond and MSP Sansom. J. Am. Chem. Soc., vol. 128, pp. 2697-2704. 2006.

import sys,commands,math,numpy,os,cPickle,string
import lib.Etc as Etc
from lib.Atom import Atom
from lib.Res import Pro
from lib.OptParse import OptionParser
from optparse import OptionGroup

# ( 0a) Option Parsing
parser = OptionParser(usage="""'%prog --help' will give you a help message with lots of useful information
        "Required options" are marked as such, and defaults appear in [brackets].""",version='%prog 0.1')
#debug
parser.add_option('--debug',
        action='store_true',default=False,
        help='turn on (lots of) debuffing information [%default]')
#input:path+name
parser.add_option('-I','--input',required=True,metavar='PATH',
        help='PATH of input .pdb file')
#output:path
parser.add_option('-O','--output',default=None,metavar='DIR',
        help='DIR where CHARMM files are written [$inputPath]')
#data:path
parser.add_option('-D','--data',default='./', metavar='DIR',
        help='DIR where contact parameter data files are stored [%default]')
#input:formatting
parser.add_option('--formatting',default='auto',choices=['auto','charmm','web'],metavar='FORM',
        help='specify the .pdb formatting to expect: ([auto],charmm,web)')
#STRIDE_path:path
parser.add_option('--stride_path',default='./',metavar='PATH',
        help='PATH to stride binary file, which is used for determining secondary structure elements [%default]')
#parm_path -- parameters.txt
parser.add_option('--parm_path',default=None,metavar='PATH',
        help='PATH to a .txt file which defines parameters to be used in building the CHARMM scripts.')
#pdb_only
parser.add_option('--pdb_only',
        action='store_true',default=False,
        help='only write a CG pdb file, without the .rtf and .prm files required to utilize it')
#charge_stream_path
parser.add_option('--charge_stream_path',default='charges.str',metavar='PATH',
        help='PATH to write a stream file containing the system charges')
################
###Parameters###
################
group = OptionGroup(parser,'Parameter Options',
        "All of the following options allow for (command line) fine tweaking of (almost) all of the parameters utilized by this script.")
   #nScale
group.add_option('--nScale',default=1.00,type='float',metavar='N',
        help='A user specified arbitrary energy scaling constant. [%default]')
   #domainScale
group.add_option('--domainScale',default=1.,type='float',metavar='N',
        help='A user specified arbitrary energy scaling factor for interactions between residues in different domains. [%default]')
   #kBond
group.add_option('--kBondHelix',default=2.988,type='float',metavar='K',
        help='Specify the spring constant for bonds. [%default]')
group.add_option('--kBondSheet',default=2.988,type='float',metavar='K',
        help='Specify the spring constant for bonds. [%default]')
group.add_option('--kBondCoil',default=2.988,type='float',metavar='K',
        help='Specify the spring constant for bonds. [%default]')
   #kAngle
group.add_option('--kAngleHelix',default=8.37,type='float',metavar='K',
        help='Specify the spring constant for angles. [%default]')
group.add_option('--kAngleSheet',default=8.37,type='float',metavar='K',
        help='Specify the spring constant for angles. [%default]')
group.add_option('--kAngleCoil',default=5.98,type='float',metavar='K',
        help='Specify the spring constant for angles. [%default]')
# Load Options
parser.add_option_group(group)
(options,args) = parser.parse_args(sys.argv)

# ( 0b) Option Processing

# debug
debug = options.debug
# input
inputPath = Etc.expandPath(options.input)
if os.path.isfile(inputPath):
    inputFile = os.path.basename(inputPath)
    if inputFile.endswith('.pdb'):
        pdbName = inputFile.split('.pdb')[0]
    else:
        raise IOError('input: inputFile %s does not appear to be a .pdb file'%inputFile)
    inputPath = os.path.dirname(inputPath) + os.sep
else:
    raise IOError('input: inputPath %s does not exist'%inputhPath)
# output
if options.output:
    outputPath = Etc.expandPath(options.output) + os.sep
else:
    outputPath = inputPath          # The Default
# directory with contact set data
contactDir = Etc.expandPath(options.data)
# formatting
if options.formatting == 'auto':
    formatting = Etc.get_ver(inputPath + inputFile)
else:
    formatting = options.formatting
# stride
stridePath = Etc.expandPath(options.stride_path)
if os.path.isfile(stridePath):
    pass
elif os.path.isdir(stridePath):
    stridePath = stridePath + os.sep + 'stride'
    if not os.path.isfile(stridePath):
        raise IOError("stride: can't find stride binary at specified location\n\t%s\n"%stridePath)
else:
    raise IOError("stride: can't find stride binary at specified location\n\t%s\n"%stridePath)

############
#PARAMETERS#
############
nScale              = options.nScale            #  0.05     # A user specified arbitrary energy scaling constant
domainScale         = options.domainScale       #  1.       # A user specified inter domain scaling factor
kBondHelix          = options.kBondHelix        #  3.5      # Spring constant for bonds
kBondSheet          = options.kBondSheet        #  3.5      # Spring constant for sheets
kBondCoil           = options.kBondCoil         #  2.5      # Spring constant for coils
kAngleHelix         = options.kAngleHelix       # Spring constant for angle bends
kAngleSheet         = options.kAngleSheet       # Spring constant for angle bends
kAngleCoil          = options.kAngleCoil        # Spring constant for angle bends

#######################
#Contact Parameter Set#
#######################
# BTMiller, for now we just assume that the user wants the bln_mvm contact file name
# ToDo, allow rescaling of arbitrary all-atom contact potentials into BLN potentials
contactFileName = "%s/bln_mvm.dat" % contactDir
if not os.path.isfile(contactFileName):
    raise IOError("contactFileName: file doesn't exist at the specified path %s"%contactFileName)

# # # # #
# Main  #
# # # # #

# ( 1 ) Initialize All-Atom
#   (a) Atoms
atomList = ( Atom(line,formatting) for line in open(inputPath+inputFile) if line.startswith ('ATOM') )
atomList = [ atom for atom in atomList if atom.element != 'H' ]
#   (b) Residues
resList = []
resDict = {}
buffer = []
for i,atom in enumerate(atomList):
    try:
        buffer.append(atom)
        if atomList[i].resid == atomList[i+1].resid:
            pass
        else:
            resName = atom.segid[0] + str(atom.resid)
            resList.append(resName)
            resDict[resName] = Pro(buffer)
            buffer = []
    except IndexError:
        resName = atom.segid[0] + str(atom.resid)
        resList.append(resName)
        resDict[resName] = Pro(buffer)
#   (c) Backbone / Sidechain
bbDict = {}
scDict = {}
for res in resList:
    bbDict[res] = [ atom for atom in resDict[res] if Etc.is_backbone(atom.atomType) ]
    if resDict[res].resName != 'GLY':
        scDict[res] = [ atom for atom in resDict[res] if not Etc.is_backbone(atom.atomType) ]

# ( 2 ) Load contact potential
#   (a) Parse Specified Contact File
resNameMap      = {}
contactMatrix   = []
for line in Etc.stripBare(open(contactFileName)):
    if line.startswith('BLN'):
        buffer = line.split(':')[1].split()
        resNameMap = dict( zip(buffer,range(len(buffer))) )
    else:
        scaledContacts = [ nScale * abs(float(contact)) for contact in line.split() ]
        contactMatrix.append(scaledContacts)
#   (b) Calculate average interaction
avgEpsilon = sum(Etc.flatten(contactMatrix))/len(Etc.flatten(contactMatrix))
#   (c) Create ResName -> ContactEnergy mapper

def getMainTypesFromHybrid(hybridClass):
    if hybridClass == 'CND':
        return ['C','ND']
    elif hybridClass == 'CBA':
        return ['C','NDA']
    elif hybridClass == 'CQA':
        return ['C','QD']
    else:
        if debug:
            print "Aieee -- %s is NOT a hybrid atom type!" % hybridClass
        raise deadlyErrorOfDeath

def get_contact_energy(atomClass1,atomClass2):
    if atomClass1 in resNameMap and atomClass2 in resNameMap:
        res1 = resNameMap[atomClass1]
        res2 = resNameMap[atomClass2]
        try:
            return float(contactMatrix[res1][res2])
        except IndexError:
            return float(contactMatrix[res2][res1])
    else:
        if atomClass1 in resNameMap:
            ac2t1, ac2t2 = getMainTypesFromHybrid(atomClass2)
            res1   = resNameMap[atomClass1]
            res2t1 = resNameMap[ac2t1]
            res2t2 = resNameMap[ac2t2]
            try:
                r = float(contactMatrix[res1][res2t1])
            except IndexError:
                r = float(contactMatrix[res2t1][res1])
            try:
                r += float(contactMatrix[res1][res2t2])
            except IndexError:
                r += float(contactMatrix[res2t2][res1])
            return r/2.0

        elif atomClass2 in resNameMap:
            ac1t1, ac1t2 = getMainTypesFromHybrid(atomClass1)
            res2   = resNameMap[atomClass2]
            res1t1 = resNameMap[ac1t1]
            res1t2 = resNameMap[ac1t2]
            try:
                r = float(contactMatrix[res1t1][res2])
            except IndexError:
                r = float(contactMatrix[res2][res1t1])
            try:
                r += float(contactMatrix[res1t2][res2])
            except IndexError:
                r += float(contactMatrix[res2][res1t2])
            return r/2.0

        else:
            ac1t1, ac1t2 = getMainTypesFromHybrid(atomClass1)
            ac2t1, ac2t2 = getMainTypesFromHybrid(atomClass2)
            res1t1 = resNameMap[ac1t1]
            res1t2 = resNameMap[ac1t2]
            res2t1 = resNameMap[ac2t1]
            res2t2 = resNameMap[ac2t2]

            try:
                r = float(contactMatrix[res1t1][res2t1])
            except IndexError:
                r = float(contactMatrix[res2t1][res1t1])
            try:
                r += float(contactMatrix[res1t1][res2t2])
            except IndexError:
                r += float(contactMatrix[res2t2][res1t1])

            try:
                r += float(contactMatrix[res1t2][res2t1])
            except IndexError:
                r += float(contactMatrix[res2t1][res1t2])
            try:
                r += float(contactMatrix[res1t2][res2t2])
            except IndexError:
                r += float(contactMatrix[res2t2][res1t2])

            return r/4.0


# ( 2a ) STRIDE 
strideOutput = commands.getstatusoutput('%s -h %s' % (stridePath,inputPath+inputFile) )[1].split('\n')
if strideOutput[0].startswith('Error'): raise IOError('stride: %s' % strideOutput[0])
if debug:
    tfp = open('stride.out', 'w')
    for line in strideOutput:
        tfp.write('%s\n' % line)
    tfp.close()

ssList = []
strideHelixOutput = [ line for line in strideOutput if line.startswith('ASG') ]
for i,line in enumerate(strideHelixOutput):
    if line[34:39] == 'Helix':
        ssList.append('helix')
    elif line[33:39] == 'Strand':
        ssList.append('sheet')
    else:
        ssList.append('other')

strideHBondOutput = [ line for line in strideOutput if line[:3] in ['ACC','DNR'] and int(line[15:20]) < int(line[35:40]) ]
#       (i) Multiplicity, ie anti-parallel beta sheets will have 'double' hbonds
hBondMult = {}
for line in strideHBondOutput:
    if debug:
        print "STRIDE> %s" % line.strip()
    res_i = int(line[15:20])
    res_j = int(line[35:40])
    if (res_i,res_j) in hBondMult:
        hBondMult[(res_i,res_j)] += 1
    else:
        hBondMult[(res_i,res_j)] = 1

if debug:
    print "resList = %s" % resList
    print "hBondMult = %s" % hBondMult
    print "hBondMult.keys = %s" % hBondMult.keys()

# ( 3 ) Map All-Atom into CG

csfp = open(options.charge_stream_path, 'w')
csfp.write("* charges\n*\n\n")

cgbbDict = {}
cgscDict = {}
cgList = []
blnResList = []
for res in resList:
    if debug:
        print "processing residue %s with resid %d" % (res,resDict[res].resid)

    cgbbDict[res] = resDict[res].get_alpha_carbon()
    # FIXME FIXME BTMILLER -- we need to determine if there are Hydrogen
    # bonds in the backbone, and if so, change this to ND, NA, or NDA
    # accordingly
    cgbbDict[res].atomType = 'B' # all backbones are mixed polar-apolar
    cgbbDict[res].mass = 72.0       # suggested Marrink value

    if resDict[res].resName == 'ALA' or resDict[res].resName == 'ILE' or resDict[res].resName == 'LEU' \
       or resDict[res].resName == 'PRO' or resDict[res].resName == 'VAL' or resDict[res].resName == 'PHE':
        cgscDict[res] = resDict[res].get_sc()
        cgscDict[res].atomType = 'S'
        cgscDict[res].mass     = 72.0
        cgscDict[res].resName  = 'BA' # BAP = Hydrophobic, Apolar
        cgbbDict[res].resName  = 'BA'
        blnResName = 'BA'

    elif resDict[res].resName == 'CYS' or resDict[res].resName == 'MET':
        cgscDict[res] = resDict[res].get_sc()
        cgscDict[res].atomType = 'S'
        cgscDict[res].mass     = 72.0
        cgscDict[res].resName  = 'B0' # B0 = Hydrophobic, no hydrogen bonding
        cgbbDict[res].resName  = 'B0'
        blnResName = 'B0'

    elif resDict[res].resName == 'ASN' or resDict[res].resName == 'GLN':
        cgscDict[res] = resDict[res].get_sc()
        cgscDict[res].atomType = 'S'
        cgscDict[res].mass     = 72.0
        cgscDict[res].resName  = 'BB' # BDA = Hydrophobic, hbond donor & acceptor
        cgbbDict[res].resName  = 'BB'
        blnResName = 'BB'

    elif resDict[res].resName == 'SER' or resDict[res].resName == 'THR':
        cgscDict[res] = resDict[res].get_sc()
        cgscDict[res].atomType = 'S'
        cgscDict[res].mass     = 72.0
        cgscDict[res].resName  = 'PL' # POL = polar
        cgbbDict[res].resName  = 'PL'
        blnResName = 'PL'

    elif resDict[res].resName == 'ASP' or resDict[res].resName == 'GLU':
        cgscDict[res] = resDict[res].get_sc()
        cgscDict[res].atomType = 'S'
        cgscDict[res].mass     = 72.0
        cgscDict[res].resName  = 'CA' # CHA = charged acceptor
        cgbbDict[res].resName  = 'CA'
        blnResName = 'CA'
        csfp.write('scalar charge set -1.0 select resid %s .and. type QA end\n' % (resDict[res].resid) )

    elif resDict[res].resName == 'GLY':
        cgbbDict[res].resName = 'BO'  # BO = backbone only, not Barack Obama ;-)
        blnResName = 'BO'

    # the following residue types are "mixed", i.e. in the original model they would
    # have two sidechain beads, here we give them one SC bead that's the average of
    # the two that they would normally receive
    elif resDict[res].resName == 'TYR':
        # Tyrosine gets C + Nd
        cgscDict[res] = resDict[res].get_sc()
        cgscDict[res].atomType = 'S'
        cgscDict[res].mass     = 72.0
        cgscDict[res].resName  = 'CN' # C + Nd
        cgbbDict[res].resName  = 'CN'
        blnResName = 'CN'

    elif resDict[res].resName == 'HSD' or resDict[res].resName == 'TRP':
        # These get C + Nda
        cgscDict[res] = resDict[res].get_sc()
        cgscDict[res].atomType = 'S'
        cgscDict[res].mass     = 72.0
        cgscDict[res].resName  = 'CB' # C + Nda
        cgbbDict[res].resName  = 'CB'
        blnResName = 'CB'

    elif resDict[res].resName == 'LYS' or resDict[res].resName == 'ARG':
        # These get C + Qd
        cgscDict[res] = resDict[res].get_sc()
        cgscDict[res].mass     = 72.0
        cgscDict[res].atomType = 'S'
        cgscDict[res].resName  = 'CQ' # C + Qd
        cgbbDict[res].resName  = 'CQ'
        blnResName = 'CQ'
        csfp.write('scalar charge set 1.0 select resid %s .and. type CQD end\n' % (resDict[res].resid) )


    # now we need to check if there's a backbone hydrogen bond that could
    # change the backbone atom or residue type.
    donor = False
    acceptor = False
    for hBondPair in hBondMult.keys():
        if hBondPair[0] == resDict[res].resid:
            donor = True
            if debug:
                print "This residue is an hbond donor."
        if hBondPair[1] == resDict[res].resid:
            acceptor = True
            if debug:
                print "This residue is an hbond acceptor."

    if donor and acceptor:
        blnResName += 'B'
        cgbbDict[res].resName += 'B'
        if resDict[res].resName != 'GLY':
            cgscDict[res].resName += 'B'
    elif donor:
        blnResName += 'D'
        cgbbDict[res].resName += 'D'
        if resDict[res].resName != 'GLY':
            cgscDict[res].resName += 'D'
    elif acceptor:
        blnResName += 'A'
        cgbbDict[res].resName += 'A'
        if resDict[res].resName != 'GLY':
            cgscDict[res].resName += 'A'

    if debug:
        print 'blnResName = %s' % blnResName
    blnResList.append(blnResName)

    # finally, get the secondary structure information for this residue and change
    # the residue and atom names accordingly.
    myResType = ssList[resDict[res].resid-1]
    if myResType == 'helix':
        if debug:
            print "Residue type = helix"
        cgbbDict[res].resName = 'H%s' % cgbbDict[res].resName
        if resDict[res].resName != 'GLY':
            cgscDict[res].resName = 'H%s' % cgscDict[res].resName
    elif myResType == 'sheet':
        if debug:
            print "Residue type = sheet"
        cgbbDict[res].resName = 'S%s' % cgbbDict[res].resName
        if resDict[res].resName != 'GLY':
            cgscDict[res].resName = 'S%s' % cgscDict[res].resName
    else:
        if debug:
            print "Residue type = coil"
        cgbbDict[res].resName = 'C%s' % cgbbDict[res].resName
        if resDict[res].resName != 'GLY':
            cgscDict[res].resName = 'C%s' % cgscDict[res].resName

    cgList.append(cgbbDict[res])
    if resDict[res].resName != 'GLY':
        cgList.append(cgscDict[res])

# close charge stream file
csfp.write("return\n")
csfp.close()

if options.pdb_only:
    outputName = 'bln_' + pdbName + '.pdb'
    writeTo = open(outputPath + outputName,'w')
    for i,line in enumerate(cgList):
        line.atomNumber = i+1
        writeTo.write(line.Print('charmm'))
    writeTo.write('TER\n')
    writeTo.close()
    sys.exit()

# ( 3 ) Define more efficient Rij machinery
Rij = {}
def get_Rij(atom_i,atom_j):
    i = atom_i.atomNumber
    j = atom_j.atomNumber
    try:
        return Rij[(i,j)]
    except KeyError:
        Rij[(i,j)] = atom_i.bond_length(atom_j)
        return Rij[(i,j)]

# ( 4 ) Reindex CG.atomNumber
for i,cg in enumerate(cgList):
    cg.atomNumber = i+1

# ( 5 ) Write .crd file
outputName = 'bln_' + pdbName + '.crd'
writeTo = open(outputPath + outputName,'w')
writeTo.write('*\n')
writeTo.write('%5d\n' % len(cgList) )
for line in cgList:
    writeTo.write(line.Print('cgcrd'))
writeTo.close()

# ( 6 ) Write a .seq file
if debug:
    print "blnReslist = %s" % blnResList
outputName = 'bln_' + pdbName + '.seq'
writeTo = open(outputPath + outputName,'w')
writeTo.write(' '.join(blnResList))
writeTo.close()

# ( 7 ) Write .pdb file
outputName = 'bln_' + pdbName + '.pdb'
writeTo = open(outputPath + outputName,'w')
for line in cgList:
    writeTo.write(line.Print('charmm'))
writeTo.write('TER\n')
writeTo.close()

# ( 8 ) Write a .prm file
#   (0) CG Parameter Values
#   (a) Bond
#   (b) Angle
#   (d) Non-Bonded
#   (e) Backbone hydrogen bonding
#   (f) Native sidechain interaction
#   (g) Backbone sidechain interaction
#   (h) END of file information & Checksum

mainAtomTypes = ['P', 'N0', 'ND', 'NA', 'NDA', 'C', 'Q0', 'QD', 'QA', 'QDA']
compAtomTypes = ['CND', 'CBA', 'CQA']
atomSuffixes = ['','H','S']

outputName = 'bln_' + pdbName + '.prm'
writeTo = open(outputPath + outputName,'w')
writeTo.write('* This CHARMM .prm file describes a Go model of %s\n' % pdbName)
writeTo.write('* nScale             = %6.3f\n' % nScale)
writeTo.write('* domainScale        = %6.3f\n' % domainScale)
writeTo.write('* kBondHelix         = %6.3f\n' % kBondHelix)
writeTo.write('* kBondSheet         = %6.3f\n' % kBondSheet)
writeTo.write('* kBondCoil          = %6.3f\n' % kBondCoil)
writeTo.write('* kAngleHelix        = %6.3f\n' % kAngleHelix)
writeTo.write('* kAngleSheet        = %6.3f\n' % kAngleSheet)
writeTo.write('* kAngleCoil         = %6.3f\n' % kAngleCoil)
writeTo.write('*\n')
writeTo.write('\n')

writeTo.write('\n')                                    # (a)
writeTo.write('BOND\n')

#   all protein bonds have lengths of 3.8 angstrom
#   Lipid bonds (not supported yet) have 4.7 angstroms
finalTypeList = []
for atomType in mainAtomTypes + compAtomTypes:
    for suffix in atomSuffixes:
        finalTypeList.append(atomType + suffix)

if debug:
    print "Final list of atom types is:\n%s\n" % (finalTypeList)

for i in range(len(finalTypeList)):
    # bond parameters are divided by two b/c the Gromacs bond potential is
    # 1/2*Kb*(r-r0)**2 whereas CHARMM's in Kb*(r-r0)**2
    for j in range(i+1):
        if finalTypeList[j].strip().endswith('H'):
            stringBuffer = '%-8s%-8s  %12.6f%12.6f\n' % (finalTypeList[i], \
              finalTypeList[j],kBondHelix/2.0,3.8)
        elif finalTypeList[j].strip().endswith('S'):
            stringBuffer = '%-8s%-8s  %12.6f%12.6f\n' % (finalTypeList[i], \
              finalTypeList[j],kBondSheet/2.0,3.8)
        else:
            stringBuffer = '%-8s%-8s  %12.6f%12.6f\n' % (finalTypeList[i], \
              finalTypeList[j],kBondCoil/2.0,3.8)
        writeTo.write(stringBuffer)

writeTo.write('\n')                                    # (b)
writeTo.write('ANGLE\n')
for k in range(len(finalTypeList)):
    middleAtom = finalTypeList[k]
    if debug:
        print "middleAtom = %s" % middleAtom
    # angle parameters are NOT divided by 2 b/c I put this into the Marrink style angles
    for i in range(len(finalTypeList)):
        for j in range(i+1):
             if middleAtom.strip().endswith('H'):
                 stringBuffer = '%-8s%-8s%-8s  %12.6f%12.6f\n' % (\
                   finalTypeList[i],middleAtom.strip(),finalTypeList[j],kAngleHelix,-270)
             elif middleAtom.strip().endswith('S'):
                 stringBuffer = '%-8s%-8s%-8s  %12.6f%12.6f\n' % (\
                   finalTypeList[i],middleAtom.strip(),finalTypeList[j],kAngleSheet,-230)
             else:
                 stringBuffer = '%-8s%-8s%-8s  %12.6f%12.6f\n' % (\
                   finalTypeList[i],middleAtom.strip(),finalTypeList[j],kAngleCoil,-240)
             writeTo.write(stringBuffer)

writeTo.write('\n')                                    # (c)
writeTo.write('DIHEDRAL\n')
writeTo.write('! The Marrinck BLN model does not use dihedrals\n')
writeTo.write('\n')                                    # (d)
writeTo.write('NONBONDED NBXMOD 4 ATOM CDIEL SHIFT VATOM VDISTANCE VSHIFT -\n')
writeTo.write('CUTNB 16 CTOFNB 12 EPS 1.0 WMIN 1.5 E14FAC 0.7\n')
writeTo.write('!atom\n')
writeTo.write('!these should NEVER be used -- see NBFIX for REAL parameters\n')
writeTo.write('!we still need them though, otherwise NBFIX will not get used\n')

for atomType in finalTypeList:
    atomClass = (atomType.strip())[:-1]
    rMinDiv2 = 4.7 # This is what Marrink et al. set this to in their paper
    if atomClass in mainAtomTypes:
        epsilon  = -1.2 # FIXME: WRONG, NEEDS TO BE SET TO PREOPER EPS.
    else:
        epsilon  = -0.8 # FIXME: WRONG, NEEDS TO BE AVERAGED
    stringBuffer = '%-8s  %3.1f%10s%12.6f ! this parameter is garbage, only NBFIX params below should be used!\n' % (\
            atomType.strip(),0,epsilon,rMinDiv2)
    writeTo.write(stringBuffer)

writeTo.write('\n')                                    # (e)
writeTo.write('NBFIX\n')

for i in range(len(finalTypeList)):
    for j in range(i+1):
        # 4 cases must be dealt with:
        # (1) both atoms are main type
        # (2) atom 1 is main, 2 is composite
        # (3) reverse of (2)
        # (4) both atoms are composite
        atomClassI = finalTypeList[i].strip()
        atomClassJ = finalTypeList[j].strip()
        if atomClassI.endswith('H') or atomClassI.endswith('S'):
            atomClassI = atomClassI[:-1]
        if atomClassJ.endswith('H') or atomClassJ.endswith('S'):
            atomClassJ = atomClassJ[:-1]
        bondLength = 4.7 * 2**(1.0/6.0) # Sigma to Rmin conversion
        comment = "! pairwise contact between class %s and %s" % (atomClassI,atomClassJ)

        if debug:
            print "atomClassI = %s and atomClassJ = %s" % (atomClassI,atomClassJ)
        # Contact potentials are already in kcal/mol so they just need to be multiplied by -1
        ljDepth = get_contact_energy(atomClassI,atomClassJ) * -1.0

        stringBuffer = '%-8s%-8s  %12.6f%12.6f %s\n' % (\
            finalTypeList[i].strip(),finalTypeList[j].strip(),ljDepth,bondLength,comment)
        writeTo.write(stringBuffer)


writeTo.write('\n')                                    # (h)
writeTo.write('! Czech Sum Info:blahblah\n') # ToDo -- fill this in
writeTo.write('END\n')
writeTo.close()

