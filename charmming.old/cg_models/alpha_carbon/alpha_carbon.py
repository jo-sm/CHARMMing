#!/usr/bin/env python
# Filename CG_alpha_carbon.py
"""
# Outline
# ( 00) Option Parsing
# ( 01) Option Processing
# ( 1 ) Initialize All-Atom
# ( 2 ) Map All-Atom into CG
# ( 3 ) Assign each residue to a domain
# ( 4 ) Load contact potential
# ( 5 ) Define more efficient Rij machinery
# ( 6 ) Determine Native Contacts
# ( 7 ) STRIDE
# ( 8 ) Reindex CG.atomNumber
# ( 9 ) Write .crd file
# (10 ) Write a .seq file
# (11 ) Write .pdb file
# (12 ) Write .rtf file
# (13 ) Write a .prm file
"""

# The original implementation of this code was authored by:
# Edward P. O'Brien Jr., D. Thirumalai, B.R. Brooks groups, NIH/UMD
# 02/07/2007

# The CHARMMING implementation of this code was authored by:
# Frank C. Pickard IV, Benjamin T. Miller, Jing-Jun Sun, Henry F. Schaefer III,
# B.R. Brooks, H. Lee Woodcock NIH/UGA
# 08/18/2009

# Please cite the following article when referencing this model:
# Title: Effects of denaturants and osmolytes on proteins are accurately predicted
# by the molecular transfer model
# Authors: O'Brien, EP; Ziv, G; Haran, G; Brooks, BR; Thirumalai, D
# Source: PNAS | September 9, 2008 | vol. 105 | no. 36 | 13404-13408

import sys,commands,math,numpy,os,cPickle,string
import lib.Etc as Etc
from lib.Atom import Atom
from lib.Res import Pro
from lib.OptParse import OptionParser
from optparse import OptionGroup

# ( 00) Option Parsing
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
#contact set:['mj'|'kgs'|'bt']
parser.add_option('-C','--contacts',required=True,choices=['bt','kgs','mj'],metavar='CONTACT',
        help='A parameter set that defines all possible pairwise interactions between protein residues.\
                options are: (bt,kgs,mj)')
#input:formatting
parser.add_option('--formatting',default='auto',choices=['auto','charmm','web'],metavar='FORM',
        help='specify the .pdb formatting to expect: ([auto],charmm,web)')
#domain_string -- '23:54:107:203'
parser.add_option('--domain_string',default=None,metavar='STRING',
        help="A colon delimited STRING of integers which defines a mapping from a protein's residues to its domains.\
                ie. '23:54:107:203' -> res(1-23) in domain 1, res(24-54) in domain 2, etc. [Domain 1]")
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
################
###Parameters###
################
group = OptionGroup(parser,'Parameter Options',
        "All of the following options allow for (command line) fine tweaking of (almost) all of the parameters utilized by this script.")
   #nScale
group.add_option('--nScale',default=0.05,type='float',metavar='N',
        help='A user specified arbitrary energy scaling constant. [%default]')
   #domainScale
group.add_option('--domainScale',default=1.,type='float',metavar='N',
        help='A user specified arbitrary energy scaling factor for interactions between residues in different domains. [%default]')
   #fnn
group.add_option('--fnn',default=0.65,type='float',metavar='N',
        help='A user specified arbitrary vdW radii of side-chains. [%default]')
   #contactRad
group.add_option('--contactRad',default=4.5,type='float',metavar='R',
        help='Specify the sidechain contact radius in angstroms. [%default]')
   #kBond
group.add_option('--kBond',default=50.,type='float',metavar='K',
        help='Specify the spring constant for bonds. [%default]')
   #kAngle
group.add_option('--kAngle',default=30.,type='float',metavar='K',
        help='Specify the spring constant for angles. [%default]')
   #kDiheHelix_1
group.add_option('--kDiheHelix_1',default=0.3,type='float',metavar='K',
        help='Specify the spring constant of multiplicity 1 for dihedral angles in an alpha helix. [%default]')
   #kDiheHelix_3
group.add_option('--kDiheHelix_3',default=0.15,type='float',metavar='K',
        help='Specify the spring constant of multiplicity 3 for dihedral angles in an alpha helix. [%default]')
   #kDiheNoHelix_1
group.add_option('--kDiheNoHelix_1',default=0.55,type='float',metavar='K',
        help='Specify the spring constant of multiplicity 1 for dihedral angles in not in an alpha helix. [%default]')
   #kDiheNoHelix_3
group.add_option('--kDiheNoHelix_3',default=0.275,type='float',metavar='K',
        help='Specify the spring constant of multiplicity 3 for dihedral angles in not in an alpha helix. [%default]')
   #hBondEnergyHelix
group.add_option('--hBondEnergyHelix',default=-0.25,type='float',metavar='E',
        help='Specify the depth of an LJ potential between two helical CA. [%default]')
   #hBondEnergyNoHelix
group.add_option('--hBondEnergyNoHelix',default=-0.5,type='float',metavar='E',
        help='Specify the depth of an LJ potential between CA where one or both are not helical. [%default]')
   #epsilonNN
group.add_option('--epsilonNN',default=1e-12,type='float',metavar='E',
        help='Specify the depth of an LJ potential between two non-native contacts. [%default]')
   #bbscInteraction
group.add_option('--bbscInteraction',default=-0.37,type='float',metavar='E',
        help='Specify the depth of an LJ potential between backbone and side-chain beads. [%default]')
# Load Options
parser.add_option_group(group)
(options,args) = parser.parse_args(sys.argv)

# ( 01) Option Processing

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
# contact set
contactParmSet = options.contacts
# directory with contact set data
contactDir = Etc.expandPath(options.data)
# formatting
if options.formatting == 'auto':
    formatting = Etc.get_ver(inputPath + inputFile)
else:
    formatting = options.formatting
# domains
if options.domain_string:
    domainString = options.domain_string
else:
    domainString = False
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
fnn                 = options.fnn               #  0.65     # A scaling constant for the vdw radius of side-chains
contactRad          = options.contactRad        #  4.5      # Sidechain contact radius (angstrom)
kBond               = options.kBond             # 50.       # Spring constant for bonds
kAngle              = options.kAngle            # 30.       # Spring constant for angle bends
kDiheHelix_1        = options.kDiheHelix_1      #  0.300    # One of the spring constants for helical dihedrals
kDiheHelix_3        = options.kDiheHelix_3      #  0.150    # The other spring constant for helical dihedrals
kDiheNoHelix_1      = options.kDiheNoHelix_1    #  0.550    # One of the spring constants for non-helical dihedrals
kDiheNoHelix_3      = options.kDiheNoHelix_3    #  0.275    # The other spring constant for non-helical dihedrals
hBondEnergyHelix    = options.hBondEnergyHelix  # -0.25     # The depth of a LJ hbond potential between two helical centers 
hBondEnergyNoHelix  = options.hBondEnergyNoHelix# -0.50     # The depth of a LJ hbond potential between non-helical centers (not-helix+not-helix & not-helix+helix)
epsilonNN           = options.epsilonNN         #  1e-12    # The well depth for the LJ potential of non-native contacts (basically zero)
bbscInteraction     = options.bbscInteraction   # -0.37     # LJ well depth between backbone and side-chain
#######################
#Contact Parameter Set#
#######################
if contactParmSet == 'kgs':
    contactFileName  = '%s/kgs_contact_potential.dat' % contactDir
    contactOffset     = 1.8
elif contactParmSet == 'bt':
    contactFileName  = '%s/bt_contact_potential.dat' % contactDir
    contactOffset     = 0.6
elif contactParmSet == 'mj':
    contactFileName  = '%s/mj_contact_potential.dat' % contactDir
    contactOffset     = 1.2
else:
    raise DeadlyErrorOfDeath

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

# ( 2 ) Map All-Atom into CG
cgbbDict = {}
cgscDict = {}
cgList = []
for res in resList:
    cgbbDict[res] = resDict[res].get_alpha_carbon()
    cgbbDict[res].resCode = res
    cgList.append(cgbbDict[res])
    if resDict[res].resName != 'GLY':
        cgscDict[res] = resDict[res].get_sc()
        cgscDict[res].resCode = res
        cgList.append(cgscDict[res])

if options.pdb_only:
    outputName = 'cg_' + pdbName + '.pdb'
    writeTo = open(outputPath + outputName,'w')
    for i,line in enumerate(cgList):
        line.atomNumber = i+1
        writeTo.write(line.Print('cgcharmm'))
    writeTo.write('TER\n')
    writeTo.close()
    sys.exit()

# ( 3 ) Assign each residue to a domain
#    (a) Assign
if domainString:    # When explicitly specified
    domains = domainString.split(':')
    domains = map(int,domains)
    assert len(resDict) == domains[-1]
    thisDomain = 0
    for i,res in enumerate(resList):
        if i+1 <= domains[thisDomain]:
            pass
        else:
            thisDomain += 1
        resDict[res].domainid       = thisDomain
        cgbbDict[res].domainid      = thisDomain
        try:
            cgscDict[res].domainid  = thisDomain
        except KeyError: pass
else:               # Default: each segment gets it's own domain
    alpha2int   = dict( [ (char,i+1) for i,char in enumerate(string.uppercase) ] )
    for res in resList:
        resDict[res].domainid       = alpha2int[res[:1]]
        cgbbDict[res].domainid        = alpha2int[res[:1]]
        try:
            cgscDict[res].domainid  = alpha2int[res[:1]]
        except KeyError: pass
#   (b) Czech that all residues have a domain assignment
for res in resList:
    try:
        dummy = resDict[res].domainid
        dummy = cgbbDict[res].domainid
        try:
            dummy = cgscDict[res].domainid
        except KeyError: pass
    except AttributeError:
        print 'domains: residue %s has no domain assignment' % res
        raise GhastlyDeath

# ( 4 ) Load contact potential
#   (a) Parse Specified Contact File
resNameMap      = {}
contactMatrix   = []
for line in Etc.stripBare(open(contactFileName)):
    if line.startswith('AA'):
        buffer = line.split(':')[1].split()
        resNameMap = dict( zip(buffer,range(len(buffer))) )
    else:
        scaledContacts = [ nScale * abs(float(contact) - contactOffset) for contact in line.split() ]
        contactMatrix.append(scaledContacts)
#   (b) Calculate average interaction
avgEpsilon = sum(Etc.flatten(contactMatrix))/len(Etc.flatten(contactMatrix))
#   (c) Create ResName -> ContactEnergy mapper
def get_contact_energy(resName1,resName2):
    res1 = resNameMap[resName1]
    res2 = resNameMap[resName2]
    try:
        return float(contactMatrix[res1][res2])
    except IndexError:
        return float(contactMatrix[res2][res1])

# ( 5 ) Define more efficient Rij machinery
Rij = {}
def get_Rij(atom_i,atom_j):
    i = atom_i.segid + str(atom_i.atomNumber)
    j = atom_j.segid + str(atom_j.atomNumber)
    try:
        return Rij[(i,j)]
    except KeyError:
        Rij[(i,j)] = atom_i.bond_length(atom_j)
        return Rij[(i,j)]

# ( 6 ) Determine Native Contacts
#   (a) SC/SC Contacts
lowerOffAlmostDiagonal = ( (resCode_i,resCode_j) for i,resCode_i in enumerate(resList) for j,resCode_j in enumerate(resList) if j - i > 2 and resDict[resCode_i].resName != 'GLY' and resDict[resCode_j].resName != 'GLY' )
scContact = []
for resCode_i,resCode_j in lowerOffAlmostDiagonal:
    contact = False
    try:
        for atom_i in scDict[resCode_i]:
            for atom_j in scDict[resCode_j]:
                if get_Rij(atom_i,atom_j) < contactRad:
                    scContact.append( (cgscDict[resCode_i],cgscDict[resCode_j]) )
                    raise AssertionError
    except AssertionError: pass

#   (b) BB/SC Contacts
lowerOffAlmostDiagonal = ( (resCode_i,resCode_j) for i,resCode_i in enumerate(resList) for j,resCode_j in enumerate(resList) if abs(j - i) > 2 and resDict[resCode_j].resName != 'GLY' )
bbscContact = []
for resCode_i,resCode_j in lowerOffAlmostDiagonal:
    contact = False
    try:
        for atom_i in bbDict[resCode_i]:
            for atom_j in scDict[resCode_j]:
                if get_Rij(atom_i,atom_j) < contactRad:
                    bbscContact.append( (cgbbDict[resCode_i],cgscDict[resCode_j]) )
                    raise AssertionError
    except AssertionError:
        pass

# ( 7 ) STRIDE
strideOutput = commands.getstatusoutput('%s -h %s' % (stridePath,inputPath+inputFile) )[1].split('\n')
if strideOutput[0].startswith('Error'): raise IOError('stride: %s' % strideOutput[0])
#   (a) Get helical secondary structure
strideHelixOutput = [ line for line in strideOutput if line.startswith('ASG') ]
for i,line in enumerate(strideHelixOutput):
    if line[34:39] == 'Helix':
        cgbbDict[resList[i]].helix = True
    else:
        cgbbDict[resList[i]].helix = False

#   (b) Get hydrogen bonds
strideHBondOutput = [ line for line in strideOutput if line[:3] in ['ACC','DNR'] and int(line[15:20]) < int(line[35:40]) ]
#       (i) Multiplicity, ie anti-parallel beta sheets will have 'double' hbonds
hBondMult = {}
for line in strideHBondOutput:
    res_i = int(line[15:20])
    res_j = int(line[35:40])
    if (res_i,res_j) in hBondMult:
        hBondMult[(res_i,res_j)] += 1
    else:
        hBondMult[(res_i,res_j)] = 1
#       (ii) Build hbond list, 'double' bonds, get double energy
orderDict = ( (i,j) for i in range(len(resList)) for j in range(len(resList)) if i < j )
bbHBond = []
for (i,j) in orderDict:
    try:
        bb_i = cgbbDict[resList[i]]
        bb_j = cgbbDict[resList[j]]
        if bb_i.helix and bb_j.helix:
            hBondEnergy = hBondEnergyHelix
        else:
            hBondEnergy = hBondEnergyNoHelix
        bbHBond.append( (bb_i,bb_j,hBondEnergy,hBondMult[(i,j)]) )
    except KeyError:
        pass

# ( 8 ) Reindex CG.atomNumber
for i,cg in enumerate(cgList):
    cg.atomNumber = i+1

# ( 9 ) Write .crd file
outputName = 'cg_' + pdbName + '.crd'
writeTo = open(outputPath + outputName,'w')
writeTo.write('*\n')
writeTo.write('%5d\n' % len(cgList) )
for line in cgList:
    writeTo.write(line.Print('cgcrd'))
writeTo.close()

# (10 ) Write a .seq file
outputName = 'cg_' + pdbName + '.seq'
writeTo = open(outputPath + outputName,'w')
writeTo.write(' '.join(resList))
writeTo.close()

# (11 ) Write .pdb file
outputName = 'cg_' + pdbName + '.pdb'
writeTo = open(outputPath + outputName,'w')
for line in cgList:
    writeTo.write(line.Print('cgcharmm'))
writeTo.write('TER\n')
writeTo.close()

# (12 ) Write .rtf file
outputName = 'cg_' + pdbName + '.rtf'
writeTo = open(outputPath + outputName,'w')
writeTo.write('* This CHARMM .rtf file describes a Go model of %s\n' % pdbName)
writeTo.write('*\n')
writeTo.write('%5d%5d\n' % (20,1) )
writeTo.write('\n')

# Mass Section
for i,cg in enumerate(cgList):
    stringBuffer = 'MASS %-5d%-8s%10.6f\n' % (i+1,cg.resCode+cg.atomType.strip(),cg.mass)
    writeTo.write(stringBuffer)
writeTo.write('\n')

# Declare statements & Defaults
writeTo.write('DECL +B\n')
writeTo.write('DECL -B\n')
writeTo.write('DECL #B\n')
writeTo.write('DEFAULT FIRST NONE LAST NONE\n')
writeTo.write('\n')

# Residue Topology Section
for res in resList:
    if resDict[res].resName == 'GLY':
        writeTo.write('RESI %-5s       0.0\n' % res )
        writeTo.write('GROUP\n')
        writeTo.write('ATOM   B  %-5s  0.0\n' % (res+'B') )
        writeTo.write('BOND   B +B\n')
        writeTo.write('ANGLE -B  B +B\n')
        writeTo.write('DIHE  -B  B +B #B\n')
        writeTo.write('\n')
    else:
        writeTo.write('RESI %-5s       0.0\n' % res )
        writeTo.write('GROUP\n')
        writeTo.write('ATOM   B  %-5s  0.0\n' % (res+'B') )
        writeTo.write('ATOM   S  %-5s  0.0\n' % (res+'S') )
        writeTo.write('BOND   B  S  B +B\n')
        writeTo.write('ANGLE -B  B  S  S  B +B -B  B +B\n')
        writeTo.write('DIHE  -B  B +B #B\n')
        writeTo.write('IMPH   B -B +B  S\n')
        writeTo.write('\n')
writeTo.write('\n')
writeTo.write('END\n')
writeTo.close()

# (13 ) Write a .prm file
#   (0) CG Parameter Values
#   (a) Bond
#   (b) Angle
#   (c) Dihedral
#   (d) Improper Dihedral
#   (e) Non-Bonded
#   (f) Backbone hydrogen bonding
#   (g) Native sidechain interaction
#   (h) Backbone sidechain interaction
#   (i) END of file information & Checksum

outputName = 'cg_' + pdbName + '.prm'
writeTo = open(outputPath + outputName,'w')
writeTo.write('* This CHARMM .param file describes a Go model of %s\n' % pdbName)
writeTo.write('* contactParmSet     = %s   \n' % contactParmSet)
writeTo.write('* nScale             = %6.3f\n' % nScale)
writeTo.write('* domainScale        = %6.3f\n' % domainScale)
writeTo.write('* fnn                = %6.3f\n' % fnn)
writeTo.write('* contactRad         = %6.3f\n' % contactRad)
writeTo.write('* kBond              = %6.3f\n' % kBond)
writeTo.write('* kAngle             = %6.3f\n' % kAngle)
writeTo.write('* kDiheHelix_1       = %6.3f\n' % kDiheHelix_1)
writeTo.write('* kDiheHelix_3       = %6.3f\n' % kDiheHelix_3)
writeTo.write('* kDiheNoHelix_1     = %6.3f\n' % kDiheNoHelix_1)
writeTo.write('* kDiheNoHelix_3     = %6.3f\n' % kDiheNoHelix_3)
writeTo.write('* hBondEnergyHelix   = %6.3f\n' % hBondEnergyHelix)
writeTo.write('* hBondEnergyNoHelix = %6.3f\n' % hBondEnergyNoHelix)
writeTo.write('* epsilonNN          = %6.3e\n' % epsilonNN)
writeTo.write('* bbscInteraction    = %6.3f\n' % bbscInteraction)
writeTo.write('*\n')
writeTo.write('\n')

writeTo.write('\n')                                    # (a)
writeTo.write('BOND\n')

#   BB(i)/BB(i+1) length
for i,res in enumerate(resList):
    try:
        bondLength = cgbbDict[res].bond_length(cgbbDict[resList[i+1]])
        stringBuffer = '%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'B',resList[i+1]+'B',kBond,bondLength)
        writeTo.write(stringBuffer)
    except IndexError:
        pass

#   BB(i)/SC(i) length
for res in resList:
    try:
        bondLength = cgbbDict[res].bond_length(cgscDict[res])
        stringBuffer = '%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'B',res+'S',kBond,bondLength)
        writeTo.write(stringBuffer)
    except KeyError:
        pass

writeTo.write('\n')                                    # (b)
writeTo.write('ANGLE\n')

#   BB(i)/BB(i+1)/BB(i+2) angle
for i,res in enumerate(resList):
    try:
        bondAngle = cgbbDict[res].bond_angle(cgbbDict[resList[i+1]],cgbbDict[resList[i+2]],'deg')
        stringBuffer = '%-8s%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'B',resList[i+1]+'B',resList[i+2]+'B',kAngle,bondAngle)
        writeTo.write(stringBuffer)
    except IndexError:
        pass

#   SC(i)/BB(i)/BB(i+1) angle
for i,res in enumerate(resList):
    try:
        bondAngle = cgscDict[res].bond_angle(cgbbDict[res],cgbbDict[resList[i+1]],'deg')
        stringBuffer = '%-8s%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'S',res+'B',resList[i+1]+'B',kAngle,bondAngle)
        writeTo.write(stringBuffer)
    except IndexError:
        pass
    except KeyError:
        pass

#   BB(i)/BB(i+1)/SC(i+1) angle
for i,res in enumerate(resList):
    try:
        bondAngle = cgbbDict[res].bond_angle(cgbbDict[resList[i+1]],cgscDict[resList[i+1]],'deg')
        stringBuffer = '%-8s%-8s%-8s  %12.6f%12.6f\n' % (\
                res+'B',resList[i+1]+'B',resList[i+1]+'S',kAngle,bondAngle)
        writeTo.write(stringBuffer)
    except IndexError:
        pass
    except KeyError:
        pass

writeTo.write('\n')                                    # (c)
writeTo.write('DIHEDRAL\n')
writeTo.write('! Backbone\n')

#   BB(i)/BB(i+1)/BB(i+2)/BB(i+3) dihedral
for i,res in enumerate(resList):
    try:
        dihedral = cgbbDict[res].bond_signed_dihedral(cgbbDict[resList[i+1]],cgbbDict[resList[i+2]],cgbbDict[resList[i+3]],'deg')
        if cgbbDict[resList[i+1]].helix and cgbbDict[resList[i+2]].helix:
            kDihe   = ['',kDiheHelix_1,'',kDiheHelix_3]
        else:
            kDihe   = ['',kDiheNoHelix_1,'',kDiheNoHelix_3]

        for multiplicity in [1,3]:
            delta       = Etc.dihedral_mod( multiplicity * dihedral - 180. )
            kDihedral   = kDihe[multiplicity]

            stringBuffer = '%-8s%-8s%-8s%-8s  %12.6f%3d%12.6f\n' % (\
                    res+'B',resList[i+1]+'B',resList[i+2]+'B',resList[i+3]+'B',kDihedral,multiplicity,delta)
            writeTo.write(stringBuffer)
    except IndexError:
        pass

writeTo.write('\n')                                    # (d)
writeTo.write('IMPROPERS\n')
writeTo.write('! Sidechain\n')

#   BB(i+1)/BB(i)/BB(i+2)/SC(i+1) dihedral
for i,res in enumerate(resList):
    try:
        dihedral = cgbbDict[resList[i+1]].bond_signed_dihedral(cgbbDict[res],cgbbDict[resList[i+2]],cgscDict[resList[i+1]],'deg')
        multiplicity = 1
        kDihedral    = 20. * abs(avgEpsilon)
        delta        = dihedral + 180.
        stringBuffer = '%-8s%-8s%-8s%-8s  %12.6f%3d%12.6f\n' % (\
                resList[i+1]+'B',res+'B',resList[i+2]+'B',resList[i+1]+'S',kDihedral,multiplicity,delta)
    except IndexError:
        pass
    except KeyError:
        pass

writeTo.write('\n')                                    # (e)
writeTo.write('NONBONDED NBXMOD 4 ATOM CDIEL SHIFT VATOM VDISTANCE VSWITCH -\n')
writeTo.write('CUTNB 23 CTOFNB 20 CTONNB 18 EPS 1.0 WMIN 1.5 E14FAC 0.7\n')
writeTo.write('!atom\n')

for res in cgList:
    if res.atomType == '   B':
        rMinDiv2 = 20.
    elif res.atomType == '   S':
        rMinDiv2 = 10. * Etc.get_aa_vdw( res.resName ) * 2.**(1./6.) * fnn
    else:
        raise DeadlyErrorOfDeath
    stringBuffer = '%-8s  %3.1f%10s%12.6f\n' % (\
            res.resCode+res.atomType.strip(),0,str(-1*epsilonNN),rMinDiv2)
    writeTo.write(stringBuffer)

writeTo.write('\n')                                    # (f)
writeTo.write('NBFIX\n')
writeTo.write('! backbone hydrogen bonding\n')

bbEnergySum = 0
for hBond in bbHBond:
    bondLength  = hBond[0].bond_length(hBond[1])
    comment     = ''
    hBondEnergy = hBond[2]
    if hBond[0].domainid != hBond[1].domainid:
        comment     += '! Interface between domains %d, %d ' % (hBond[0].domainid,hBond[1].domainid)
        hBondEnergy *= domainScale
    if hBond[3] != 1:
        comment     += '! hBond multiplcity is %d ' % hBond[3]
        hBondEnergy *= hBond[3]
    bbEnergySum += hBondEnergy
    stringBuffer = '%-8s%-8s  %12.6f%12.6f %s\n' % (\
            hBond[0].resCode+'B',hBond[1].resCode+'B',hBondEnergy,bondLength,comment)
    writeTo.write(stringBuffer)

writeTo.write('! native side-chain interactions\n')    # (g)
scEnergySum = 0
for contact in scContact:
    bondLength  = contact[0].bond_length(contact[1])
    comment     = ''
    ljDepth     = -1. * get_contact_energy( contact[0].resName,contact[1].resName )
    if contact[0].domainid != contact[1].domainid:
        comment     += '! Interface between domains %d, %d ' % (contact[0].domainid,contact[1].domainid)
        ljDepth     *= domainScale/nScale
    scEnergySum += ljDepth
    stringBuffer = '%-8s%-8s  %12.6f%12.6f %s\n' % (\
            contact[0].resCode+'S',contact[1].resCode+'S',ljDepth,bondLength,comment)
    writeTo.write(stringBuffer)

writeTo.write('! backbone side-chain interactions\n')  # (h)
bbscEnergySum = 0
for contact in bbscContact:
    bondLength  = contact[0].bond_length(contact[1])
    comment     = ''
    ljDepth     = bbscInteraction
    if contact[0].domainid != contact[1].domainid:
        comment     += '! Interface between domains %d, %d ' % (contact[0].domainid,contact[1].domainid)
        ljDepth     *= domainScale
    bbscEnergySum += ljDepth
    stringBuffer = '%-8s%-8s  %12.6f%12.6f %s\n' % (\
            contact[0].resCode+'B',contact[1].resCode+'S',ljDepth,bondLength,comment)
    writeTo.write(stringBuffer)

writeTo.write('\n')                                    # (i)
writeTo.write('! Czech Sum Info:%5d,%8.2f,%8.2f\n' % (bbEnergySum,scEnergySum,bbscEnergySum) )
writeTo.write('END\n')
writeTo.close()

