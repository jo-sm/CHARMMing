#!/usr/bin/env python
# Filename: Atom.py

import numpy,math,string
import Etc

class AtomError(Exception):
    """Exception to raise when errors occur involving the Atom class."""
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Atom(object):
    """An Atom object is parsed from a line of a pdb file begining with
    'ATOM' or 'HETATM'. The variables are parsed from their character offset
    relative to the begining of the line."""
    def __init__(self,line,ver='web',debug=False):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            pass
        else:
            raise AtomError('constructor: Atom line should start with "ATOM" or "HETATM"')
        # Raw Data
        self.line       = line.upper()
        self.ver        = ver
        self.debug      = debug
        # Parsed Data
        self.set_vars()
        # Initialize
        if not self.element: self.set_element()
        self.model = 0
        self.set_segid()
        self.set_mass()
        self.set_tag()
        self.cart       = numpy.array([self.xCart,self.yCart,self.zCart])
        self.resIndex   = self.resid
        del self.xCart
        del self.yCart
        del self.zCart

    def set_vars(self):
        """Parses the raw string input into instance variables."""
        if self.ver == 'web':
            self.tag        =       self.line[ 0: 6].strip()
            self.atomNumber = int(  self.line[ 6:11]    )
            self.atomType   =       self.line[12:16]
            self.resName    =       self.line[16:20].strip()
            self.segid      =       self.line[21:22]
            self.resid      = int(  self.line[22:26]    )
            self.xCart      = float(self.line[30:38]    )
            self.yCart      = float(self.line[38:46]    )
            self.zCart      = float(self.line[46:54]    )
            try:
                self.weight = float(self.line[55:60]    )
            except ValueError:
                self.weight = float(0)
            try:
                self.bFactor= float(self.line[61:66]    )
            except ValueError:
                self.bFactor= float(0)
            try:
                self.element=       self.line[66:  ].strip()
            except ValueError:
                pass

        elif self.ver == 'charmm':
            self.tag        =       self.line[ 0: 6].strip()
            self.atomNumber = int(  self.line[ 6:11]    )
            self.atomType   =       self.line[12:16]
            self.resName    =       self.line[17:21].strip()
            self.resid      = int(  self.line[22:26]    )
            self.xCart      = float(self.line[30:38]    )
            self.yCart      = float(self.line[38:46]    )
            self.zCart      = float(self.line[46:54]    )
            try:
                self.weight = float(self.line[55:60]    )
            except ValueError:
                self.weight = float(0)
            try:
                self.bFactor= float(self.line[61:66]    )
            except ValueError:
                self.bFactor= float(0)
            try:
                self.segid  =       self.line[72:76].strip()
            except ValueError:
                pass
            self.element    = ''
        else:
            raise AtomError('set_vars: self.ver == %s'%self.ver)
        return

    def get_resAddr(self):
        return '%d%s%s%d' % (self.model,self.type,self.segid,self.resid)

    def get_atomAddr(self):
        return '%d%s%s%d%d' % (self.model,self.type,self.segid,self.resid,self.atomNumber)

    def set_segid(self):
        """If an Atom's segid is blank, set it to 'A'."""
        if self.segid in [' ','']:
            self.segid = 'A'
        elif self.segid in string.digits:
            self.segid = string.uppercase[int(self.segid)-1]
        return

    def set_element(self):
        """Verify the Atom's element, or derive it from the atomType."""
        if Etc.isNuc( self.resName ) or Etc.isPro( self.resName ):  # If it is a pro or a nuc, grab the first letter from the type
            self.element = self.atomType.strip()[:1]
        else:
            self.element = self.atomType[:2].strip()                # Else grab the first two
        return

    def set_mass(self):
        """Take an Atom's element, and set it's mass if possible."""
        if self.element in Etc.atomMass:
            self.mass = Etc.atomMass[self.element]
        else:
            self.mass = None
            if self.debug:
                print 'set_mass: self.mass == None (segid = %s, resid = %5i, atomNumber = %6i)' % (self.segid,self.resid,self.atomNumber)
        return

    def make_compliant(self):
        """Change resName and atomType to be CHARMM compliant."""
        # FCP: This function handles charmm renaming for everything _except_
        # terminal residue renaming. This is because this class doesn't know if its
        # members are in the terminal residue. This is ugly, we may want to change it.
        if self.resName in ['HOH','TIP3']:                              # Water renaming
           self.resName = 'TIP3'
           self.atomType = 'OH2 '
        elif self.resName == 'HIS':                                     # HIS -> HSE
           self.resName = 'HSD'
        elif self.resName in ['A','DA']:
           self.resName = 'ADE'
        elif self.resName in ['T','DT']:
           self.resName = 'THY'
        elif self.resName in ['C','DC']:
           self.resName = 'CYT'
        elif self.resName in ['G','DG']:
           self.resName = 'GUA'
        elif self.resName in ['U','DU']:
           self.resName = 'URA'
        elif (self.resName == 'ILE' and self.atomType == ' CD1'):       # CD1 -> CD for res ILE
           self.atomType = ' CD '
        elif (self.atomType == 'ZN  ' and self.resName == 'ZN'):        # Ion renaming scheme
           self.resName = 'ZN2'
        elif (self.atomType == 'NA  ' and self.resName == 'NA'):
           self.resName = 'SOD'
           self.atomType = 'SOD '
        elif (self.atomType == 'CS  ' and self.resName == 'CS'):
           self.resName = 'CES'
           self.atomType = 'CES '
        elif (self.atomType == 'CL  ' and self.resName == 'CL'):
           self.resName = 'CLA'
           self.atomType = 'CLA '
        elif (self.atomType == 'CA  ' and self.resName == 'CA'):
           self.resName = 'CAL'
           self.atomType = 'CAL '
        elif (self.atomType ==  ' K  ' and self.resName ==  'K'):
           self.resName = 'POT'
           self.atomType = 'POT '
        return
    
    def set_tag(self):
        """Check each Atom's resName, and set tag to 'nuc','pro','good','bad'."""
        if Etc.isNuc(self.resName):
            self.type = 'nuc'
        elif Etc.isPro(self.resName):
            self.type = 'pro'
        elif Etc.isLip(self.resName):
            self.type = 'lip'
        elif Etc.isCrb(self.resName):
            self.type = 'crb'
        elif Etc.isGood(self.resName):
            self.type = 'good'
        else:
            self.type = 'bad'
        return

    def Print(self,format):
        if format == 'web':
            taco = '%-6s%5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n' % \
                (self.tag,self.atomNumber,self.atomType,self.resName,self.segid,self.resid,self.cart[0],self.cart[1],self.cart[2],self.weight,self.bFactor,self.element)
        elif format == 'charmm':
            taco = '%-6s%5i %4s %4s %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n' % \
                (self.tag,self.atomNumber,self.atomType,self.resName,self.resid,self.cart[0],self.cart[1],self.cart[2],self.weight,self.bFactor,self.segid)
        elif format == 'debug':
            taco = '%-6s%5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s' % \
                (self.type,self.atomNumber,self.atomType,self.resName,self.segid,self.resid,self.cart[0],self.cart[1],self.cart[2],self.weight,self.bFactor,self.element)
        elif format in ['crd','cor','card']:
            taco = '%5i%5i %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4i%10.5f\n' % \
                (self.atomNumber,self.resIndex,self.resName,self.atomType,self.cart[0],self.cart[1],self.cart[2],self.segid,self.resid,self.weight)
        elif format in ['xcrd','xcor','xcard']:
            taco = '%10i%10i  %-4s      %-4s    %20.10f%20.10f%20.10f  %-5s    %-4i    %20.10f\n' % \
                (self.atomNumber,self.resIndex,self.resName,self.atomType,self.cart[0],self.cart[1],self.cart[2],self.segid,self.resid,self.weight)
        elif format in ['cgcrd','cgcor','cgcard']:
            taco = '%5i%5i %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4i%10.5f\n' % \
                (self.atomNumber,self.resIndex,self.segid[0]+str(self.resid),self.atomType,self.cart[0],self.cart[1],self.cart[2],self.segid,self.resid,self.weight)
        elif format == 'cgcharmm':
            taco = '%-6s%5i %4s %4s %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n' % \
                (self.tag,self.atomNumber,self.atomType,self.segid[0]+str(self.resid),self.resid,self.cart[0],self.cart[1],self.cart[2],self.weight,self.bFactor,self.segid)
        else:
            raise AtomError('Print: unknown format == %s'%format)
        return taco

    def __repr__(self):
        return self.Print('debug')

    def __str__(self):
        return self.Print(self.ver)

    def bond_length(self,atom):
        """Returns the cartesian distance between two pdb_line objects."""
        temp = self.cart - atom.cart
        return math.sqrt( numpy.dot(temp,temp) )

    def bond_angle(self,atom1,atom2,units='deg'):
        """Returns the bond angle between three pdb_line objects.  You can optionally
        specify the units as the 3rd argument with 'rad', default is 'deg'."""
        # for 3-atoms (self,atom1,atom2) -> (i,j,k)
        Eji = (atom1.cart -  self.cart) /  self.bond_length(atom1)
        Ejk = (atom1.cart - atom2.cart) / atom1.bond_length(atom2)
        if   units == 'deg':
            return math.acos(numpy.dot(Eji,Ejk)) * Etc.RAD2DEG
        elif units in ['rad','au']:
            return math.acos(numpy.dot(Eji,Ejk))

    def bond_dihedral(self,atom1,atom2,atom3,units='deg'):
        """Returns the bond torsion between four pdb_line objects.  You can optionally
        specify the units as the 4th argument with 'rad', default is 'deg'.
        # for 4-atoms (self,atom1,atom2,atom3) -> (i,j,k,l)
        """
        # Unit vectors defining the 3 bonds
        Eij = (self.cart  - atom1.cart) /  self.bond_length(atom1)
        Ejk = (atom1.cart - atom2.cart) / atom1.bond_length(atom2)
        Ekl = (atom2.cart - atom3.cart) / atom2.bond_length(atom3)
        # Unit vectors normal to the two planes   
        Nijk = numpy.cross(Eij,Ejk) / math.sin(  self.bond_angle(atom1,atom2,'rad') )
        Njkl = numpy.cross(Ejk,Ekl) / math.sin( atom1.bond_angle(atom2,atom3,'rad') )
        if   units == 'deg':
            return math.acos( numpy.dot(Nijk,Njkl) ) * Etc.RAD2DEG
        elif units in ['rad','au']:
            return math.acos( numpy.dot(Nijk,Njkl) )

    def bond_signed_dihedral(self,atom1,atom2,atom3,units='deg'):
        """Returns the signed bond torsion between four pdb_line objects.  You can optionally
        specify the units as the 4th argument with 'rad', default is 'deg'.
        # for 4-atoms (self,atom1,atom2,atom3) -> (i,j,k,l)
        """
        # Unit vectors defining the 3 bonds
        Eij = (self.cart  - atom1.cart) /  self.bond_length(atom1)
        Ejk = (atom1.cart - atom2.cart) / atom1.bond_length(atom2)
        Ekl = (atom2.cart - atom3.cart) / atom2.bond_length(atom3)
        # Unit vectors normal to the two planes   
        Nijk = numpy.cross(Eij,Ejk) / math.sin(  self.bond_angle(atom1,atom2,'rad') )
        Njkl = numpy.cross(Ejk,Ekl) / math.sin( atom1.bond_angle(atom2,atom3,'rad') )
        # Determine sign
        if numpy.dot( numpy.cross(Nijk,Njkl), Ejk) > 0:
            sign = -1
        else:
            sign =  1
        if   units == 'deg':
            return sign * math.acos( numpy.dot(Nijk,Njkl) ) * Etc.RAD2DEG
        elif units in ['rad','au']:
            return sign * math.acos( numpy.dot(Nijk,Njkl) )

