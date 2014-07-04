# charmming v1 rough cut
# fcp
# 06/20/2010

import math,numpy,string
from copy import copy
import lib.Etc as Etc

RAD2DEG = 180/math.pi
DEG2RAD = math.pi/180

class AtomError(Exception):
    """
    Exception to raise when errors occur involving the Atom class.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class BaseAtom(object):
    """
    This class defines data, meta data and methods common to all atom like classes.
    # Metadata
        text        = string
        _strict     = Bool (True)
        _debug      = Bool (False)
    # Methods
        get_addr
        get_hash
        czech_all
        set_tag         (6*char)
        set_atomNum     (0 <= int < 100000)
        set_atomType    (4*char)
        set_resName     (4*char)
        set_segid       (char)
        set_resid       (0 <= int < 10000)
        set_xCart       (float)
        set_yCart       (float)
        set_xCart       (float)
        set_cart        (float,float,float)
        Print
        calc_length     (other)
        calc_angle      (other1,other2)
        calc_dihedral   (other1,other2,other3)
        calc_signedDihedral (other1,other2,other3)
    """
    def __init__(self,text='',index=0,strict=True,autoFix=True,debug=False,CC='#'):
        # Metadata
        self.text = text.lower().split(CC)[0]
        self._strict = strict
        self._autoFix = True
        self._debug = debug
        self._index = index
        # Data Defaults
        self._tag       = 'TAG'
        self._atomNum   = -1
        self._atomType  = 'aTyp'
        self._resName   = 'RES'
        self._segid     = '!'
        self._resid     = -1
        self._xCart     = 9999.999
        self._yCart     = 9999.999
        self._zCart     = 9999.999
        self._cart      = numpy.array([self._xCart,self._yCart,self._zCart])
        self._mass      = 0.

    def get_addr(self):
        return '%s%d%d' % (self._segid,self._resid,self._atomNum)

    def get_hash(self):
        try:
            return '%s%d%d' % (self._segid0,self._resid0,self._atomNum0)
        except AttributeError:
            return 'index-%d' % self._index

    def czech_tag(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._tag
        if not len(self._tag) <= 6:
            msg = 'czech_tag %s: tag %s is longer than 6 characters' % (self.get_hash(),testMe)
            if self._strict:
                raise AtomError(msg)
            else:
                print msg

    def czech_atomNum(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._atomNum
        if not 0 <= testMe < 100000:
            msg = 'czech_atomNum %s: atomNum %d is out of range (0,100000)' % (self.get_hash(),testMe)
            if self._strict:
                raise AtomError(msg)
            else:
                print msg

    def czech_atomType(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._atomType
        if len(testMe) > 4:
            msg = 'czech_atomType %s: atomType %s is longer than 4 characters' % (self.get_hash(),testMe)
            if self._strict:
                raise AtomError(msg)
            else:
                print msg

    def czech_resName(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._resName
        if len(testMe) > 4:
            msg = 'czech_resName %s: resName %s is longer than 4 characters' % (self.get_hash(),testMe)
            if self._strict:
                raise AtomError(msg)
            else:
                print msg

    def czech_segid(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._segid
        if len(testMe) > 1:
            msg = 'czech_segid %s: segid %s is longer than 1 character' % (self.get_hash(),testMe)
            if self._strict:
                raise AtomError(msg)
            else:
                print msg

    def czech_resid(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._resid
        if not -1000 <= testMe < 10000:
            msg = 'czech_resid %s: resid %d is out of range [-1000,10000)' % (self.get_hash(),testMe)
            if self._strict:
                raise AtomError(msg)
            else:
                print msg

    def czech_weight(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._weight
        if not 0. <= testMe <= 1.:
            msg = 'czech_weight %s: weight %4.2f is out of range [0,1]' % (self.get_hash(),testMe)

    def czech_bFactor(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._bFactor
        if not 0. <= testMe <= 100.:
            msg = 'czech_weight %s: weight %5.2f is out of range [0,100]' % (self.get_hash(),testMe)

    def czech_mass(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._mass
        if not 0. < testMe < 1000.:
            msg = 'czech_mass %s: mass %6.2f is out of range [0,1000)' % (self.get_hash(),testMe)
            if self._strict:
                raise AtomError(msg)
            else:
                print msg

    def czech_element(self,arg=None):
        if arg:
            testMe = arg
        else:
            testMe = self._element
        if not testMe.upper() in Etc.atomMass.keys():
            msg = 'czech_element %s: element name %s is not valid' % (self.get_hash(),testMe)
            if self._strict:
                raise AtomError(msg)
            else:
                print msg

    def set_tag(self,arg):
        arg = str(arg).strip()
        self.czech_tag(arg)
        self._tag = arg

    def set_atomNum(self,arg):
        arg = int(arg)
        self.czech_atomNum(arg)
        self._atomNum = arg

    def set_atomType(self,arg):
        arg = str(arg)
        self.czech_atomType(arg)
        self._atomType = arg

    def set_resName(self,arg):
        arg = str(arg).strip()
        self.czech_resName(arg)
        self._resName = arg

    def set_segid(self,arg):
        arg = str(arg)
        try:
            self.czech_segid(arg)
        except AtomError:
            arg = arg[0]
        self._segid = arg

    def set_resid(self,arg):
        arg = int(arg)
        self.czech_resid(arg)
        self._resid = arg

    def set_mass(self,arg):
        arg = float(arg)
        self.czech_mass(arg)
        self._mass = arg

    def set_xCart(self,arg):
        self._xCart = float(arg)
        self._cart  = numpy.array([self._xCart,self._yCart,self._zCart])

    def set_yCart(self,arg):
        self._yCart = float(arg)
        self._cart  = numpy.array([self._xCart,self._yCart,self._zCart])

    def set_zCart(self,arg):
        self._zCart = float(arg)
        self._cart  = numpy.array([self._xCart,self._yCart,self._zCart])

    def set_cart(self,arg1,arg2,arg3):
        self._xCart = float(arg1)
        self._yCart = float(arg2)
        self._zCart = float(arg3)
        self._cart  = numpy.array([self._xCart,self._yCart,self._zCart])

    def set_weight(self,arg):
        try:
            arg = float(arg)
        except ValueError:
            arg = 0.
        if self._autoFix:
            if arg >= 1.:
                arg = 1.
            elif arg < 0.:
                arg = 0.
        self.czech_weight(arg)
        self._weight = arg

    def set_bFactor(self,arg):
        try:
            arg = float(arg)
        except ValueError:
            arg = 99.9
        if self._autoFix:
            if arg >= 100.:
                arg = 99.9
            elif arg < 0.:
                arg = 0.
        self.czech_bFactor(arg)
        self._bFactor = arg

    def set_element(self,arg):
        arg = str(arg).strip()
        try:
            self.czech_element(arg)
        except AtomError:
            if self._autoFix:
                self._element = 'X'
            else:
                self.czech_element(arg)
        self._element = arg

    def Print(self,format='debug'):
        if format == 'debug':
            taco = 'debug: %d %4s%4s %1s%4d   %10.3f %10.3f %10.3f' % \
                (self._atomNum,self._atomType,self._resName,self._segid,self._resid,self._cart[0],self._cart[1],self._cart[2])
        else:
            raise TypeError('Print %s: unknown format' % self.get_hash())
        return taco

    def calc_length(self,other):
        """
        Returns the cartesian distance between two baseAtom objects.
        """
        if not issubclass(type(other),baseAtom):
            raise TypeError('calc_length: Second argument must be derived from baseAtom class.')
        temp = self._cart - other._cart
        return math.sqrt( numpy.dot(temp,temp) )

    def calc_angle(self,other1,other2,units='deg'):
        """
        Returns the bond calc_angle between three baseAtom objects.  You can optionally
        specify the units as the 3rd argument with 'rad', default is 'deg'.
        """
        if not issubclass(type(other1),baseAtom):
            raise TypeError('calc_angle: Second argument must be derived from baseAtom class.')
        if not issubclass(type(other2),baseAtom):
            raise TypeError('calc_angle: Third argument must be derived from baseAtom class.')
        # for 3-others (self,other1,other2) -> (i,j,k)
        Eji = (other1._cart -   self._cart) /   self.calc_length(other1)
        Ejk = (other1._cart - other2._cart) / other1.calc_length(other2)
        if   units == 'deg':
            return math.acos(numpy.dot(Eji,Ejk)) * RAD2DEG
        elif units in ['rad','au']:
            return math.acos(numpy.dot(Eji,Ejk))

    def calc_dihedral(self,other1,other2,other3,units='deg'):
        """
        Returns the torsion calc_angle between four baseAtom objects.  You can optionally
        specify the units as the 4th argument with 'rad', default is 'deg'.
        """
        if not issubclass(type(other1),baseAtom):
            raise TypeError('calc_dihedral: Second argument must be derived from baseAtom class.')
        if not issubclass(type(other2),baseAtom):
            raise TypeError('calc_dihedral: Third argument must be derived from baseAtom class.')
        if not issubclass(type(other3),baseAtom):
            raise TypeError('calc_dihedral: Forth argument must be derived from baseAtom class.')
        # for 4-others (self,other1,other2,other3) -> (i,j,k,l)
        # Unit vectors defining the 3 bonds
        Eij = (  self._cart - other1._cart) /   self.calc_length(other1)
        Ejk = (other1._cart - other2._cart) / other1.calc_length(other2)
        Ekl = (other2._cart - other3._cart) / other2.calc_length(other3)
        # Unit vectors normal to the two planes   
        Nijk = numpy.cross(Eij,Ejk) / math.sin(   self.calc_angle(other1,other2,'rad') )
        Njkl = numpy.cross(Ejk,Ekl) / math.sin( other1.calc_angle(other2,other3,'rad') )
        if   units == 'deg':
            return math.acos( numpy.dot(Nijk,Njkl) ) * RAD2DEG
        elif units in ['rad','au']:
            return math.acos( numpy.dot(Nijk,Njkl) )

    def calc_signedDihedral(self,other1,other2,other3,units='deg'):
        """
        Returns the signed torsion calc_angle between four baseAtom objects.  You can optionally
        specify the units as the 4th argument with 'rad', default is 'deg'.
        """
        if not issubclass(type(other1),baseAtom):
            raise TypeError('calc_signedDihedral: Second argument must be derived from baseAtom class.')
        if not issubclass(type(other2),baseAtom):
            raise TypeError('calc_signedDihedral: Third argument must be derived from baseAtom class.')
        if not issubclass(type(other3),baseAtom):
            raise TypeError('calc_signedDihedral: Forth argument must be derived from baseAtom class.')
        # for 4-others (self,other1,other2,other3) -> (i,j,k,l)
        # Unit vectors defining the 3 bonds
        Eij = (  self._cart - other1._cart) /   self.calc_length(other1)
        Ejk = (other1._cart - other2._cart) / other1.calc_length(other2)
        Ekl = (other2._cart - other3._cart) / other2.calc_length(other3)
        # Unit vectors normal to the two planes   
        Nijk = numpy.cross(Eij,Ejk) / math.sin(   self.calc_angle(other1,other2,'rad') )
        Njkl = numpy.cross(Ejk,Ekl) / math.sin( other1.calc_angle(other2,other3,'rad') )
        # Determine sign
        if numpy.dot( numpy.cross(Nijk,Njkl), Ejk) > 0:
            sign = -1
        else:
            sign =  1
        if   units == 'deg':
            return sign * math.acos( numpy.dot(Nijk,Njkl) ) * RAD2DEG
        elif units in ['rad','au']:
            return sign * math.acos( numpy.dot(Nijk,Njkl) )

###################
# Special Methods #
###################

    def __str__(self):
        return self.Print('debug')

    def __hash__(self):
        """
        _segid > _resid > _atomNum
        """
        segidScore = dict([(char,i+1) for i,char in enumerate(string.lowercase)])
        return segidScore[self._segid] * 1E9 + self._resid * 1E5 + self._atomNum

    def __eq__(self,other):
        return self.__hash__() == other.__hash__()

    def __ne__(self,other):
        return self.__hash__() != other.__hash__()

    def __lt__(self,other):
        return self.__hash__() < other.__hash__()

    def __le__(self,other):
        return self.__hash__() <= other.__hash__()

    def __gt__(self,other):
        return self.__hash__() > other.__hash__()

    def __ge__(self,other):
        return self.__hash__() >= other.__hash__()

###################
# Private Methods #
###################

    def _save_initialData(self):
        """
        Store initial values for crucial atom metadata.
        """
        self._segid0 = copy(self._segid)
        self._resid0 = copy(self._resid)
        self._atomNum0 = copy(self._atomNum)
