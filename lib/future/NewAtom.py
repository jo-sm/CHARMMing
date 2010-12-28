#!/usr/bin/env python

from lib.future.BaseAtom import BaseAtom
from lib.future.BaseAtom import AtomError
import string,math
from copy import copy
import lib.Etc as Etc

class NewAtom(BaseAtom):
    """
    """
    def __init__(self,text='',format='web',index=0,modelNum=0,strict=True,CC='!'):
        # Setup
        BaseAtom.__init__(self,text,index,strict,CC)
        self._format    = format
        self.note      = ''
        # Defaults
        self._modelNum  = modelNum
        self._weight    = 0.
        self._bFactor   = 0.
        self._element   = '$'
        # Parse
        if text:
            self.parse()
            self.set_mass()
        self.set_segType()
        if not self._element: self._fix_element()
        self._resIndex = self._resid
        self._save_initialData()

##################
# Public Methods #
##################

    def parse(self):
        """
        Parses _text into instance variables.
        """
        if self._format == 'web':
            self.set_tag(self.text[0:6])
            self.set_atomNum(self.text[ 6:11])
            self.set_atomType(self.text[12:16])
            self.set_resName(self.text[16:20])
            self.set_segid(self.text[21])
            self.set_resid(self.text[22:26])
            self.set_cart(self.text[30:38],self.text[38:46],self.text[46:54])
            self.set_weight(self.text[55:60])
            self.set_bFactor(self.text[61:66])
            self.set_element(self.text[66:])
        elif self._format == 'charmm':
            self.set_tag(self.text[0: 6])
            self.set_atomNum(self.text[ 6:11])
            self.set_atomType(self.text[12:16])
            self.set_resName(self.text[17:21])
            self.set_segid(self.text[72:76])
            self.set_resid(self.text[22:26])
            self.set_cart(self.text[30:38],self.text[38:46],self.text[46:54])
            self.set_weight(self.text[55:60])
            self.set_bFactor(self.text[61:66])
            self.set_element()
        elif self._format == 'crd':
            pass
        elif self._format == 'amber':
            pass
        else:
            raise AtomError('parse: unknown _format == %s' % self._format)

# Inherited
    def get_hash(self):
        try:
            return '%d%s%s%d%d' % (self._modelNum0,self._segType0,self._segid0,self._resid0,self._atomNum0)
        except AttributeError:
            return 'index-%d' % self._index

    def get_addr(self):
        return '%d%s%s%d%d' % (self._modelNum,self._segType,self._segid,self._resid,self._atomNum)

# Inherited
#   def czech_tag(self,arg=None):
#   def czech_atomNum(self,arg=None):
#   def czech_atomType(self,arg=None):
#   def czech_resName(self,arg=None):
#   def czech_segid(self,arg=None):
#   def czech_resid(self,arg=None):
#   def czech_weight(self,arg=None):
#   def czech_bFactor(self,arg=None):
#   def czech_mass(self,arg=None):
#   def czech_element(self,arg=None):

# Inherited
#   def set_tag(self,arg):
#   def set_atomNum(self,arg):
#   def set_atomType(self,arg):
#   def set_resName(self,arg):
#   def set_segid(self,arg):
#   def set_resid(self,arg):
#   def set_mass(self,arg):
#   def set_xCart(self,arg):
#   def set_yCart(self,arg):
#   def set_zCart(self,arg):
#   def set_cart(self,arg1,arg2,arg3):
#   def set_weight(self,arg):
#   def set_bFactor(self,arg):
#   def set_element(self,arg):

    def set_segType(self):
        """
        Lookup _resName and set tag to: 'pro', 'nuc', 'good' or 'bad'
        """
        if Etc.isPro(self._resName.upper()):
            self._segType = 'pro'
        elif Etc.isNuc(self._resName.upper()):
            self._segType = 'nuc'
        elif Etc.isGood(self._resName.upper()):
            self._segType = 'good'
        else:
            self._segType = 'bad'

    def set_mass(self):
        """
        Lookup _element and set _mass in AMU
        """
        if self._element:
            BaseAtom.set_mass(self,Etc.atomMass[self._element.upper()])
        else:
            BaseAtom.set_mass(self,0.)

    def Print(self,format):
        if format == 'web':
            taco = '%-6s%5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s\n' % \
                (self._tag,self._atomNum,self._atomType,self._resName,self._segid,self._resid,\
                self._cart[0],self._cart[1],self._cart[2],self._weight,self._bFactor,self._element)
        elif format == 'charmm':
            taco = '%-6s%5i %4s %4s %4i    %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s\n' % \
                (self._tag,self._atomNum,self._atomType,self._resName,self._resid,\
                self._cart[0],self._cart[1],self._cart[2],self._weight,self._bFactor,self._segid)
        elif format == 'debug':
            taco = '%15s -> %15s: %5i %4s%4s %1s%4i    %8.3f%8.3f%8.3f%6.2f%6.2f%12s' % \
                (self.get_hash(),self.get_addr(),self._atomNum,self._atomType,self._resName,self._segid,self._resid,\
                self._cart[0],self._cart[1],self._cart[2],self._weight,self._bFactor,self._element)
        elif format in ['crd','cor','card']:
            taco = '%5i%5i %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4i%10.5f\n' % \
                (self._atomNum,self._resIndex,self._resName,self._atomType,\
                self._cart[0],self._cart[1],self._cart[2],self._segid,self._resid,self._weight)
        elif format in ['xcrd','xcor','xcard']:
            taco = '%10i%10i  %-4s      %-4s    %20.10f%20.10f%20.10f  %-5s    %-4i    %20.10f\n' % \
                (self._atomNum,self._resIndex,self._resName,self._atomType,\
                self._cart[0],self._cart[1],self._cart[2],self._segid,self._resid,self._weight)
        else:
            raise AtomError('Print: unknown format == %s'%self._format)
        return taco.upper()

###################
# Special Methods #
###################

# Inherited
#   def __eq__(self,other):
#   def __ne__(self,other):
#   def __lt__(self,other):
#   def __le__(self,other):
#   def __gt__(self,other):
#   def __ge__(self,other):

    def __repr__(self):
        return self.Print('debug')

    def __str__(self):
        return self.Print(self._format)

    def __hash__(self):
        """
        modelNum > _segType > _segid > _resid > _atomNum
        """
        typeScore = {'pro':1,'dna':2,'rna':3,'nuc':4,'good':5,'bad':6}
        segidScore = dict([(char,i+1) for i,char in enumerate(string.lowercase)])
        return (self._modelNum) * 1E12 + typeScore[self._segType] * 1E11 + segidScore[self._segid] * 1E9 + self._resid * 1E5 + self._atomNum

    def _unhash(self,hash):
        unType = {1:'pro',2:'dna',3:'rna',4:'nuc',5:'good',6:'bad'}
        unSegid = dict([(i+1,char) for i,char in enumerate(string.lowercase)])
        modelNum = int(math.floor(hash/1E12))
        hash = int(math.floor(hash%1E12))
        type = int(math.floor(hash/1E11))
        hash = int(math.floor(hash%1E11))
        segid = int(math.floor(hash/1E9))
        hash = int(math.floor(hash%1E9))
        resid = int(math.floor(hash/1E5))
        hash = int(math.floor(hash%1E5))
        atomNum = int(math.floor(hash))
        print '(modelNum,type,segid,resid,atomNum) = (%d,%s,%s,%4d,%5d)' % (modelNum,unType[type],unSegid[segid],resid,atomNum)

###################
# Private Methods #
###################

    def _save_initialData(self):
        """
        Store initial values for crucial atom metadata.
        """
        BaseAtom._save_initialData(self)
        self._modelNum0 = copy(self._modelNum)
        self._segType0  = copy(self._segType)

    def _fix_element(self):
        print self._index
        if self._segType in ['pro','nuc','rna','dna']:
            self._element = self._atomType.strip()[0]
        else:
            self._element = self._atomType.strip()[:2]

# CHARMM Compliance
    def _compliance_resName(self):
        if self._resName in ['hoh','tip3']:
           self._resName = 'tip3'
        elif self._resName == 'his':
           self._resName = 'hsd'
        elif self._resName in ['a','da']:
           self._resName = 'ade'
        elif self._resName in ['t','dt']:
           self._resName = 'thy'
        elif self._resName in ['c','dc']:
           self._resName = 'cyt'
        elif self._resName in ['g','dg']:
           self._resName = 'gua'
        elif self._resName in ['u','du']:
           self._resName = 'ura'
        elif (self._atomType == 'zn  ' and self._resName == 'zn'):
           self._resName = 'zn2'
        elif (self._atomType == 'na  ' and self._resName == 'na'):
           self._resName = 'sod'
        elif (self._atomType == 'cs  ' and self._resName == 'cs'):
           self._resName = 'ces'
        elif (self._atomType == 'cl  ' and self._resName == 'cl'):
           self._resName = 'cla'
        elif (self._atomType == 'ca  ' and self._resName == 'ca'):
           self._resName = 'cal'
        elif (self._atomType ==  ' k  ' and self._resName ==  'k'):
           self._resName = 'pot'
    
    def _compliance_atomType(self):
        if self._resName in ['hoh','tip3']:
           self._atomType = 'oh2 '
        elif (self._resName == 'ile' and self._atomType == ' cd1'):
           self._atomType = ' cd '
        elif (self._atomType == 'na  ' and self._resName == 'na'):
           self._atomType = 'sod '
        elif (self._atomType == 'cs  ' and self._resName == 'cs'):
           self._atomType = 'ces '
        elif (self._atomType == 'cl  ' and self._resName == 'cl'):
           self._atomType = 'cla '
        elif (self._atomType == 'ca  ' and self._resName == 'ca'):
           self._atomType = 'cal '
        elif (self._atomType ==  ' k  ' and self._resName ==  'k'):
           self._atomType = 'pot '
