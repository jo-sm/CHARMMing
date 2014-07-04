#!/usr/bin/env python

from lib.future.BaseAtom import BaseAtom

class SegError(Exception):
    """
    Exception to raise when errors occur involving the Seg class."""
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Seg(list):
    """
        self._modelNum  = int
        self._segType   = string    # [ 'pro' | 'good' | 'bad' | 'nuc' | 'rna' | 'dna' ]
        self._segid     = char
    """
    def __init__(self,atomList=[]):
        list.__init__(self,atomList)
        self.czech_Atom()
        self._modelNum  = -1
        self._segType   = 'unk'
        self._segid     = '?'
        if atomList:
            self._modelNum  = self[0]._modelNum
            self._segType   = self[0]._segType
            self._segid     = self[0]._segid
        self.czech_allAttr()

    def czech_Atom(self):
        """
        Raises an error if all of the objects inside the Seg object are not Atom objects.
        """
        try:
            for i,atom in enumerate(self):
                assert issubclass(atom.__class__,BaseAtom)
        except AssertionError:
            raise SegError('czech_Atom: %dth item in the Segment is not of _segType: Atom' % i)

    def czech_attr(self,attr):
        """
        Raises an error if specified attr is not the same for all Atom objects inside of the Seg.
        """
        attrSet = set([ getattr(atom,attr) for atom in self ])
        chk_len = len(attrSet)
        if chk_len > 1:
            raise SegError('czech_%s: %d different Atom attributes found in Segment' % (attr,chk_len))

    def czech_allAttr(self):
        self.czech_attr('_modelNum')
        self.czech_attr('_segType')
        self.czech_attr('_segid')

    def append(self,atom):
        try:
            assert issubclass(atom.__class__,BaseAtom)
        except AssertionError:
            raise SegError('append: object appended was not derived from _segType: BaseAtom')
        atom._modelNum  = self._modelNum
        atom._segType   = self._segType
        atom._segid     = self._segid
        list.append(self,atom)

    def extend(self,iterable):
        for atom in iterable:
            self.append(atom)

    def __add__(self,other):
        buffer = deepcopy(self)
        for atom in other:
            buffer.append(atom)
        return buffer
