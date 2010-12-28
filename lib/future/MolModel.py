#!/usr/bin/env python

from lib.future.BaseAtom import BaseAtom
from copy import deepcopy

class MolModelError(Exception):
    """
    Exception to raise when errors occur involving the Mol class."""
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class MolModel(list):
    """
        self._modelNum  = int
    """
    def __init__(self,atomList=[],modelNum=None):
        list.__init__(self,atomList)
        if modelNum:
            self._modelNum  = modelNum
        else:
            if atomList:
                self._modelNum      = self[0]._modelNum
            else:
                self._modelNum      = 0
        for atom in self: atom._modelNum = self._modelNum

    def czech_Atom(self):
        """
        Raises an error if all of the objects inside the Seg object are not Atom objects.
        """
        try:
            for i,atom in enumerate(self):
                assert issubclass(atom.__class__,BaseAtom)
        except AssertionError:
            raise MolError('czech_Atom: %dth item in the MolModel is not of type: Atom' % i)

    def append(self,atom):
        try:
            assert issubclass(atom.__class__,BaseAtom)
        except AssertionError:
            print 'append: object appended was not derived from type: BaseAtom'
        atom._modelNum = self._modelNum
        list.append(self,atom)

    def extend(self,iterable):
        for atom in iterable:
            self.append(atom)

    def __add__(self,other):
        buffer = deepcopy(self)
        for atom in other:
            buffer.append(atom)
        return buffer
