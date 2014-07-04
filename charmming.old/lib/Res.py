#!/usr/bin/env python
# Filename Res.py

import numpy,copy
import Etc

class ResError(Exception):
    """Exception to raise when errors occur involving the Res class."""
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class ProError(Exception):
    """Exception to raise when errors occur involving the Pro class."""
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Res(list):
    """The Res class is derived from the __builtin__.list class.  Additional
    methods and metadata are provided to assist with the manipulation of residues."""
    def __init__(self,atom_list=[]):
        """The constructor is given a list of Atom objects and then builds
        a Res container object."""
        list.__init__(self,atom_list)       # A list containing the atom list
        # Initialization
        self.chk_resName()
        self.chk_resid()
        self.chk_segid()
        self.chk_tag()
        self.chk_type()
        # Variables
        self.resName    = self[0].resName
        self.resid      = self[0].resid
        self.segid      = self[0].segid
        self.tag        = self[0].tag
        self.type       = self[0].type
    
    def chk_resName(self):
        """Verifies all the members of the residue have the same resName.
        Raises 'ResError' otherwise."""
        buffer = [ atom.resName[-3:] for atom in self ]
        chk_len = len(set(buffer))
        if chk_len != 1:
            raise ResError('chk_resName: Atom objects from %i residues are being fed into this Res object'%chk_len)
        return

    def chk_resid(self):
        """Verifies all the members of the residue have the same resid. Raises
        'ResError" otherwise."""
        buffer = [ atom.resid for atom in self ]
        chk_len = len(set(buffer))
        if chk_len != 1:
            raise ResError('chk_resid: Atom objects from %i residues are being fed into this Res object'%chk_len)
        return

    def chk_segid(self):
        """Verifies all the members of the residue have the same segid. Raises
        'ResError" otherwise."""
        buffer = [ atom.segid for atom in self ]
        chk_len = len(set(buffer))
        if chk_len != 1:
            raise ResError('chk_segid: Atom objects from %i segments are being fed into this Res object'%chk_len)
        return

    def chk_tag(self):
        """Verifies all the members of the residue have the same tag."""
        buffer = [ atom.tag for atom in self ]
        chk_len = len(set(buffer))
        if chk_len != 1:
            raise ResError('chk_type: Atom objects with i% tags are being fed into this Res object'%chk_len)
        return

    def chk_type(self):
        """Verifies all the members of the residue have the same type."""
        buffer = [ atom.type for atom in self ]
        chk_len = len(set(buffer))
        if chk_len != 1:
            raise ResError('chk_type: Atom objects of i% types are being fed into this Res object'%chk_len)
        return

    def get_mass(self):
        """Return the residue's mass in AMU."""
        sum = 0
        for atom in self:
            sum += atom.mass
        return sum

    def get_com(self):
        """For a given residue, calculate the cartesian coordinates of the
        center of mass of _all_ of the atoms (including hydrogen), and return 
        them as a numpy array [x,y,z]."""
        mass_times_coord_sum    = numpy.array([0.,0.,0.])
        mass_sum                = numpy.array([0.,0.,0.])
        for atom in self:
            mass_times_coord_sum    += atom.mass * atom.cart
            mass_sum                += atom.mass
        return mass_times_coord_sum/mass_sum

    def get_heavy_com(self):
        """For a given residue, calculate the cartesian coordinates of the
        center of mass of the heavy atoms, and return them as a numpy
        array [x,y,z]."""
        mass_times_coord_sum    = numpy.array([0.,0.,0.])
        mass_sum                = numpy.array([0.,0.,0.])
        for atom in self:
            if   atom.element == 'H':               # Ignore hydrogen atoms
                pass
            else:
                mass_times_coord_sum    += atom.mass * atom.cart
                mass_sum                += atom.mass
        return mass_times_coord_sum/mass_sum

class Pro(Res):
    """The Pro class is derived from the more general Res class.  Additional protein
    specific methods and metadata are provided to assist with the manipulation of protein
    residues."""
    def __init__(self,atom_list=[]):
        """The constructor is given a list of pdb_atom objects and then builds
        a Pro container object."""
        Res.__init__(self,atom_list)
        # Initialization
        self.chk_pro()

    def chk_pro(self):
        """Verifies Atom's of the protein type are being fed into the Pro object."""
        if not self[0].type == 'pro':
            raise ProError('chk_pro: Non protein Atom objects are being fed into this Pro object.')
        return

    def get_alpha_carbon(self):
        """For a given protein residue, return the alpha carbon pdb_atom object."""
        for atom in self:
            if atom.atomType == ' CA ':
                backbone               = copy.deepcopy(atom)
                backbone.tag           = 'ATOM'
                backbone.atom_number   = backbone.resid
                backbone.atomType      = '   B'
                backbone.segid         = backbone.segid[0] 
                backbone.weight        = 0
                backbone.b_factor      = 0
                backbone.mass          = Etc.get_aa_mass('GLY')
                return backbone
        raise NoAlphaCarbonError

    def get_sc_com(self):
        """For a given protein residue, calculate the cartesian coordinates of 
        the center of mass of the heavy atoms in the sidechain.  If the residue 
        is Glycine, a None will be returned."""
        if self.resName == 'GLY':                  # Catch Glycine (which has no heavy atoms in its SC)
            return None
        mass_times_coord_sum    = numpy.array([0.,0.,0.])
        mass_sum                = numpy.array([0.,0.,0.])
        for atom in self:
            if   atom.element == 'H':                   # Ignore hydrogen atoms
                pass
            elif Etc.is_backbone( atom.atomType ):   # Ignore Backbone atoms
                pass
            else:
                mass_times_coord_sum    += atom.mass * atom.cart
                mass_sum                += atom.mass
        return mass_times_coord_sum/mass_sum

    def get_sc(self):
        """For a given protein residue, return a pdb_atom object with coordinates
        at the side chain center of mass."""
        if self.resName == 'GLY':                  # Catch Glycine (which has no heavy atoms in its SC)
            return None
        side_chain              = copy.deepcopy(self[0])
        # Variables of returned pdb_atom 
        side_chain.tag          = 'ATOM'
        side_chain.atom_number  = self.resid
        side_chain.atomType     = '   S'
        side_chain.segid          = side_chain.segid[0] 
        side_chain.cart         = self.get_sc_com()
        side_chain.weight       = 0
        side_chain.b_factor     = 0
        side_chain.mass         = Etc.get_aa_mass(side_chain.resName) - Etc.get_aa_mass('GLY')
        return side_chain

