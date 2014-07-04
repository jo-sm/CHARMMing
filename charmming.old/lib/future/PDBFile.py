#!/usr/bin/env python

import string,copy,os.path
from lib.future.MolModel import MolModel
from lib.future.NewAtom import NewAtom
from lib.future.Seg import Seg
import lib.Etc as Etc
import lib.myStats as Stats
#import lib.Etc as Etc
#import lib.Atom as Atom
#from lib.Seg import Seg
#import lib.myStats as Stats

ALPHAnum = string.uppercase + string.digits

def get_formatting(input):
    """
    Takes a pdbFile or an iterator of strings, and returns a string indicating if it is 'web', 'charmm'
    or 'unknown' format.
    """
    try:
        iterator = ( line.lower() for line in open(input) if line.startswith(('atom','hetatm')) )
    except IOError:
        iterator = ( line.lower() for line in input if line.startswith(('atom','hetatm')) )
    for line in iterator:
        if line[21:22] == ' ' and line[72:73] in ALPHAnum:
            return 'charmm'
        elif ( line[21:22] in ALPHAnum ) and ( line[12:14].strip() == line[66:].strip() ):
            return 'web'
    return 'unknown'

class PDBFileError(Exception):
    """
    Exception to raise when the parser encounters _really_ big problems.
    """
    def __init__(self,value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class PDBFile(object):
    """
    #
        _pdbFileName        = string
        _pdbCode            = string
        _basePath
        _warnings           = string*
        _debug              = bool      [ False | True ]
    # Text Processing
        _format             = string    [ 'web' | 'charmm' ]
        _text               = string*
        _headerText         = string*
        _crdText            = string*
        _footerText         = string*
    #
        _models             = MolModel*
        _segments           = Dict
    """
    _defaultFormatting = 'web'
    def __init__(self,pdbFileName,charmm=True,debug=False):
        #
        self._pdbFileName   = Etc.expandPath(pdbFileName)
        self._pdbCode       = os.path.basename(self._pdbFileName).split('.pdb')[0]
        self._basePath      = os.path.dirname(self._pdbFileName)
        self._warnings      = []
        self._debug         = debug
        # Text Processing
        self._formatting()              # self._format
        self._text = [ line.lower() for line in open(pdbFileName) ]
        self._partition_text()
        self._fix_segids()
        # Date Structures
        self._build_models()            # self._models = []
        self._build_segments()          # self._segments = {}
        # Cleaning
        self._fix_multi_model()
        self._diff_nucThyUra()
        self._diff_nucNames()
        self.sort()
        self.reindex()
        if charmm:
            self._compliance_charmmTerminalOxygen()
            self._compliance_resName()
            self._compliance_atomType()

##################
# Public Methods #
##################

# General
    def get_segKeys(self,*args):
        keys = map(set,self._segments.keys())
        if not args:
            pass
        else:
            args = set(Etc.flatten(args))
            keys = filter(args.issubset,keys)
        keys = map(tuple,keys)
        def order(Tuple):
            buffer = [None,None,None]
            for item in Tuple:
                if type(item) == int:
                    buffer[0] = item
                elif item in string.lowercase:
                    buffer[2] = item
                else:
                    buffer[1] = item
            return tuple(buffer)
        keys = map(order,keys)
        def _hash(Tuple):
            typeScore   = {'pro':1,'dna':2,'rna':3,'nuc':4,'good':5,'bad':6}
            segidScore  = dict( [ (char,i+1) for i,char in enumerate(string.lowercase) ] )
            return Tuple[0] * 1000 + typeScore[Tuple[1]] * 100 + segidScore[Tuple[2]] * 1
        keys.sort(key=_hash)
        return keys
    
    def get_segValues(self,*args):
        """
        Wrapper method for self._segments
        args should be 0 or more of the following:
            _modelNum            = int
            _segid               = char
            _segType             = string    [ 'pro' | 'nuc' | 'good' | 'bad' | 'rna' | 'dna' ]
        examples:
            taco.get_segs()                 -> returns all segments
            taco.get_segs('a')              -> returns all segments whose _segid == 'a'
            taco.get_segs(0,'d','pro')      -> returns all protein segments in model 0 whose _segid == 'd'
        """
        return [ self._segments[key] for key in self.get_segKeys(args) ]

    def sort(self,*args):
        """
        Sort by scoring each atom line using the following priority
        _modelNum > _segType > _segid > _resid > _atomNum
        """
        for model in self._models:
            model.sort()
        for segment in self._segments.itervalues():
            segment.sort()

    def reindex(self):
    # atomNumber
        for segment in self.get_segValues():
            for i,atom in enumerate(segment):
                atom._atomNum = i+1
    # resid
        for segment in self.get_segValues():
            oldResid = segment[0]._resid
            newResid = 1
            for atom in segment:
                if not atom._resid == oldResid:
                    oldResid = atom._resid
                    newResid += 1
                atom._resid = newResid
    # resIndex
        for model in self._models:
            model[0].resIndex   = 1
            this_addr           = model[0].get_addr()
            new_resIndex        = 1
            for atom in model:
                if not atom.get_addr() == this_addr:
                    this_addr = atom.get_addr()
                    new_resIndex += 1
                atom.resIndex = new_resIndex

    def WritePDB(self,style='charmming',model=0):
        typeMap = {'nuc':'nuc','pro':'pro','good':'goodhet','bad':'het','dna':'dna','rna':'rna'}
        tagMap  = {'nuc':'ATOM','pro':'ATOM','good':'ATOM','bad':'HETATM','dna':'ATOM','rna':'ATOM'}
        to_stdOut = []
        if style == 'charmming':
            for segment in self.get_segValues(model):
                # File Writing
                map(lambda x:tagMap[x._segType],segment)
                outputName = 'new_%s-%s-%s-%d.pdb' % (self._pdbCode,segment._segid,typeMap[segment._segType],segment._modelNum)
                outputName = outputName.lower()
                to_stdOut.append( '%s-%s' % (segment._segid,typeMap[segment._segType]) )
                write_to = open(outputName,'w')
                outputString = ''.join([ atom.Print('charmm') for atom in segment ])
                write_to.write(outputString)
                write_to.write('TER\n')
                write_to.close()
        elif style == 'generic':
            outputName = '%s_%d.pdb' % (self._pdbCode,segment._modelNum)
            outputName = outputName.lower()
            write_to = open(outputName,'w')
            outputString = ''.join([ atom.Print('charmm') for atom in self._models[model] ])
            write_to.write(outputString)
            write_to.write('TER\n')
            write_to.close()
        # stdout
        print 'natom=%d' % len(self._models[model])
        print 'nwarn=%d' % len(self._warnings)
        if self._warnings:
            print 'warnings=', str(self._warnings)
        print 'seg=', str(to_stdOut)

###################
# Private Methods #
###################

# Constructor
    def _formatting(self):
        """
        Wrapper method to detect PDB text formatting, if it can't tell it sets _format = 'web'
        """
        self._format = get_formatting(self._pdbFileName)
        if self._format == 'unknown':
            self._format = PDBFile._defaultFormatting
            if self._debug: print 'PDB formatting not detected, continuing with %s formatting...' % self._format

    def _partition_text(self):
        """
        Partitions the PDB file text into: _headerText, _footerText and crdText
        """
        atomBool = []
        for line in self._text:
            if line.startswith(('atom','hetatm')):
                atomBool.append(1)
            else:
                atomBool.append(0)
        firstAtom = atomBool.index(1)
        atomBool.reverse()
        lastAtom = atomBool.index(1)
        self._headerText    = self._text[:firstAtom]
        self._crdText       = self._text[firstAtom:-lastAtom]
        self._footerText    = self._text[-lastAtom:]

    def _fix_segids(self):
        # Detect missing seg info
        models = Etc.paragraphs(self._crdText,['model'])
        for model in models:
            badSegs = False
            if self._format == 'web':
                testMe = set([ line[21:22] for line in self._crdText if line.startswith(('atom','hetatm')) ])
                if ' ' in testMe:
                    badSegs = True
            elif self._format == 'charmm':
                testMe = set([ line[72:76].strip() for line in self._crdText if line.startswith(('atom','hetatm')) ])
                if '' in testMe:
                    badSegs = True
        # Regenerate missing seg info
        if badSegs:
            print '%s: One or more segids appear to be missing, attempting to fix them...' % self._pdbCode
            nTH_seg = 0
            if self._format == 'web':
                repPlace = 21
            elif self._format == 'charmm':
                repPlace = 72
            for i,line in enumerate(self._crdText):
                if line.startswith(('atom','hetatm')):
                    self._crdText[i] = Etc.replaceChar(line,repPlace,string.uppercase[nTH_seg])
                elif line.startswith('ter'):
                    nTH_seg += 1
                elif line.startswith('model'):
                    nTH_seg = 0

    def _build_models(self):
        self._models    = []
        models = Etc.paragraphs(self._crdText,['model'])
        for i,model in enumerate(models):
            iterator = ( NewAtom(line,self._format,y,i) for y,line in enumerate(model) if line.startswith(('atom','hetatm')) )
            self._models.append(MolModel(iterator,i))

    def _build_segments(self):
        self._segments  = {}
        for model in self._models:
            for atom in model:
                try:
                    self._segments[(atom._modelNum,atom._segType,atom._segid)].append(atom)
                except KeyError:
                    self._segments[(atom._modelNum,atom._segType,atom._segid)] = Seg([atom])

    def _fix_multi_model(self):
        """
        Tag redundant atoms from each multi-model section for removal.
          (a) First check if a line belongs to a multi-model section
          (b) Dump the current line to the buffer if we:
                > Start a new MM section (empty buffer)
                > Stay in a MM section (leading character is 'bigger' ie, 'BMET' > 'AMET')
          (c) Empty the buffer and clean up if we leave the current MM section:
                > The leading character is 'smaller' ie, 'AMET' > 'BMET'
                >
                Now we have something to compare the following lines in the multi-
                model section to (the buffer is filled)
          (d) We need to determine if we are still in the current multi-model
                section, and when the section is terminated
          (e) Check the weight of the current line vs. the weight of the buffered
                line
          (f) If the buffered line has a bigger weight than the current line, tag
                the current line for removal
          (g) If the buffered line has a smaller weight than the current line, tag
                the buffered line for removal and buffer the current line
          (h) We are in the next multi-model section, fill the buffer with the current line.
        Multi-Model clean up
          (i) Delete atoms tagged for removal by the previous routine
          (j) Relabel residue names from nRES -> RES (remove multi-model prefixes)
        """
        for key,segment in self._segments.iteritems():
            Buffer = []
            keep = False
            for i,atom in enumerate(segment):
                try:
                    if len(atom._resName) == 4 and ( (i == 0 and segment[0]._resName[-3:] == segment[1]._resName[-3:]) or (i != 0 and segment[i]._resName[-3:] == segment[i-1]._resName[-3:]) ):
                        if not Buffer:
                            Buffer.append(atom)
                            continue
                        if atom._resName[0] > Buffer[-1]._resName[0]:
                            Buffer.append(atom)
                        else:
                            AuxBuffer = atom
                            raise AssertionError
                    else:
                        if Buffer:
                            AuxBuffer = ''
                            raise AssertionError
                except AssertionError:
                    maxWeight = Stats.max( [ atom1._weight for atom1 in Buffer ] )
                    for atom2 in Buffer:
                        if keep == False:
                            if atom2._weight == maxWeight:
                                atom2.note = 'rename'
                                keep = True
                            else:
                                atom2.note = 'remove'
                        else:
                            atom2.note = 'remove'
                    if AuxBuffer:
                        Buffer = [AuxBuffer]
                    else:
                        Buffer = []
                    keep = False
            # Clean up
            self._segments[key] = Seg(( atom for atom in self._segments[key] if atom.note != 'remove' ))
            for atom in segment:
                if atom.note == 'rename':
                    atom._resName = atom._resName[-3:]
        for i,model in enumerate(self._models):
            self._models[i] = MolModel(( atom for atom in model if atom.note != 'remove' ))

    def _diff_nucThyUra(self):
        """
        Differentiate dna/rna using Thymine/Uracil
        """
        for segment in self.get_segValues('nuc'):
            thymine = False
            uracil  = False
            for atom in segment:
                if atom._resName in ['T','THY','DT']:
                    thymine = True
                elif atom._resName in ['U','URA','DU']:
                    uracil  = True
        #WARN# Throw a warning when Uracil and Thymine are found in the same segment
            if thymine and uracil:
                self._warnings.append('_diff_nucThyUra: (%d,%s,%s) URA & THY in same segment' % (segment._modelNum,segment._segType,segment._segid))
            elif thymine:
                for atom in segment:
                    atom._segType = 'dna'
            elif uracil:
                for atom in segment:
                    atom._segType = 'rna'
        self._fix_nuc()

    def _diff_nucNames(self):
        """
        Differentiate dna/rna using pdb residue names
        """
        for segment in self.get_segValues('nuc'):
            ribo = False
            deoxy = False
            for atom in segment:
                if atom._resName in ['A','C','G','U','I']:
                    ribo = True
                elif atom._resName in ['DA','DC','DG','DT','DI']:
                    deoxy = True
        #WARN# Throw a warning when ribo- and deoxy- are found in the same segment
            if ribo and deoxy:
                self._warnings.append('_diff_nucNames: (%d,%s,%s) ribo- & deoxy- in same segment' % (segment._modelNum,segment._segType,segment._segid))
            elif ribo:
                for atom in segment:
                    atom._segType = 'rna'
            elif deoxy:
                for atom in segment:
                    atom._segType = 'dna'
        self._fix_nuc()

    def _fix_nuc(self):
        for segment in self.get_segValues('nuc'):
            nucType = list(set([atom._segType for atom in segment]))
            if len(nucType) == 1:
                if nucType[0] == 'rna':
                    segment._segType = 'rna'
                    self._segments[(segment._modelNum,'rna',segment._segid)] = segment
                    del self._segments[(segment._modelNum,'nuc',segment._segid)]
                elif nucType[0] == 'dna':
                    segment._segType = 'dna'
                    self._segments[(segment._modelNum,'dna',segment._segid)] = segment
                    del self._segments[(segment._modelNum,'nuc',segment._segid)]
                elif nucType[0] == 'nuc':
        #WARN# Throw a warning when a nucleic acid segment isn't identified as DNA or RNA
                    self._warnings.append('_fix_nuc: (%d,%s,%s) Parser cant ID nuc segment' % (segment._modelNum,segment._segType,segment._segid))
                else:
                    raise PDBParserError('_fix_nuc: (%d,%s,%s) unexpected segment _type: %s encountered' % (segment._modelNum,segment._segType,segment._segid))
            else:
        #WARN# Throw a warning when ribo- and deoxy- are found in the same segment
                    self._warnings.append('_fix_nuc: (%d,%s,%s) ribo- & deoxy- in same segment' % (segment._modelNum,segment._segType,segment._segid))

# CHARMM compliance
    def _compliance_charmmTerminalOxygen(self):
        for segment in self.get_segValues('pro'):
            last_resid = segment[-1]._resid
            for atom in segment:
                if atom._resid == last_resid and atom._atomType == ' o  ':
                    atom._atomType = ' ot1'
                elif atom._resid == last_resid and atom._atomType == ' oxt':
                    atom._atomType = ' ot2'

    def _compliance_resName(self):
        for model in self._models:
            for atom in model:
                atom._compliance_resName()

    def _compliance_atomType(self):
        for model in self._models:
            for atom in model:
                atom._compliance_atomType()
