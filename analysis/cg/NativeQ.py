#!/usr/bin/env python

from analysis.cg.BaseCGAnalysis import load_correlOutput
from analysis.cg.BaseCGAnalysis import BaseCGAnalysis
from lib.Atom import Atom
from lib.Res import Pro
import lib.Etc as Etc
import cPickle as pickle
import os,numpy
import matplotlib.pyplot as pyplot
import lib.myStats as Stats

class NativeQ(BaseCGAnalysis):
    """
    generic vars are:
        _aaCrdFileName          = string
        _basePath               = string
        _aaCrdFileFormatting    = string
        _outputPath             = string
    native contact vars are:
        contactRad              = float
        _nativeContacts         = dict
        _selection              = string
    correl vars are:
        _nScale                 = float
        _temp                   = int
        _maxTimeSteps           = int
        _skipSteps              = int
        _maxCorrelArray         = int
    plotting vars are:
        _rawData                = float**
        _QofT                   = float*
        _nativeRad              = float
        _mu                     = float
        _sigma                  = float
        xAxisScaleFactor        = float
    """
    def __init__(self,aaCrdFileName):
        BaseCGAnalysis.__init__(self,aaCrdFileName)
        self.contactRad         = 4.5
        # Native Contacts
        self._selection         = 'all'
        self._pickle_nativeContacts('r')
        # Plotting
        self.xAxisScaleFactor   = 1./self._maxCorrelArray   # 1 \mu S

##################
# Public Methods #
##################

# Inherited
#   def rm_pickles(self,nscale=None,temp=None)
#   def set_maxTimeSteps(self,Int)
#   def set_outputPath(self,arg=None)

# General
    def set_vars(self,nscale,temp,selection='all'):
        """
        Public method for setting critical instance variables.
        """
        BaseCGAnalysis.set_vars(self,nscale,temp)
        self._selection = selection

    def set_nativeRad(self,arg1,arg2=-1):
        """
        Public method for setting nativeRad, it can be used one of two ways...

        1 arg:    taco.set_nativeRad(8.01)
        2 args:   taco.set_nativeRad(1,2)

        The former method is used to set the native radius to 8.01 angstrom. The
        latter method is used to set the native radius to one mean + 2 std of
        the normal native radius.
        """
        if arg2 == -1:
            assert arg1 > self.contactRad
            self._nativeRad = arg1
        else:
            assert arg1 > 0
            assert arg2 >= 0
            self._mu = arg1
            self._sigma = arg2
            self._nativeRad = self._mu * 5.4 + self._sigma * .85

    def get_nativeContacts(self):
        """
        Read/write self._nativeContacts
        """
        try:
            return self._nativeContacts[self._selection]
        except KeyError:
            self._nativeContacts[self._selection] = self._calc_nativeContacts()
            self._pickle_nativeContacts('w')
            return self._nativeContacts[self._selection]

    def macro_charmm(self,nscale,temp,bunchSize=1):
        """
        Wrapper method for writing/running CHARMM inputs.
        """
        self.set_vars(nscale,temp)
        print 'writing files...'
        self._write_allCorrelInputs(bunchSize)
        print 'running files...'
        self._run_allCorrelInputs(bunchSize)

# Plotting
    def plot(self,nscale,temp,selection='all',write=None,sigma=2):
        """
        Plots a figure containing native contact fraction as a function of time.
        nscale      = float
        temp        = int
        selection   = [ 'all' | 'a' | ... | 'ab' | ... etc. ] -- segid 'a' etc...
        write       = [ None | 'png' ] -- Untested: [ 'ps' | 'svg' | 'pdf' ]
        sigma       = int -- nativeRad specification in standard deviations
        """
        # Setup
        self.set_vars(nscale,temp,selection)
        self.set_nativeRad(1,sigma)
        self._pickle_QofT()
        # Plotting
        pyplot.plot(map(lambda x: self.xAxisScaleFactor*x,range(len(self._QofT))),self._QofT,'lightblue')
        pyplot.xlabel(r'$time\ (ns)$')
        pyplot.ylabel(r'$fraction\ of\ native\ contacts\ (Q)$')
        pyplot.axis([0,self.xAxisScaleFactor*len(self._QofT),0,1])
        pyplot.title('nScale: %4.2f, temp: %d K' % (self._nScale,self._temp))
        # Show/Write
        pyplot.show()
        if write:
            FileName    = '%s/plot_%4.2f_%s_%4.2f.%s' % (self._basePath,self._nScale,self._selection,self._nativeRad,write)
            Format      = write
            pyplot.savefig(FileName,format=Format)

    def plotTempFig(self,nscale,selection='all',write=False,sigma=2):
        """
        Plots a figure containing a 3x3 grid of native contact fraction as a
        function of time subplots.  Each subplot corresponds to a single temp.
        nscale      = float
        selection   = [ 'all' | 'a' | ... | 'ab' | ... etc. ] -- segid 'a' etc...
        write       = [ None | 'png' ] -- Untested: [ 'ps' | 'svg' | 'pdf' ]
        sigma       = int -- nativeRad specification in standard deviations
        """
        fig = pyplot.figure(facecolor='grey')
        fig.text(0.2,0.92,'nScale: %4.2f, selection: %s' % (nscale,selection))
        xCrd = [.15,.425,.7]*3
        yCrd = [.68]*3+[.4]*3+[.12]*3
        # Subplots
        for i,temp in enumerate(range(275,476,25)):
            self.set_vars(nscale,temp,selection)
            self.set_nativeRad(1,sigma)
            self._pickle_QofT()
            sub = fig.add_subplot(3,3,i+1)
            sub.plot(map(lambda x: self.xAxisScaleFactor*x,range(len(self._QofT))),self._QofT,'lightblue')
            sub.set_ylim(0,1)
            sub.set_xticklabels(('','0.2','0.4','0.6','0.8',''),size='x-small')
            sub.set_yticklabels(('0','0.2','0.4','0.6','0.8','1.0'),size='x-small')
            fig.text(xCrd[i],yCrd[i],r'$(T=%3dK,\mu=%4.2f,\sigma=%4.2f)$' % (self._temp,Stats.mean(self._QofT),Stats.stdDev(self._QofT)),size='xx-small')
        # Show/Write
        pyplot.show()
        if write:
            FileName    = '%s/plotTempFig_%4.2f_%s_%4.2f.%s' % (self._basePath,self._nScale,self._selection,self._nativeRad,write)
            Format      = write
            pyplot.savefig(FileName,format=Format)

    def plotCutoffFig(self,nscale,temp,selection='all',write=False):
        """
        Plots a figure containing a 3x3 grid of native contact fraction as a
        function of time subplots.  Each subplot corresponds to a different 
        nativeRad cutoff.
        nscale      = float
        temp        = int
        selection   = [ 'all' | 'a' | ... | 'ab' | ... etc. ] -- segid 'a' etc...
        write       = [ None | 'png' ] -- Untested: [ 'ps' | 'svg' | 'pdf' ]
        """
        fig = pyplot.figure(facecolor='grey')
        fig.text(0.2,0.92,'nScale: %4.2f, temp: %dK, selection: %s' % (nscale,temp,selection))
        xCrd = [.15,.425,.7]*3
        yCrd = [.68]*3+[.4]*3+[.12]*3
        # Subplots
        for i,sigma in enumerate(map(lambda x: x*.5,range(9))):
            self.set_vars(nscale,temp,selection)
            self.set_nativeRad(1,sigma)
            self._pickle_QofT()
            sub = fig.add_subplot(3,3,i+1)
            sub.plot(map(lambda x: self.xAxisScaleFactor*x,range(len(self._QofT))),self._QofT,'lightblue')
            sub.set_ylim(0,1)
            sub.set_xticklabels(('','0.2','0.4','0.6','0.8',''),size='x-small')
            sub.set_yticklabels(('0','0.2','0.4','0.6','0.8','1.0'),size='x-small')
            fig.text(xCrd[i],yCrd[i],r'$(cutoff=%4.2f\AA,\sigma=%3.2f)$' % (self._nativeRad,self._sigma),size='xx-small')
        # Show/Write
        pyplot.show()
        if write:
            FileName    = '%s/plotCutoffFig_%4.2f_%dK_%s.%s' % (self._basePath,self._nScale,self._temp,self._selection,write)
            Format      = write
            pyplot.savefig(FileName,format=Format)

###################
# Private Methods #
###################

# Inherited
#   def _get_correlInputHeader(self)

# Native Contacts
    def _calc_allNativeContacts(self):
        """
        Parse the all-atom .pdb file specified in the constructor and generate a
        list of native contacts. This is determined from inter-atomic distances and
        self.contactRad
        """
        # (1) Load Atoms
        atomList = ( Atom(line,self._aaCrdFileFormatting) for line in open(self._aaCrdFileName) if line.startswith ('ATOM') )
        atomList = [ atom for atom in atomList if atom.element != 'H' ]
        # (2) Load Residues
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
        # (3) Determine Backbone / Sidechain
        bbDict = {}
        scDict = {}
        for res in resList:
            bbDict[res] = [ atom for atom in resDict[res] if Etc.is_backbone(atom.atomType) ]
            if resDict[res].resName != 'GLY':
                scDict[res] = [ atom for atom in resDict[res] if not Etc.is_backbone(atom.atomType) ]
        # (4) Map All-Atom into CG
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
        # (5) Define more efficient Rij machinery
        Rij = {}
        def get_Rij(atom_i,atom_j):
            i = atom_i.segid + str(atom_i.atomNumber)
            j = atom_j.segid + str(atom_j.atomNumber)
            try:
                return Rij[(i,j)]
            except KeyError:
                Rij[(i,j)] = atom_i.bond_length(atom_j)
            return Rij[(i,j)]
        # (6) Determine Native SC Contacts
        lowerOffAlmostDiagonal = ( (resCode_i,resCode_j) for i,resCode_i in enumerate(resList) for j,resCode_j in enumerate(resList) if j - i > 2 and resDict[resCode_i].resName != 'GLY' and resDict[resCode_j].resName != 'GLY' )
        scContact = []
        for resCode_i,resCode_j in lowerOffAlmostDiagonal:
            contact = False
            try:
                for atom_i in scDict[resCode_i]:
                    for atom_j in scDict[resCode_j]:
                        if get_Rij(atom_i,atom_j) < self.contactRad:
                            scContact.append( (cgscDict[resCode_i],cgscDict[resCode_j]) )
                            raise AssertionError
            except AssertionError: pass
        # (7) Reindex CG.atomNumber
        for i,cg in enumerate(cgList):
            cg.atomNumber = i+1
        # (8) Write .pdb file
        writeTo = open('%s_cg.pdb'%self._aaCrdFileName,'w')
        for line in cgList:
            writeTo.write(line.Print('cgcharmm'))
        writeTo.close()
        ###
        return [ (contact[0].atomNumber,contact[1].atomNumber) for contact in scContact ]

    def _calc_nativeContacts(self):
        """
        Most confusing function ever.  Bonus points if you can follow it all.
        Allows the selection of a subset of the entire protein.  Selection is
        based upon segid.  For example selection ='abc' will return all native
        contacts between residues in segments A and B and C.
        """
        # Default behavior, calculates all the native contacts
        if self._selection == 'all':
            return self._calc_allNativeContacts()
        else:
        # Translate selection into a set, then use this set as a filter.
            selection = self._selection.upper()
            filter = set([ char for char in selection if char in ['A','B','C','D'] ])
        # Load up the pdb file containing the canonical mapping from atomNumber to segid.
            segids = [ line.split()[-1] for line in open('%s_cg.pdb'%self._aaCrdFileName) if line.startswith('ATOM') ]
            def map_segids(atomNumber):
                # This dummy function allows us to 'map' from atomNumber to segid
                segidDict = dict(zip(range(len(segids)),segids))
                return segidDict[atomNumber]
        # Go from a list of tuples with atomNumbers to a list of tuples with segids
            mapped_segids = [ map(map_segids,tuple) for tuple in self._nativeContacts['all'] ]
        # Change them to sets so we can test with subsets
            mapped_segid_sets = [ set(taco) for taco in mapped_segids ]
        # Test if each native contact is a subset of the filter.
            truth_array = []
            for taco in mapped_segid_sets:
                if taco.issubset(filter):
                    truth_array.append(1)
                else:
                    truth_array.append(0)
        # Filter away
            return [ tuple for i,tuple in enumerate(self._nativeContacts['all']) if truth_array[i] ]

    def _pickle_nativeContacts(self,mode='r'):
        """
        Pickle/unpickle self._nativeContacts
        """
        fileName = '%s_nativeContacts.pickle' % self._aaCrdFileName
        if mode == 'r':
            try:
                self._nativeContacts = pickle.load(open(fileName))
                print 'found pickled native contacts...'
            except IOError:
                self._nativeContacts = {}
                print 'computing all native contacts...'
                self._nativeContacts['all'] = self._calc_allNativeContacts()
        elif mode == 'w':
            pickle.dump(self._nativeContacts,open(fileName,'w'))
        else:
            raise TypeError('pickle_nativeContacts: unknown mode %s'%mode)

# Charmm input creation
    def _write_singleCorrelInput(self,index,*contacts):
        """
        Called by the _write_allCorrelInputs() method, this method writes a single
        charmm .inp file.
            index       = int
            *contacts   = tuple(int,int)
        """
        # Make sure all the contacts are tuples
        for contact in contacts: assert type(contact) == tuple
        # Make String
        String  = self._get_correlInputHeader()
        String += '!open files for writing\n'
        for i,contact in enumerate(contacts):
            String += '    open unit 1%02d write card name nat%03d.anl\n' % (i,index*len(contacts)+i)
        String += '\n'
        String += 'traj query unit 10\n'
        String += 'correl maxtimesteps %d maxatom 250 maxseries %d\n' % (self._maxTimeSteps,len(contacts))
        for i,contact in enumerate(contacts):
            String += 'enter n%02d bond bynum %d bynum %d geometry\n' % (i,contact[0],contact[1])
        String += 'traj firstu 10 nunit 1 skip %d\n' % self._skipSteps
        String += '\n'
        for i,contact in enumerate(contacts):
            String += 'write n%02d card unit 1%02d\n' % (i,i)
            String += '*Native contact %d: between cgAtoms %d and %d\n' % (index*len(contacts)+i,contact[0],contact[1])
            String += '*\n'
            String += '\n'
        String += 'stop\n'
        # Write file
        write_to = open('%s/nat%03d.inp' % (self._outputPath,index),'w')
        write_to.write(String)
        write_to.close()

    def _write_allCorrelInputs(self,bunchSize=1):
        """
        For a given nScale, temp and atom selection: write all correl inputs. The
        bunchSize parameter allows control of how many correl outputs are produced
        for every correl input.  Larger bunches require more main memory (linearly)
        but speed up correl calls.
        """
        for i,bunch in enumerate(Etc.buncher(self.get_nativeContacts(),bunchSize)):
            self._write_singleCorrelInput(i,*bunch)

    def _run_allCorrelInputs(self,bunchSize=1):
        """
        For a given nScale, temp and atom selection: run all correl inputs. The
        bunchSize parameter allows control of how many correl outputs are produced
        for every correl input.  Larger bunches require more main memory (linearly)
        but speed up correl calls.
        """
        os.chdir(self._outputPath)
        nCorrelInputs = len(list(Etc.buncher(self.get_nativeContacts(),bunchSize)))
        iterator = ( 'nat%03d' % i for i in range(nCorrelInputs) )
        print '  nScale = %4.2f, temp = %d' % (self._nScale,self._temp)
        for i,file in enumerate(iterator):
            print '  %3d of %3d inputs' % (i,nCorrelInputs)
            try:
                os.remove('%s.out'%file)
            except OSError: pass
            os.system('%s < %s.inp > %s.out' % (self.charmmBin,file,file))
        os.chdir(self._basePath)

# Charmm output processing
    def _calc_rawData(self):
        """
        Parse correl output files for nativeContact distance information, and 
        build a matrix of this data.  Finally, transpose this matrix so each
        row corresponds to all data for a single time step.
        """
        def map_selection(tuple):
            allContactsDict = dict(zip(self._nativeContacts['all'],range(len(self._nativeContacts['all']))))
            return allContactsDict[tuple]
        iterator = ( '%s/nat%03d.anl' % (self._outputPath,i) for i in map(map_selection,self.get_nativeContacts()) )
        matrix = numpy.array([ load_correlOutput(file) for file in iterator ])
        print 'transposing data matrix...'
        self._rawData = numpy.transpose(matrix)
        #self._rawData = map(None,*matrix) # transpose the matrix

    def _pickle_rawData(self):
        """
        Pickle/unpickle self._rawData
        """
        fileName = '%s/correlOutputs_%s.pickle' % (self._outputPath,self._selection)
        try:
            self._rawData = pickle.load(open(fileName))
            print 'found pickled correl outputs in %s ...' % fileName
        except IOError:
            print 'processing correl outputs in %s ...' % self._outputPath
            self._calc_rawData()
            pickle.dump(self._rawData,open(fileName,'w'))

# Native contact timeseries building
    def _calc_QofT(self):
        """
        Load self._rawData and compute Q fraction of native contacts as a function
        of time.  self._nativeRad is used in this determination
        """
        print 'getting time series...'
        self._pickle_rawData()
        # nativeRad # Less than this distance, we consider two CG beads to be in native contact
        contactVal = 1./len(self.get_nativeContacts())   # Each contact increases Q by this fraction
        QofT = []                                   # Time series to be filled
        for timeStep in self._rawData:
            fraction = 0.
            for contact in timeStep:
                if contact <= self._nativeRad:
                    fraction += contactVal
            QofT.append(fraction)
        self._QofT = QofT

    def _pickle_QofT(self):
        """
        Pickle/unpickle self._QofT
        """
        fileName = '%s/QofT_%s_%5.2f.pickle' % (self._outputPath,self._selection,self._nativeRad)
        try:
            self._QofT = pickle.load(open(fileName))
            print 'found pickled QofT data in %s ...' % fileName
        except IOError:
            print 'processing QofT data in %s ...' % self._outputPath
            self._calc_QofT()
            pickle.dump(self._QofT,open(fileName,'w'))
