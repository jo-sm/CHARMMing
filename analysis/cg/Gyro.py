#!/usr/bin/env python

from analysis.cg.BaseCGAnalysis import load_correlOutput
from analysis.cg.BaseCGAnalysis import BaseCGAnalysis
import cPickle as pickle
import matplotlib.pyplot as pyplot
import lib.myStats as Stats
import os,math

class Gyro(BaseCGAnalysis):
    """
    generic vars are:
        _aaCrdFileName          = string
        _basePath               = string
        _aaCrdFileFormatting    = string
        _outputPath             = string
    correl vars are:
        _nScale                 = float
        _temp                   = int
        _maxTimeSteps           = int
        _skipSteps              = int
        _maxCorrelArray         = int
    plotting vars are:
        _RgOfT                  = float*
        xAxisScaleFactor        = float
    """
    def __init__(self,aaCrdFileName):
        BaseCGAnalysis.__init__(self,aaCrdFileName)
        # Plotting
        self._selection         = 'all'
        self.xAxisScaleFactor   = 1./self._maxCorrelArray 

##################
# Public Methods #
##################

# Inherited
#   def rm_pickles(self,nscale=None,temp=None)
#   def set_outputPath(self,arg=None)
#   def set_maxTimeSteps(self,Int)

# General
    def set_vars(self,nscale,temp,selection='all'):
        """
        Public method for setting critical instance variables.
        """
        BaseCGAnalysis.set_vars(self,nscale,temp)
        self._selection = selection

    def macro_charmm(self,nscale,temp,selection='all'):
        """
        Wrapper method for writing/running CHARMM inputs.
        """
        self.set_vars(nscale,temp,selection)
        print 'writing file...'
        self._write_singleCorrelInput()
        print 'running file...'
        self._run_SingleCorrelInput()

# Plotting
    def plot(self,nscale,temp,selection='all',write=False):
        """
        Plots a figure containing radius of gyration as a function of time.
        nscale      = float
        temp        = int
        selection   = [ 'all' | 'a' | ... | 'ab' | ... etc. ] segid 'a' etc...
        write       = [ None | 'png' ] -- Untested: [ 'ps' | 'svg' | 'pdf' ]
        """
        # Setup
        self.set_vars(nscale,temp,selection)
        self._pickle_RgOfT()
        # Plotting
        pyplot.plot(map(lambda x: self.xAxisScaleFactor*x,range(len(self._RgOfT))),self._RgOfT,'darkblue')
        pyplot.xlabel(r'$time\ (ns)$')
        pyplot.ylabel(r'$Radius\ of\ Gyration\ (\AA)$')
        pyplot.axis([0,self.xAxisScaleFactor*len(self._RgOfT),math.floor(Stats.min(self._RgOfT)),math.ceil(Stats.max(self._RgOfT))])
        pyplot.title('nScale: %4.2f, temp: %dK, selection: %s' % (self._nScale,self._temp,self._selection))
        # Show/Write
        pyplot.show()
        if write:
            FileName    = '%s/RgPlot_%4.2f_%d_%s.%s' % (self._basePath,self._nScale,self._temp,self._selection,write)
            Format      = write
            pyplot.savefig(FileName,format=Format)

    def plotTempFig(self,nscale,temps,selection='all',start=0.1,stop=1,write=False):
        """
        Plots average Radius of Gyration as a function of temperature.
        nscale      = float
        temps       = tuple(int)
        selection   = [ 'all' | 'a' | ... | 'ab' | ... etc. ] segid 'a' etc...
        start       = 0 <= float <= 1 -- Start averaging from this fraction of the trajectory
        stop        = 0 <= start < float <= 1 -- Stop averaging from this faction of the trajectory
        write       = [ None | 'png' ] -- Untested: [ 'ps' | 'svg' | 'pdf' ]
        """
        assert 0. <= float(start) <= 1.
        assert 0. <= float(start) < float(stop) <= 1.
        # Setup
        temps = map(int,temps)
        meanRgs = []
        for temp in temps:
            self.set_vars(nscale,temp,selection)
            self._pickle_RgOfT()
            startFrame = int(len(self._RgOfT)*float(start))
            stopFrame  = int(len(self._RgOfT)*float(stop))
            meanRgs.append(Stats.mean(self._RgOfT[startFrame:stopFrame]))
        # Plotting
        pyplot.plot(temps,meanRgs,'g^')
        pyplot.xlabel(r'$temperature\ (K)$')
        pyplot.ylabel(r'$Radius\ of\ Gyration\ (\AA)$')
        pyplot.title('nScale: %4.2f, selection: %s' % (self._nScale,self._selection))
        # Show/Write
        pyplot.show()
        if write:
            FileName    = '%s/RgTempPlot_%4.2f_%s_%3.1f_%3.1f.%s' % (self._basePath,self._nScale,self._selection,start,stop,write)
            Format      = write
            pyplot.savefig(FileName,format=Format)

###################
# Private Methods #
###################

# Inherited
#   def _get_correlInputHeader(self)

# Charmm input creation
    def _write_singleCorrelInput(self):
        """
        Called by macro_charmm() method, this method writes a single charmm .inp
        file for a radius of gyration calculation.
        """
        String  = self._get_correlInputHeader()
        String += '!open files for writing\n'
        String += '    open unit 100 write card name gyro_%s.anl\n' % self._selection
        String += '\n'
        if self._selection == 'all':
            String += 'defi taco select all end\n'
        else:
            String += 'defi taco select '
            for i,letter in enumerate(self._selection):
                String += 'segid %s ' % letter.upper()
                if i != len(self._selection)-1: String += '.or. '
            String += 'end\n'
        String += '\n'
        String += 'traj query unit 10\n'
        String += 'correl maxtimesteps %d maxatom 250 maxseries 1\n' % self._maxTimeSteps
        String += 'enter gyro gyration\n'
        String += 'traj firstu 10 nunit 1 skip %d taco\n' % self._skipSteps
        String += '\n'
        String += 'write gyro card unit 100 \n'
        String += '* Radius of Gyration - Selection: %s\n' % self._selection
        String += '* nScale: %4.2f\n' % self._nScale
        String += '* temp: %dK\n' % self._temp
        String += '*\n'
        String += '\n'
        String += 'stop\n'
        # Write file
        write_to = open('%s/gyro_%s.inp' % (self._outputPath,self._selection),'w')
        write_to.write(String)
        write_to.close()

    def _run_SingleCorrelInput(self):
        """
        Called by macro_charmm() method, this method runs a single charmm .inp
        file for a radius of gyration calculation.
        """
        os.chdir(self._outputPath)
        try:
            os.remove('gyro_%s.out' % self._selection)
        except OSError: pass
        os.system('%s < gyro_%s.inp > gyro_%s.out' % (self._charmmBin,self._selection,self._selection))
        os.chdir(self._basePath)

    def _pickle_RgOfT(self):
        """
        Wrapper method for Radius of Gyration as a function of time data.
        """
        fileName = '%s/RgOfT_%s.pickle' % (self._outputPath,self._selection)
        try:
            self._RgOfT = pickle.load(open(fileName))
            print 'found pickled RgOfT data in %s ...' % fileName
        except IOError:
            print 'processing RgOfT data in %s ...' % self._outputPath
            self._RgOfT = load_correlOutput('%s/gyro_%s.anl' % (self._outputPath,self._selection))
            pickle.dump(self._RgOfT,open(fileName,'w'))

