#!/usr/bin/env python

import os
import lib.Etc as Etc

def load_correlOutput(fileName):
    """
    Reads a charmm formatted correl output reporting a bondlength as a function
    of time.  Data is returned in a list.
    """
    assert fileName.endswith('anl')
    buffer = []
    for line in open(fileName):
        try:
            dummy = int(line.split()[0])
            buffer.append(float(line.split()[1]))
        except ValueError: pass
    return buffer

class BaseCGAnalysis(object):
    """
    Base class from which CG Analysis scripts inherit. Common functionalities
    such as trajectory processing and directory structure are described herein.
    """
    def __init__(self,aaCrdFileName):
        self._charmmBin         = 'cg_charmm'
        self._aaCrdFileName     = Etc.expandPath(aaCrdFileName)
        self._basePath          = os.path.dirname(Etc.expandPath(self._aaCrdFileName))
        self._aaCrdFileFormatting = Etc.get_ver(self._aaCrdFileName)
        self._outputPath        = None
        # Correl Vars
        self._nScale            = None
        self._temp              = None
        self._maxCorrelArray    = 10000
        self.set_maxTimeSteps(200000000)

##################
# Public Methods #
##################

    def rm_pickles(self,nscale=None,temp=None):
        """
        Recursively searches from the active outputPath for pickled data and
        deletes it.  Defaults to all data for a given AA .pdb file.
        """
        # Set Path
        if nscale and temp:
            self.set_outputPath()
            path = self._outputPath
        if nscale and not temp:
            path = '%s/nscale_%d' % (self._basePath,nscale)
        else:
            path = self._basePath
        # Purgation
        pickleFiles = ( os.path.abspath('/'.join((root,file))) for root,dirs,files in os.walk(path) for file in files if file.endswith('pickle') )
        for pickleFile in pickleFiles:
            print 'removing %s' % pickleFile
            os.remove(pickleFile)

    def set_vars(self,nscale,temp):
        """
        Public method for setting self._nScale and self._temp variables.
        """
        self._nScale = nscale
        self._temp = temp
        self.set_outputPath()

    def set_maxTimeSteps(self,Int):
        """
        Public method for setting length of trajectory files to process.
        """
        Int=int(Int)
        self._maxTimeSteps  = Int
        self._skipSteps     = self._maxTimeSteps / self._maxCorrelArray

    def set_outputPath(self,arg=None):
        """
        Public method for setting the active directory that will be read and
        written to.  There are two ways to use this method...
            0 args: taco.set_outputPath()
            1 arg:  taco.set_outputPath(PATH)
        The former syntax is the default and sets self._outputPath to the following:
        '%s/nscale_%4.2f/%dk' % (self._basePath,self._nScale,self._temp)
        The latter syntax simply sets self._outPath to PATH.
        """
        if arg: # Explicit Path Assignment
            arg = Etc.expandPath(arg)
            try:
                os.mkdir(arg)
            except OSError: pass
            self._outputPath = arg
        else: # Default
            try:
                os.mkdir('%s/nscale_%4.2f' % (self._basePath,self._nScale))
            except OSError: pass
            try:
                os.mkdir('%s/nscale_%4.2f/%dk' % (self._basePath,self._nScale,self._temp))
            except OSError: pass
            self._outputPath = '%s/nscale_%4.2f/%dk' % (self._basePath,self._nScale,self._temp)

###################
# Private Methods #
###################

    def _get_correlInputHeader(self):
        String  = '* A first attempt at native contact scripting\n'
        String += '*\n'
        String += 'bomlev -1\n'
        String += 'wrnlev 5\n'
        String += '\n'
        String += '! toppar\n'
        String += '    read rtf  card name %s/ld-nscale%4.2f/merged.rtf\n' % (self._basePath,self._nScale)
#       String += '    read rtf  card name /v/bigbox12/home/tim/projects/cg_model/1prb/ld-nscale%4.2f/merged.rtf\n' % self._nScale
#       String += '    read rtf  card name /v/bigbox12/home/tim/projects/cg_model/2qmt/ld-nscale%4.2f/merged.rtf\n' % self._nScale
        String += '    read para card name %s/ld-nscale%4.2f/merged.prm\n' % (self._basePath,self._nScale)
#       String += '    read para card name /v/bigbox12/home/tim/projects/cg_model/1prb/ld-nscale%4.2f/merged.prm\n' % self._nScale
#       String += '    read para card name /v/bigbox12/home/tim/projects/cg_model/2qmt/ld-nscale%4.2f/merged.prm\n' % self._nScale
        String += '\n'
        String += '! psfcor\n'
        String += '    read psf  card name %s/ld-nscale%4.2f/cg-vacuum.psf\n' % (self._basePath,self._nScale)
#       String += '    read psf  card name /v/bigbox12/home/tim/projects/cg_model/1prb/ld-nscale%4.2f/cg-vacuum.psf\n' % self._nScale
#       String += '    read psf  card name /v/bigbox12/home/tim/projects/cg_model/2qmt/ld-nscale%4.2f/cg-vacuum.psf\n' % self._nScale
        String += '    read coor card name %s/ld-nscale%4.2f/cg-vacuum.crd\n' % (self._basePath,self._nScale)
#       String += '    read coor card name /v/bigbox12/home/tim/projects/cg_model/1prb/ld-nscale%4.2f/cg-vacuum.crd\n' % self._nScale
#       String += '    read coor card name /v/bigbox12/home/tim/projects/cg_model/2qmt/ld-nscale%4.2f/cg-vacuum.crd\n' % self._nScale
        String += '\n'
        String += 'bomlev -2\n'
        String += '\n'
        String += '! open trajectory for reading\n'
        String += '    open unit 10 read unform name %s/ld-nscale%4.2f/run-ld-%d.dcd\n' % (self._basePath,self._nScale,self._temp)
#       String += '    open unit 10 read unform name /v/bigbox12/home/tim/projects/cg_model/1prb/ld-nscale%4.2f/run-ld-%d.dcd\n' % (self._nScale,self._temp)
#       String += '    open unit 10 read unform name /v/bigbox12/home/tim/projects/cg_model/2qmt/ld-nscale%4.2f/run-ld-%d.dcd\n' % (self._nScale,self._temp)
        String += '\n'
        return String
