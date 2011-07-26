#
#                            PUBLIC DOMAIN NOTICE
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the authors' official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software is freely available
#  to the public for use.  There is no restriction on its use or
#  reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, NIH and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. NIH, NHLBI, and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any
#  particular purpose.
from django.db import models
import structure.Task, structure.WorkingFile
import os, copy

class minimizeTask(structure.Task):
    selected = models.CharField(max_length=1)

    sdsteps = models.PositiveIntegerField(default=0)
    abnrsteps = models.PositiveIntegerField(default=0)
    tolg = models.DecimalField(max_digits=8,decimal_places=5)
    usepbc = models.CharField(max_length=1)
    useqmmm = models.CharField(max_length=1)
    qmmmsel = models.CharField(max_length=250)

    def finish(self,sstring):
        """test if the job suceeded, create entries for output"""

        loc = self.workstruct.structure.location
        bnm = 'mini-' + self.workstruct.identifier

        # Check if an output file was created and if so create
        # a WorkingFile for it.
        try:
            os.stat(loc + '/mini-' + bnm + '.out')
        except:
            self.status = 'F'
            return

        wfout = structure.WorkingFile()
        wfout.task = self
        wfout.path = loc + '/mini-' + bnm + '.out'
        wfout.canonPath = wf.path
        wfout.type = 'out'
        wfout.description = 'minimization script output'
        wfout.save()

        if sstring == 'failed':
            self.status = 'F'
            return

        # check and make sure that the output PSF/CRD were 
        # created
        try:
            os.stat(loc + '/mini-' + bnm + '.crd')
        except:
            self.status = 'F'
            return


        # create Working files for PDB, CRD, and PSF.
        wf = structure.WorkingFile()
        wf.task = self
        wf.path = loc + '/mini-' + bnm + '.crd'
        wf.canonPath = wf.path
        wf.type = 'crd'
        wf.description = 'minimized structure'
        wf.pdbkey = 'mini_' + self.workstruct.identifier
        wf.save()
        self.workstruct.addCRDToPickle(wf.path,'mini_' + self.workstruct.identifier)

        wfpsf = copy.deepcopy(wf)
        wfpsf.path = loc + '/mini-' + bnm + '.psf'
        wfpsf.canonPath = wf.path
        wfpsf.type = 'psf'
        wfpsf.save()

        wfpdb = copy.deepcopy(wf)
        wfpdb.path = loc + '/mini-' + bnm + '.pdb'
        wfpdb.canonPath = wf.path
        wfpdb.type = 'pdb'
        wfpdb.save()

        self.status = 'C'
