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
from structure.models import Task, WorkingFile
import os, copy

class minimizeTask(Task):
    sdsteps = models.PositiveIntegerField(default=0)
    abnrsteps = models.PositiveIntegerField(default=0)
    tolg = models.FloatField(null=True)
    usepbc = models.CharField(max_length=1,null=True)
    useqmmm = models.CharField(max_length=1,null=True)
    qmmmsel = models.CharField(max_length=250,null=True)

    def finish(self):
        """test if the job suceeded, create entries for output"""

        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier

        # There's always an input file, so create a WorkingFile
        # for it.
        wfinp = WorkingFile()
        wfinp.task = self
        wfinp.path = loc + '/' + bnm + '-minimize' + '.inp'
        wfinp.canonPath = wfinp.path
        wfinp.type = 'inp'
        wfinp.description = 'minimization script output'
        wfinp.save()


        # Check if an output file was created and if so create
        # a WorkingFile for it.
        try:
            os.stat(loc + '/' + bnm + '-minimize' + '.out')
        except:
            self.status = 'F'
            return

        wfout = WorkingFile()
        wfout.task = self
        wfout.path = loc + '/' + bnm + '-minimize' + '.out'
        wfout.canonPath = wfout.path
        wfout.type = 'out'
        wfout.description = 'minimization script output'
        wfout.save()

        if self.status == 'F':
            return

        # check and make sure that the output PSF/CRD were 
        # created
        try:
            os.stat(loc + '/' + bnm + '-minimization.crd')
        except:
            self.status = 'F'
            self.save()
            return


        # create Working files for PDB, CRD, and PSF.
        wf = WorkingFile()
        wf.task = self
        wf.path = loc + '/' + bnm + '-minimization.crd'
        wf.canonPath = wf.path
        wf.type = 'crd'
        wf.description = 'minimized structure'
        wf.pdbkey = 'mini_' + self.workstruct.identifier
        wf.save()
        self.workstruct.addCRDToPickle(wf.path,'mini_' + self.workstruct.identifier)

        wfpsf = WorkingFile()
        wfpsf.task = self
        wfpsf.path = loc + '/' + bnm + '-minimization.psf'
        wfpsf.canonPath = wfpsf.path
        wfpsf.type = 'psf'
        wfpsf.description = 'minimized structure'
        wfpsf.save()

        wfpdb = WorkingFile()
        wfpdb.task = self
        wfpdb.path = loc + '/' + bnm + '-minimization.pdb'
        wfpdb.canonPath = wfpdb.path
        wfpdb.type = 'pdb'
        wfpdb.description = 'minimized structure'
        wfpdb.save()

        self.status = 'C'
