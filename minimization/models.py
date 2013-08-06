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
from structure.models import Task, WorkingFile, CGWorkingStructure
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
        basepath = loc + '/' + bnm + "-" + self.action
        # There's always an input file, so create a WorkingFile
        # for it.

        path = basepath + ".inp"
        wfinp = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfinp.task = self
            wfinp.path = path
            wfinp.canonPath = wfinp.path
            wfinp.type = 'inp'
            wfinp.description = 'minimization script input'
            wfinp.save()


        # Check if an output file was created and if so create
        # a WorkingFile for it.
        try:
            os.stat(basepath + '.out')
        except:
            self.status = 'F'
            return

        path = basepath + ".out"
        wfout = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfout.task = self
            wfout.path = path
            wfout.canonPath = wfout.path
            wfout.type = 'out'
            wfout.description = 'minimization script output'
            wfout.save()

        if self.status == 'F':
            return

        # check and make sure that the output PSF/CRD were 
        # created
        try:
            os.stat(basepath + '.crd')
        except:
            self.status = 'F'
            self.save()
            return


        path = basepath + ".crd"
        # create Working files for PDB, CRD, and PSF.
        wf = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wf.task = self
            wf.path = loc + '/' + bnm + '-minimization.crd'
            wf.canonPath = wf.path
            wf.type = 'crd'
            wf.description = 'minimized structure'
            wf.pdbkey = 'mini_' + self.workstruct.identifier
            wf.save()
            self.workstruct.addCRDToPickle(wf.path,'mini_' + self.workstruct.identifier)

        path = basepath + ".psf"
        wfpsf = WorkingFile()
        inp_file = path
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfpsf.task = self
            wfpsf.path = loc + '/' + bnm + '-minimization.psf'
            wfpsf.canonPath = wfpsf.path
            wfpsf.type = 'psf'
            wfpsf.description = 'minimized structure'
            wfpsf.save()

        path = basepath + ".pdb"
        out_file = path
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfpdb = WorkingFile()
            wfpdb.task = self
            wfpdb.path = loc + '/' + bnm + '-minimization.pdb'
            wfpdb.canonPath = wfpdb.path
            wfpdb.type = 'pdb'
            wfpdb.description = 'minimized structure'
            wfpdb.save()

        #Generic coarse-grain code goes here.
        try:
            cgws = CGWorkingStructure.objects.get(workingstructure_ptr=self.workstruct)
        except CGWorkingStructure.MultipleObjectsReturned: #Uh oh. This MAY be alright if AA/CG...
            self.status = "F"
            return
        except: #Catch everything else
            pass

        if cgws:
            cgws.addBondsToPDB(inp_file,out_file) #Sometimes the DB is too slow to catch up and returns blank.
        self.status = 'C'
