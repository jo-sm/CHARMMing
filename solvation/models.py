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
from django.contrib.auth.models import User
from django.contrib import admin
from scheduler.schedInterface import schedInterface
from structure.models import WorkingFile, Task
import math, os

class solvationTask(Task):
    # these contain the shape of the unit cell and its dimensions
    solvation_structure = models.CharField(max_length=50,null=True) # null=True to allow early save()ing of the model
    xtl_x = models.FloatField(default=0.0)
    xtl_y = models.FloatField(default=0.0)
    xtl_z = models.FloatField(default=0.0)
    spradius = models.FloatField(default=0.0)

    # The below are used for neutralization
    salt = models.CharField(max_length=5,null=True)
    concentration = models.FloatField(default=0.0)
    ntrials = models.PositiveIntegerField(default=0)

    def finish(self):
        """test if the job suceeded, create entries for output"""

        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier
        basepath = loc + "/" + bnm + "-" + self.action

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
            wfinp.description = 'solvation input script'
            wfinp.save()

        if self.concentration > 0.0001:
            path = loc + "/" + bnm + "-neutralize.inp"
            wfninp = WorkingFile()
            try:
                wftest = WorkingFile.objects.get(task=self,path=path)
            except:
                wfninp.task = self
                wfninp.path = path
                wfninp.canonPath = wfninp.path
                wfninp.type = 'inp'
                wfninp.description = 'neutralization input script'
                wfninp.save()

        # Check if an output file was created and if so create
        # a WorkingFile for it.
        path = basepath + ".out"

        try:
            os.stat(path)
        except:
            self.status = 'F'
            return

        wfout = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfout.task = self
            wfout.path = path
            wfout.canonPath = wfout.path
            wfout.type = 'out'
            wfout.description = 'solvation script output'
            wfout.save()

        if self.concentration > 0.0001:
            path = loc + "/" + bnm + "-neutralize.out"
            wfnout = WorkingFile()
            try:
                wftest= WorkingFile.objects.get(task=self,path=path)
            except:
                wfnout.task = self
                wfnout.path = path
                wfnout.canonPath = wfnout.path
                wfnout.type = 'out'
                wfnout.description = 'neutralization script output'
                wfnout.save()

        # now check if all the expected psf/crd/pdb files exist
        path = basepath + '.crd'
        try:
            os.stat(path)
        except:
            self.status = 'F'
            self.save()
            return

        wf = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wf.task = self
            wf.path = path
            wf.canonPath = wf.path
            wf.type = 'crd'
            wf.description = 'solvated structure'
            wf.pdbkey = 'solv_' + self.workstruct.identifier
            wf.save()
            self.workstruct.addCRDToPickle(wf.path, 'neut_' + self.workstruct.identifier)

        path = basepath + ".psf"
        wfpsf = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfpsf.task = self
            wfpsf.path = path
            wfpsf.canonPath = wfpsf.path
            wfpsf.type = 'psf'
            wfpsf.description = 'solvated structure'
            wfpsf.save()

        path = basepath + ".pdb"
        wfpdb = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfpdb.task = self
            wfpdb.path = path
            wfpdb.canonPath = wfpdb.path
            wfpdb.type = 'pdb'
            wfpdb.description = 'solvated structure'
            wfpdb.save()

        path = loc + "/" + bnm + "-neutralization.crd"
        if self.concentration > 0.0001:
            try:
                os.stat(path)
            except:
                self.status = 'F'
                self.save()
                return


            wfn = WorkingFile()
            try:
                wftest = WorkingFile.objects.get(task=self,path=path)
            except:
                wfn.task = self
                wfn.path = path
                wfn.canonPath = wfn.path
                wfn.type = 'crd'
                wfn.description = 'neutralized structure'
                wfn.pdbkey = 'neut_' + self.workstruct.identifier
                wfn.save()
                self.workstruct.addCRDToPickle(wf.path,'neut_' + self.workstruct.identifier)

            path = loc + "/" + bnm + "-neutralization.psf"
            wfnpsf = WorkingFile()
            try:
                wftest = WorkingFile.objects.get(task=self,path=path)
            except:
                wfnpsf.task = self
                wfnpsf.path = path
                wfnpsf.canonPath = wfnpsf.path  
                wfnpsf.type = 'psf'
                wfnpsf.description = 'neutralized structure'
                wfnpsf.save()

            path = loc + "/" + bnm + "-neutralization.pdb"
            wfnpdb = WorkingFile()
            try:
                wftest = WorkingFile.objects.get(task=self,path=path)
            except:
                wfnpdb.task = self
                wfnpdb.path = path
                wfnpdb.canonPath = wfnpdb.path  
                wfnpdb.type = 'pdb'
                wfnpdb.description = 'neutralized structure'
                wfnpdb.save()


    def calcEwaldDim(self):
        """
        Find the nearest integer that is greater than two times the longest
        crystal box length that is a multiple of ONLY 2, 3, and 5.
        """

        logfp = open('/tmp/foo', 'w')

        def prime_factorize(x,li=[]):
            until = int(math.sqrt(x))+1
            for i in xrange(2,until):
                if not x%i:
                    li.append(i)
                    break
            else:                      #This else belongs to for
                    li.append(x)
                    return li
            prime_factorize(x/i,li)
            return li

        trgt = int(2 * max(self.xtl_x,self.xtl_y,self.xtl_z))

        while True:
           trgt += 1
           if trgt % 2 != 0: continue # ewald factors must be even

           good = True
           flst = prime_factorize(trgt,li=[])
           if not flst: continue
           logfp.write('trgt %d facts %s\n' % (trgt,flst))
           for f in flst:
               if not f in [2,3,5]:
                   good = False
                   break
           if good:
               logfp.write('returning %d\n' % trgt)
               logfp.close()
               return trgt
           

    def check_valid(self):
        if self.solvation_structure == 'cubic' or self.solvation_structure == 'rhdo':
            return self.xtl_x == self.xtl_y == self.xtl_z
        elif self.solvation_structure == 'tetra' or self.solvation_structure == 'hexa':
            return self.xtl_x == self.xtl_y
        else:
            return True

    @property
    def angles(self):
        if self.solvation_structure == 'cubic' or self.solvation_structure == 'tetra':
            return (90., 90., 90.)
        elif self.solvation_structure == 'rhdo':
            return (60., 90., 60.)
        elif self.solvation_structure == 'hexa':
            return (90., 90., 120.)
        else:
            return (-1., -1., -1.)

