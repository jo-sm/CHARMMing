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
from scheduler.schedInterface import schedInterface
from structure.models import WorkingFile, Task

class nmodeTask(Task):
    # type 1 = all atom, 2 - ENM
    type = models.PositiveIntegerField(default=0)

    nmodes = models.PositiveIntegerField(default=0)
    rcut = models.FloatField(default=0.0)
    kshort = models.FloatField(default=0.0)
    klong = models.FloatField(default=0.0)
    selected = models.CharField(max_length=1)
    nma_movie_status = models.CharField(max_length=250,default=None,null=True)
    make_nma_movie = models.BooleanField(default=False)
    nma_movie_req = models.BooleanField(default=False)

    def finish(self):
        """test if the job suceeded, create entries for output"""

        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier

        # There's always an input file, so create a WorkingFile
        # for it.
        wfinp = WorkingFile()
        wfinp.task = self
        wfinp.path = loc + '/' + bnm + '-nmodes.inp'
        wfinp.canonPath = wfinp.path
        wfinp.type = 'inp'
        wfinp.description = 'normal mode input script'
        wfinp.save()


        # Check if an output file was created and if so create
        # a WorkingFile for it.
        try:
            os.stat(loc + '/' + bnm + '-nmodes.out')
        except:
            self.status = 'F'
            return

        wfout = WorkingFile()
        wfout.task = self
        wfout.path = loc + '/' + bnm + '-nmodes.out'
        wfout.canonPath = wfout.path
        wfout.type = 'out'
        wfout.description = 'normal modes script output'
        wfout.save()

        if self.status == 'F':
            return

        # if there is a movie, get the trajectory
        if self.make_nma_movie:
            pass

        self.status = 'C'
