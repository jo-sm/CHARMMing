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
from structure.models import Task, WorkingFile

class apbsParams(models.Model):
    gdim_x = models.PositiveIntegerField(default=0)
    gdim_y = models.PositiveIntegerField(default=0)
    gdim_z = models.PositiveIntegerField(default=0)
    npoint_x = models.PositiveIntegerField(default=0)
    npoint_y = models.PositiveIntegerField(default=0)
    npoint_z = models.PositiveIntegerField(default=0)

class redoxTask(Task):
    redoxsite = models.CharField(max_length=250) # segid + 1-/2- or 2-/3-
    apbsparams = models.ForeignKey(apbsParams,null=True)

    def finish(self):
        # This is a pain, since we need to get all of the input, output, and
        # PDB files from the REDOX calculation.

        self.status = 'C'
