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
from structure.models import WorkingFile
import structure
import django.forms

class nmodeParams(models.Model):
    structure = models.ForeignKey(structure.models.Structure,null=True)
    inpStruct = models.ForeignKey(WorkingFile,null=True)

    statusHTML = models.CharField(max_length=250)
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

# Create your models here.
