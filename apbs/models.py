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
import django.forms
import structure

class apbsParams(models.Model):
    gdim_x = models.PositiveIntegerField(default=0)
    gdim_y = models.PositiveIntegerField(default=0)
    gdim_z = models.PositiveIntegerField(default=0)
    npoint_x = models.PositiveIntegerField(default=0)
    npoint_y = models.PositiveIntegerField(default=0)
    npoint_z = models.PositiveIntegerField(default=0)

class redoxParams(models.Model):
    selected = models.CharField(max_length=1)
    pdb = models.ForeignKey(structure.models.Structure,null=True)
    statusHTML = models.CharField(max_length=250)
    segments = models.CharField(max_length=250)
    redoxsite = models.CharField(max_length=250) # segid + 1-/2- or 2-/3-
    apbsparams = models.ForeignKey(apbsParams,null=True)
