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


#Replica exchange parameters
class rexParams(models.Model):

    firstStep = models.PositiveIntegerField(default=1)
    lastStep = models.PositiveIntegerField(default=2)
    cleanStep = models.PositiveIntegerField(default=50)
    npref = models.DecimalField(default=0,null=False,max_digits=8,decimal_places=3)
    nbath = models.PositiveIntegerField(default=4)
    temperatures = models.CharField(max_length=250)

class mdParams(models.Model):

    pdb = models.ForeignKey(structure.models.Structure,null=True)
    statusHTML = models.CharField(max_length=250)
    selected = models.CharField(max_length=1)
    sequence = models.PositiveIntegerField(default=1)
    type = models.CharField(max_length=50)
    nstep = models.PositiveIntegerField(default=1000)
    firstt = models.DecimalField(default=310.15,max_digits=8,decimal_places=3,null=True)

    #temp will represent the temperature Kelvin if "type" is heat
    temp = models.DecimalField(default=410.15,null=True,max_digits=8,decimal_places=3)
    finalt = models.DecimalField(default=410.15,null=True,max_digits=8,decimal_places=3)
    teminc = models.DecimalField(default=10.0,null=True,max_digits=8,decimal_places=3)
    ihtfrq = models.DecimalField(default=100.0,null=True,max_digits=8,decimal_places=3)
    tbath = models.DecimalField(default=410.15,null=True,max_digits=8,decimal_places=3)
    scipism = models.BooleanField(default=False)
    make_md_movie = models.BooleanField(default=False)
    md_movie_status = models.CharField(max_length=250,null=True)
    md_movie_req = models.BooleanField(default=False)

    replica_exchange = models.ForeignKey(rexParams,null=True)
# Create your models here.

class ldParams(models.Model):

    pdb = models.ForeignKey(structure.models.Structure,null=True)
    statusHTML = models.CharField(max_length=250)
    selected = models.CharField(max_length=1)
    nstep = models.PositiveIntegerField(default=1000)
    fbeta = models.DecimalField(default=60.0,null=True,max_digits=8,decimal_places=3)
    scipism = models.BooleanField(default=False)
    sgld = models.BooleanField(default=False)
    make_ld_movie = models.BooleanField(default=False)
    ld_movie_status = models.CharField(max_length=250,null=True)
    ld_movie_req = models.BooleanField(default=False)
    replica_exchange = models.ForeignKey(rexParams,null=True)
# Create your models here.

class sgldParams(models.Model):

    pdb = models.ForeignKey(structure.models.Structure,null=True)
    statusHTML = models.CharField(max_length=250)
    selected = models.CharField(max_length=1)
    nstep = models.PositiveIntegerField(default=1000)
    fbeta = models.DecimalField(default=60.0,null=True,max_digits=8,decimal_places=3)
    tsgavg = models.DecimalField(default=0.5,null=True,max_digits=8,decimal_places=6)
    tempsg = models.DecimalField(default=1.0,null=True,max_digits=8,decimal_places=6)
    scipism = models.BooleanField(default=False)
    sgld = models.BooleanField(default=False)
    make_sgld_movie = models.BooleanField(default=False)
    sgld_movie_status = models.CharField(max_length=250,null=True)
    sgld_movie_req = models.CharField(max_length=250,null=True)
    sgld_movie_req = models.BooleanField(default=False)
    replica_exchange = models.ForeignKey(rexParams,null=True)

# Create your models here.
