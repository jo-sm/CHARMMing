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
from django.template import *
from django.http import HttpResponseRedirect, HttpResponse
from django.template.loader import get_template
from django.db import models
from django import forms
from django.contrib.auth.models import User
from structure.models import WorkingStructure
from structure.models import Task

class AtomSelection(models.Model):

    workstruct = models.ForeignKey(WorkingStructure)
    task = models.ForeignKey(Task)
    selectionstring = models.CharField(max_length=4000,default=None)
    linkatom_num = models.PositiveIntegerField(default=0)

class LonePair(models.Model):
    divid = models.PositiveIntegerField(default=0) #This is the same for both - it's linkqm + (divid) or linkmm + (divid), since they're pairs, and same for qmatomtype/mmatomtype
    qmresid = models.PositiveIntegerField(default=0) #These two can be very different.
    mmresid = models.PositiveIntegerField(default=0)
    qmatomname = models.CharField(max_length=10,default=None)
    mmatomname = models.CharField(max_length=10,default=None)
    qmsegid = models.CharField(max_length=40,default=None)
    mmsegid = models.CharField(max_length=40,default=None) #These will usually be the same, but I'm going for the super-generic case
    selection = models.ForeignKey(AtomSelection)

