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

    #The NULLs are all present because MM selections don't have any of those fields. However, they much be attached to a task and a workingstructure and have a selection
    workstruct = models.ForeignKey(WorkingStructure)
    task = models.ForeignKey(Task) #HOlds the parent task for this atom selection. NOT the one that was performed on it, but the one where the original coordinates come from.
    #Not sure why we have that field.
    selectionstring = models.CharField(max_length=4000,default="all") #Default "all" to make outer-layer selections easier; i.e. sele all 
    linkatom_num = models.PositiveIntegerField(default=None,null=True)
    selection_type = models.CharField(max_length=100,default="qmmm") #"qmmm" for standard, "oniom" for MSCALE, open to additional fields like semi_emp
    #Following are the other parameters so the user doesn't have to fill them out
    exchange = models.CharField(max_length=30,default=None,null=True)
    charge = models.IntegerField(default=None,null=True)
    correlation = models.CharField(max_length=10,default=None,null=True)
    basis_set = models.CharField(max_length=30,default=None,null=True)
    multiplicity = models.PositiveIntegerField(default=None,null=True)
    time_created = models.DateTimeField(auto_now=True)

class OniomSelection(AtomSelection): #Holds the MSCALE bits
    total_layers = models.PositiveIntegerField(default=2) #Holds how many layers this ONIOM model has in total
    layer_num = models.PositiveIntegerField(default=1) #Holds which layer (1,2,3,etc.) this selection is
    isQM = models.BooleanField(default=False) #Holds whether this is a QM layer.

class LonePair(models.Model):
    divid = models.CharField(max_length=50,default="0") #This is the same for both - it's linkqm + (divid) or linkmm + (divid), since they're pairs, and same for qmatomtype/mmatomtype
    qmresid = models.PositiveIntegerField(default=0) #These two can be very different.
    mmresid = models.PositiveIntegerField(default=0)
    qmatomname = models.CharField(max_length=10,default=None)
    mmatomname = models.CharField(max_length=10,default=None)
    qmsegid = models.CharField(max_length=40,default=None)
    mmsegid = models.CharField(max_length=40,default=None) #These will usually be the same, but I'm going for the super-generic case
    selection = models.ForeignKey(AtomSelection)
    qqh = models.PositiveIntegerField(default=1) #Simplifies accounting for these things in MSCALE scripts

