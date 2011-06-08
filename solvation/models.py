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
import structure

class solvationParams(models.Model):
    structure = models.ForeignKey(structure.models.WorkingStructure)

    solvation_structure = models.CharField(max_length=50)
    solv_pref = models.CharField(max_length=50)

    #If the user has no preference on the structure size ("Let the GUI take care of it for me" option)
    #Then no_pref_radius stores the distance from structure to edge of solvation structure
    no_pref_radius = models.DecimalField(max_digits=8,decimal_places=4,default=0)

    #If the user chooses to manually set preferences ("I want to set my own dimensions" option)
    #pref_x,y, and z will store that information
    xtl_x = models.DecimalField(max_digits=8,decimal_places=4,default=0)
    xtl_y = models.DecimalField(max_digits=8,decimal_places=4,default=0)
    xtl_z = models.DecimalField(max_digits=8,decimal_places=4,default=0)

    #The below are used for neutralization
    salt = models.CharField(max_length=5,null=True)
    concentration = models.DecimalField(max_digits=8,decimal_places=6,default=0)
    ntrials = models.PositiveIntegerField(default=0)
