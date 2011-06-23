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
import structure.models
import math

class solvationParams(models.Model):
    structure = models.ForeignKey(structure.models.WorkingStructure)
    inpStruct = models.ForeignKey(WorkingFile,null=True)

    # these contain the shape of the unit cell and its dimensions
    solvation_structure = models.CharField(max_length=50)
    xtl_x = models.DecimalField(max_digits=8,decimal_places=4,default=0)
    xtl_y = models.DecimalField(max_digits=8,decimal_places=4,default=0)
    xtl_z = models.DecimalField(max_digits=8,decimal_places=4,default=0)
    spradius = models.DecimalField(max_digits=8,decimal_places=4,default=0)
    selected = models.CharField(max_length=1)
    statusHTML = models.CharField(max_length=250)

    #The below are used for neutralization
    salt = models.CharField(max_length=5,null=True)
    concentration = models.DecimalField(max_digits=8,decimal_places=6,default=0)
    ntrials = models.PositiveIntegerField(default=0)

    def calcEwaldDim(self):
        """
        Find the nearest integer that is greater than two times the longest
        crystal box length that is a multiple of ONLY 2, 3, and 5.
        """

        def facts(n):
            r = []

            mn = int(math.sqrt(n))+1
            for i in range(2,mn):
                if n % i == 0: r.append(i)
            return r

        trgt = int(2 * max(self.xtl_x,self.xtl_y,self.xtl_z))

        while True:
           trgt += 1
           good = True
           for f in facts(trgt):
               if not f in [2,3,5]:
                   good = False
                   break
           if good:
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
