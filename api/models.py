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
import pychm
from django.db import models
from minimization.models import MinimizeTask
from structure.models import EnergyTask,WorkingFile

class APIUser(models.Model):
    callerName = models.CharField(max_length=40,null=False,default=None)
    callerKey = models.CharField(max_length=30,null=False,default=None)

class APIJob(models.Model):
    jobType = models.PositiveIntegerField(default=0)
    user = models.ForeignKey(APIUser)
    timestamp = models.DateField(auto_now_add=True)
    directory = models.CharField(max_length=30,null=False)

class APIEnergy(APIJob):
    implicitSolvent = models.CharField(max_length=8,null=False,default=None)
    EnergyValue = models.FloatField(null=True)
    task = models.ForeignKey(EnergyTask,null=True)

    def run(self,req):
        response = {}
        
        response['errcode'] = -999
        return response

class APIOptimization(APIJob):
    implicitSolvent = models.CharField(max_length=8,null=False,default=None)
    nSDStep = models.PositiveIntegerField(default=0)
    nABNRStep = models.PositiveIntegerField(default=0)
    FinalEnergyValue = models.FloatField(null=True)
    OptimizedPDB = models.Foreignkey(WorkingFile)
    task = models.ForeignKey(MinimizeTask,null=True)

    def run(self,req):
        response = {}

        return response
