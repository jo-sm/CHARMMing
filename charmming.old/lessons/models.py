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

class Lesson(models.Model):
    # data for lessons (should not be overridden by subclasses)
    # specifying both user and PDBFile is redundant (since the PDBFile references the user),
    # but makes it easier to look up all lessons being done by a particular user.
    user = models.ForeignKey(User)
    nSteps = models.PositiveIntegerField(default=0)
    curStep = models.PositiveIntegerField(default=0)

    # lesson methods (subclasses should override as needed)
    def onFileUpload(self,postdata):
        return True

    def onEditPDBInfo(self,postdata):
        return True

    def onMinimizeSubmit(self,postdata):
        return True

    def onMinimizeDone(self,file):
        return True

    def onSolvationSubmit(self,postdata):
        return True

    def onSolvationDone(self,file):
        return True

    def onEnergySubmit(self,postdata):
        return True

    def onEnergyDone(self,file):
        return True

    def onNMASubmit(self,postdata):
        return True

    def onNMADone(self,file):
        return True

    def onMDSubmit(self,postdata):
        return True

    def onMDDone(self,file):
        return True

    def onLDSubmit(self,postdata):
        return True

    def onLDDone(self,file):
        return True
    
    def onSGLDSubmit(self,postdata):
        return True

    def onSGLDDone(self,file):
        return True

    def onBuildStructureSubmit(self,postdata):
        return True

    def onBuildStructureDone(self,file):
        return True

    def onRedoxSubmit(self,postdata):
        return True

    def onRedoxDone(self,file):
        return True

    def onRMSDSubmit(self,postdata):
        return True

"""
Eventually this should be included as its own model but temporarily 
this will be left here:
class Lesson1
def onFileUpload(self,postdata):

"""
class LessonProblem(models.Model):
    class Admin:
        pass
    # data model
    lesson_type = models.CharField(max_length=50)
    lesson_id = models.PositiveIntegerField(default=0)
    errorstep = models.PositiveIntegerField(default=0)
    severity = models.PositiveIntegerField(default=0)
    description = models.CharField(max_length=250)

    # only method right now is the constructor...
  #  def __init__(self,type,id,severity,description):
  #      self.lesson_type = type
  #      self.lesson_id = id
  #      self.severity = severity
  #      self.description = description
  #      self.save()

