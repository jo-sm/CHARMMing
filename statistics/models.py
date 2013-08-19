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
import structure.models #For all its database goodness

#Statistics are made up of data points...thus the name
class DataPoint(models.Model): #I don't have a good name for this yet.
    #NONE of these fields must be ForeignKeys!
    #All this data MUST persist in the DB permanently. Otherwise
    #it becomes worthless.
    task_id = models.PositiveIntegerField(default=0)
    task_action = models.CharField(max_length=100)
    user = models.CharField(max_length=30)
    structure_name = models.CharField(max_length=100)
    struct_id = models.CharField(max_length=100) #This is here because otherwise it is very difficult to avoid multiple-counting structures
    success = models.BooleanField(default=False)
    date_created = models.DateTimeField(auto_now=True)
