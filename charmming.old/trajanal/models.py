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
from django import forms

# trajectory analysis parameters
class trajAnalParams(models.Model):
    selected      = models.CharField(max_length=1)
    trajFileName  = models.CharField(max_length=250)
    trajStartStep = models.PositiveIntegerField(default=1000)
    trajStopStep  = models.PositiveIntegerField(default=10000)
    trajSkip      = models.PositiveIntegerField(default=100)
    atomSelect    = models.CharField(max_length=50)

class trajanalFileForm(forms.Form):
    startstep = forms.CharField(max_length=5)
    stopstep = forms.CharField(max_length=5)
    skip = forms.CharField(max_length=5)
    atomselection = forms.CharField(max_length=50)
    trjfile = forms.Field( widget = forms.FileInput() )
