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
from django import forms
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from account.views import isUserTrustworthy
from structure.models import Structure
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from minimization.models import minimizeTask
import re
import copy
import os
import lessonaux
import lessons

   
def displayLessonStatus(request):
    try:
        file = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse('')

    try:
        lesson_obj = lessonaux.doLessonAct(file,"getObject")
    except:
        return HttpResponse('')
    try:
        lessonprob = lessons.models.LessonProblem.objects.get(lesson_type = file.lesson_type,lesson_id = file.lesson_id)
    except:
        lessonprob = None
    step_status_list = lesson_obj.generateStatusHtml(file)

    #Gets a list of numbers that corresponds to whether a lesson step was successfully done, failed,
    #or has not been performed yet
    return render_to_response('html/lesson1report.html', {'file':file,'lesson':lesson_obj, 'lessonprob' : lessonprob,'status_list':step_status_list})

