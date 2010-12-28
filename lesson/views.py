from django import forms
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from account.views import isUserTrustworthy
from pdbinfo.models import PDBFile, PDBFileForm
from lessons.models import LessonProblem
from lesson.models import Lesson
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
import re
import copy
import os

#Displays lesson  page
def lessonDisplay(request):
    PDBFile().checkRequestData(request)
    try:
        file = PDBFile.objects.filter(owner=request.user,selected='y')[0]
    except:
        return render_to_response('html/lesson.html')
    #If its a lesson object, get the id by the file id
    if file.lesson_type == 'lesson':
        lesson_obj = Lesson.objects.filter(user=request.user,id=file.lesson_id)[0]
        html_step_list = lesson_obj.getHtmlStepList()
    else:
        lesson_obj = None
        html_step_list = None
    try:
        lessonprob_obj = LessonProblem.objects.filter(lesson_type='lesson',lesson_id=lesson_obj.id)[0]
    except:
        lessonprob_obj = None
    return render_to_response('html/lesson.html',{'lesson':lesson_obj,'lessonproblem':lessonprob_obj,'html_step_list':html_step_list})
