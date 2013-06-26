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
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
import lesson1, charmming_config
import mimetypes
import re, os, copy

def downloadLessonFiles(request,lessonfolder, filename,mimetype = None):
    username = request.user.username

    fnarr = filename.split("/")
    filename = fnarr[-1]

    path = "%s/mytemplates/lessons/%s" % (charmming_config.charmming_root, lessonfolder)
    if mimetype is None:
        mimetype,encoding = mimetypes.guess_type(path + "/" + filename)
    response = HttpResponse(mimetype=mimetype)
    response['Content-Disposition'] = 'attachment; filename=%s' %filename
    response.write(file(path + "/" + filename, "rb").read())
    return response
   
def doLessonAct(file,postdata,function,filename=None):
     lessontype = file.lesson_type
     lessonid = file.lesson_id
     if lessontype == 'lesson1':
         lesson_obj = lesson1.models.Lesson1.objects.filter(id=lessonid)[0]

     if function == 'onFileUpload':
         lesson_obj.onFileUpload(postdata)
         return True
     elif function == 'onMinimizeSubmit':
         lesson_obj.onMinimizeSubmit(postdata,filename)
         return True
     elif function == 'onSolvationSubmit':
         lesson_obj.onSolvationSubmit(postdata,filename)
         return True
     elif function == 'onMDSubmit':
         lesson_obj.onMDSubmit(postdata,filename)
         return True
     elif function == 'onRedoxSubmit':
         lesson_obj.onRedoxSubmit(postdata,filename)
         return True
     elif function == 'onRedoxDone':
         lesson_obj.onRedoxDone(postdata,filename)
         return True

