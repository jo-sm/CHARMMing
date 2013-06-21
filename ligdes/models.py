
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
# Create your views here.
from django import forms
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from account.models import *
from structure.models import Task 
from django.contrib.auth.models import User
from django.template import *
from account.views import checkPermissions
from lessons.models import LessonProblem
from mutation.models import mutateTask
import output, lesson1, lesson2, lesson3, lesson4, lessonaux
import structure.models, input
from lesson_config import *
import os, copy, json, mimetypes, string, re
import pychm.io, charmming_config
from pychm.const.bio import aaAlphabet
import cPickle
