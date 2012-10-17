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
# Create your views here.
from httplib import HTTPConnection
from django import forms
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from pdbinfo.models import PDBFile, PDBFileForm, ParseException, energyParams
from minimization.models import minimizeParams
from minimization.views import append_tpl
from dynamics.models import mdParams, ldParams, sgldParams
from solvation.models import solvationParams
from account.models import *
from dynamics.views import combinePDBsForMovie
from normalmodes.views import combineNmaPDBsForMovie
from normalmodes.aux import getNormalModeMovieNum
from normalmodes.models import nmodeParams
from apbs.models import redoxParams
from pdbinfo.qmmm import makeQChem_tpl, handleLinkAtoms, writeQMheader
from django.contrib.auth.models import User
from django.core import validators 
from django.core.mail import mail_admins
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from account.views import isUserTrustworthy
from pdbinfo.editscripts import generateHTMLScriptEdit
from pdbinfo.aux import checkNterPatch
#dd stuff
from dd_target.models import proteins, protein_conformations 
from dd_infrastructure.models import files,files_objects,file_types
###
import output
import lesson1
import lesson2
import lesson3
import lesson4
import lessonaux
# import all created lessons by importing lesson_config
# also there is a dictionary called 'file_type' in lesson_config.py specififying the file type of files uploaded by the lessons  
from lesson_config import *
import os
import re
import copy
import datetime
import time
import mimetypes
import time
import stat
import string
import random
import glob
import sys, traceback
import commands
import charmming_config
import common


def newUpload(request, template="html/ddtargetfileupload.html"):

    command_list=[]
    uploadtargetlog=open("/tmp/targetupload.log",'w')
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    newtarget = proteins()
    newfile = files()
    newfileobject= files_objects()

    try:
        request.FILES['target_file']['filename']
        target_uploaded_by_user = 1
    except:
        target_uploaded_by_user = 0

    if target_uploaded_by_user==1:
       
	try: 
	    newtarget=proteins.objects.filter(owner=u,id=request.POST['existingtarget'])[0]
	except:
            newtarget.owner=u
	    newtarget.protein_name=request.POST['targetname']
            newtarget.pdb_code=request.POST['pdbcode']
            newtarget.protein_owner_index=str(common.getNewProteinOwnerIndex(u))
	    newtarget.description=request.POST['targetdescription']
           
            newtarget.save()
            newtarget = target.objects.filter(owner=u).order_by("-id")[0]

        try: 
	    newconformation=protein_conformations.objects.filter(owner=u,id=request.POST['existingconformation'])[0]
	except:
            newconformation.owner=u
            newconformation.protein=newtarget
	    newconformation.conformation_name=request.POST['conformationname']
            newconformation.conformation_protein_index=str(common.getNewConformationProteinIndex(newtarget))
	    newconformation.description=request.POST['targetdescription']
            
            newconformation.save()
            newconformation = target_conformations.objects.filter(owner=u).order_by("-id")[0]

        location = charmming_config.user_dd_targets_home + '/' + username + '/' + 'target_' + str(newtarget.protein_owner_index)
    	filename='target_' + str(newtarget.protein_owner_index) + '_' + str(newconformation.conformation_protei_index)+'.mol2'
    	destination = open(location + filename,'w')
        destination.write(request.FILES['target_file']['content'])
	destination.close()
	target_file_type=file_types.objects.filter(name="Target Structure File")[0]
	newfile.owner=u
	newfile.filename=filename
	newfile.location=location
	newfile.file_type=target_file_typ
        newfile.description=request.POST['filedesc']
	newfile.save()
	newfile = Files.objects.filter(owner=u).order_by("-id")[0]

	newfileobject.file=newfile
	newfileobject.object_table_name='dd_target_protein_conformations'
	newfileobject.object_id=newconformation.id
	newfileobject.save()
    
    return render_to_response('html/ddtargetfileupload.html')

def targetFileUploadForm(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    
    form = TargetFileForm()
    return render_to_response('html/ddtargetfileuploadform.html', {'form':form})

def targetDropdown(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    
    targetlist=proteins.objects.filter(owner=request.user).order_by("protein_name")
    return render_to_response('html/ddtargetdropdown.html', {'targetlist':targetlist})

class TargetFileForm(forms.Form):
    target_file = forms.FileField()

	
