
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
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from account.models import *
from structure.models import Task
from django.contrib.auth.models import User
from django.template import *
from account.views import checkPermissions
import structure.models
from lesson_config import *
import pychm.io, charmming_config
from selection.models import AtomSelection

def selectstructure(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    if request.POST:
        postdata=request.POST
        source = postdata['source']
        dest = source + "/" #Yes this is a bad naming convention. However, it works.
        try:
            struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
        except:
            return HttpResponse("No structure selected.")

        try:
            ws = structure.models.WorkingStructure.objects.get(structure=struct,selected='y')
        except:
            return HttpResponse("Please build a working structure.")

        task_id = postdata['task_id']

        try:
            task = Task.objects.filter(workstruct=ws,id=task_id)[0] #This might break with build tasks?
        except:
            return HttpResponse("Task " + task_id + " does not exist.")
        if postdata.has_key('atomselection'):
            try: #First see if there's already an atomselection present
                oldselect = AtomSelection.objects.filter(workstruct=ws)[0]
                oldselect.task = task
                oldselect.workstruct = ws #Just in case...
                oldselect.selectionstring = postdata['atomselection'] #replace it
                oldselect.save()
            except: #There's two exceptions possible, so first we try to make a new atomselect
                atomselect = AtomSelection()
                atomselect.task = task
                atomselect.workstruct = ws
                try:
                    atomselect.selectionstring = postdata['atomselection']
                    atomselect.save()
                except: #THis is if you have so many atoms you break the postdata request size
                    return HttpResponse("You have selected too many atoms.")
            return HttpResponseRedirect('/charmming/'+ dest)
            #Do stuff with the database and return redirect to source
        tdict = {}
        filepath = struct.location.replace(charmming_config.user_home,'') + "/" + ws.identifier + "-" + task.action + ".pdb"
        #for example 1yjp-ws1-build.pdb
        tdict['filepath'] = filepath
        tdict['source'] = source
        tdict['task_id'] = task_id
        return render_to_response("html/selection.html", tdict)
    else:
        return HttpResponse("No source.")
