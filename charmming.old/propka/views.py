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

from django.template.loader import get_template
from django.contrib import messages as messages
from django.contrib.messages import get_messages
from django.shortcuts import render_to_response
from django.http import HttpResponseRedirect
from django.http import HttpResponse
from structure.models import Structure, WorkingStructure, CGWorkingStructure, WorkingFile, Task
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from django.template import *
import output, charmming_config, lessonaux, input, lessons, lesson1, lesson2, lesson3, lesson4, lesson5
from account.views import checkPermissions
import structure.models
import os
import traceback
import subprocess
import mimetypes #for download page

def propka_form(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    #chooses the file based on if it is selected or not
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please submit a structure.")
        return HttpResponseRedirect("/charmming/fileupload/")
#        return HttpResponse("Please submit a structure first.")
    try:
         ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "You need a working structure to get to this page. If you have reached this page in error, please inform the system administrator.")
        return HttpResponseRedirect("/charmming/buildstruct/")
#        return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before attempting an propka calculation")

    cg_model = False
    try:
        cgws = structure.models.CGWorkingStructure.objects.get(workingstructure_ptr=ws.id)
        cg_model = True
    except:
        pass
    #we probably can't do this with CG cause of the way they work...leaving them out for now ~VS
    if cg_model:
        messages.error(request, "Our PROPKA installation does not currently support coarse-grained working structures. Please build a non-coarse-grained working structure first.")
        return HttpResponseRedirect("/charmming/buildstruct/")

    tdict = {}
    tdict['ws_identifier'] = ws.identifier
    tasks = Task.objects.filter(workstruct=ws,status='C',active='y',modifies_coordinates=True)
    #we need to get the tasks to find out the files we need to put up for download
    propka_files = []
    file_loc = ws.structure.location #to make the loop a little faster
    for task in tasks:
        curr_path = "%s/%s-%s.pka" % (file_loc,ws.identifier,task.action)
        try:
            os.stat(curr_path)
            propka_files.append((task.action,("%s-%s.pka" % (ws.identifier,task.action))))
        except:
            continue
    tdict['propka_files'] = propka_files
    #we don't actually need the following logic except for the display page.
#    propka_lines = []
#    if len(propka_files) > 0:
#        for file in propka_files:
#            try:
#                propkafp = open('%s-%s-propka.pka' % (ws.identifier,task.action),'r')
#                propka_lines = propkafp.readlines()
#                propkafp.close()
#            except: #it's unlikely anything will ever happen here, so
#                logfp = open("/tmp/propka-errors.txt","w")
#                logfp.write("Something has gone horribly wrong with the PROPKA loading procedure.\n")
#                traceback.print_exc(file=logfp)
#                continue
#
#    tdict['propka_lines'] = propka_lines

    if request.POST.has_key('form_submit'):
        #We do not in fact need a Task object. WHile it would be nice in the sense of attaching
        #a file to our propka task, the problem is that we don't have input/output CHARMM stuff for this
        #PROPKA is totally CHARMM-independent, so sending a task is kind of pointless
        #by way of the task's "start()" routine
        #we'd need to do some major overriding of standard methods to get this to work
        #and it's probably better to get the page running first, and then worry about the
        #job scheduling stuff, given it's nonstandard
#        try:
#            oldtsk = propkaTask.objects.filter(workstruct=ws,active='y')[0]
#            oldtsk.active = 'n'
#            oldtsk.save()
#        except:
#            pass
#
#        et = propkaTask()
#        et.setup(ws)
#        et.active = 'y'
#        et.action = 'propka'
#        et.save()
#        if ws.isBuilt != 't':
#            try:
#                pTask = ws.build(et)
#            except noNscaleFound, e:
#                return HttpResponse('The nScale parameterization process has not yet completed. It may take 1-2 hours.')
#            pTaskID = pTask.id
#        else:
    #if we don't have a Task, how can we do lessons based on PROPKA?:
        pTaskID = int(request.POST['ptask'])
        pTask = Task.objects.get(id=pTaskID)
#        if request.POST.has_key('useqmmm'):
#            saveAtomSelections(request, ws, pTask)
        return calcpropka(request,ws,pTask)
    else:
        # get all workingFiles associated with this struct
        tdict['tasks'] = tasks #tasks was assigned before
        #we do need these tasks so we can get the right coordinates
        lesson_ok, dd_ok = checkPermissions(request)
        tdict['messages'] = get_messages(request)
        tdict['lesson_ok'] = lesson_ok
        tdict['dd_ok'] = dd_ok
        return render_to_response('html/propka_task.html', tdict)

def calcpropka(request, ws, pTask):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    tdict = {}
    #set up files for PROPKA calculation, store file for later display
    os.chdir(ws.structure.location) #else PROPKA can't get permissions right
    inp_file_location = "%s-%s.pdb" % (ws.identifier,pTask.action)
    out_file_location = inp_file_location.replace(".pdb",".pka")
    logfp = open("/tmp/propka-error.txt","w")
    try:
        subprocess.check_call(["propka31",inp_file_location],stderr=logfp)
        #PROPKA thankfully has normal exit calls
    except:
#        traceback.print_exc(file=logfp)
        messages.error(request, "PROPKA failed to process your structure. Please check your output for problems.")
    logfp.close()
    output_lines = []
    try:
        propka_file = open(out_file_location,"r")
        for line in propka_file:
            output_lines.append(line)
        propka_file.close()
        tdict['output_lines'] = output_lines
    except:
        pass
    lesson_ok, dd_ok = checkPermissions(request)
    tdict['messages'] = get_messages(request)
    tdict['lesson_ok'] = lesson_ok
    tdict['dd_ok'] = dd_ok
    return render_to_response("html/propka_finished.html",tdict) #I would redirect but we need output_lines

def display(request,filename):
    #displays PROPKA result, allows display without download of output from page
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    username = request.user.username
    tdict = {}
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please upload a structure first.")
        return HttpResponseRedirect("/charmming/fileupload/")
    try:
         ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "A working structure needs to be uploaded before PROPKA analysis can be performed on your system.")
        return HttpResponseRedirect("/charmming/buildstruct/")

    output_lines = []
    try:
        propka_file = open(("%s/%s" % (ws.structure.location,filename)),"r")
        for line in propka_file:
            output_lines.append(line)
        propka_file.close()
        tdict['output_lines']=output_lines
    except:
        logfp = open("/tmp/propka-error.txt","w")
        traceback.print_exc(file=logfp)
        logfp.close()
        messages.error(request, "Your PROPKA output could not be displayed. Perhaps the file is no longer there. Please re-do your PROPKA calculations, then try again.")
        return HttpResponseRedirect("/charmming/anaylsis/propka/")
    lesson_ok, dd_ok = checkPermissions(request)
    tdict['messages'] = get_messages(request)
    tdict['lesson_ok'] = lesson_ok
    tdict['dd_ok'] = dd_ok
    return render_to_response("html/propka_finished.html",tdict) #I would redirect but we need output_lines

def download(request, filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please upload a structure first.")
        return HttpResponseRedirect("/charmming/fileupload/")
    try:
         ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "A working structure needs to be uploaded before PROPKA analysis can be performed on your system.")
        return HttpResponseRedirect("/charmming/buildstruct/")

    mimetype,encoding = mimetypes.guess_type("%s/%s" % (struct.location,filename))

    try:
        sval = os.stat("%s/%s" % (struct.location,filename))
    except:
        messages.error("The file " + filename + " does not appear to exist anymore.")
        return HttpResponseRedirect("/charmming/fileupload/")

    response = HttpResponse(mimetype=mimetype)
    response['Content-Disposition']= 'attachment; filename=%s' % filename
    response['Content-length'] = sval.st_size
    response.write(file("%s/%s" % (struct.location,filename), "rb").read())
    return response
