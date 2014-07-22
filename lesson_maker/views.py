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
from django.contrib.messages import get_messages
import django.contrib.messages as messages
from django.template import *
from django.template import RequestContext
from account.views import checkPermissions
from mutation.models import mutateTask
import structure.models
import lesson_maker.models
import os, copy
import charmming_config
import traceback

def lessonmaker_display(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    if not request.user.is_superuser:
        messages.error(request,"You are not authorized to access this page.")
        return HttpResponseRedirect("/charmming/")
    tdict = {}
    try:
        struct=structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request,"No structure has been uploaded. Please upload a structure to begin recording lessons for the Lesson Maker.")
        return HttpResponseRedirect("/charmming/fileupload/")

    if not struct.lessonmaker_active:
        messages.error(request, "Your structure is not currently recording your actions for the Lesson Maker. Please upload a structure that does so before using the Lesson Maker.")
        return HttpResponseRedirect("/charmming/fileupload/")

    #now that we've verified that things are doing fine we can keep going
    #we should not assume that we have a workstruct or ANYTHING. Lessons can be just an upload
    #in the end. They're not very useful lessons, clearly, but we have to leave the
    #possibility open.

    recorded_lesson = lesson_maker.models.Lesson.objects.filter(structure=struct)[0]
    #everything else is attached to Structure - if it's selected, then we don't have to check again, etc.

    steps = lesson_maker.models.Step.objects.filter(lesson=recorded_lesson)
    #now we have all the steps that are in that lesson, so we can make a list of steps

    step_list = []
    #let's make a list!

    for step in steps:
        step_list.append((step.step,step.type))
    if len(step_list) < 1:
        #uh oh
        messages.error(request, "Your lesson has less than one step. This should not happen. Please report this bug to your system administrator, and try to build your structure for recording lessons again.")
        return HttpResponseRedirect("/charmming/fileupload/")
    elif len(step_list) < 3:
        messages.error(request, "Your lesson does not perform any calculations. Please perform some calculation on your structure for a proper lesson.")
        return HttpResponseRedirect("/charmming/")
    new_step_list = [] #we need to hold references to previous and next, because the template can't read them
    for i in range(0,len(step_list)):
        if i < len(step_list) - 2:
            new_step_list.append((step_list[i],step_list[i+1])) #this, next, we don't need prev
        else:
            new_step_list.append((step_list[i],None))

    tdict['step_list'] = new_step_list
    tdict['messages'] = get_messages(request)
    #now that we have a list we can make a bunch of text boxes for the user to put stuff in to
    #in theory we can display the parameters there too but that's a lot more work and
    #involves iterating through all of the class properties
    #it would be useful to display pdb results? FOr now let's make simple boxes like these
    #since we don't need much more.
    return render_to_response('html/lesson_maker_create.html',tdict)


def lessonmaker_done(request):
    #if the user puts in a lesson that already exists, we overwrite, because deleting
    #is a pain and they should deal with that, not us
    logfp = open("/tmp/create_lesson_errors.txt","w")
    objectives = []
    steps = []
    for key,value in request.POST.iteritems():
        #since we need to iterate through to get the objectives stuff anyway
        #we might as well fill in the rest of our variables here
        if "lesson_name" in key:
            lesson_name = request.POST['lesson_name']
        elif "intro_text" in key:
            intro_text = request.POST['intro_text']
        elif "objective" in key:
            objectives.append((key.split(" ")[1],value))
        elif "title" in key:
            lesson_title = request.POST['title']
        elif "done" in key or "running" in key:
            new_key = key.strip().split(" ")
            steps.append((value,float(new_key[1]),new_key[2])) #this is text, step number, status (running vs done)
    #ok, now that we did all we needed to do with the postdata, let's get our Lesson, Step and Objective objects
    try:
        curr_struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
        lesson_obj = lesson_maker.models.Lesson.objects.filter(structure=curr_struct)[0]
        step_objects = lesson_maker.models.Step.objects.filter(lesson=lesson_obj) #don't take 0, we need all of them 
    except:
        traceback.print_exc(file=logfp)
        logfp.close()
        messages.error(request, "Something has gone terribly wrong when trying to create database objects for the Lesson Maker. Please contact your system administrator.")
        return HttpResponseRedirect("/charmming/fileupload/")
    if lesson_obj.finished: #get out - you shouldn't do this twice. This object should have been wiped.
        messages.error(request, "Your current lesson recording has already been turned into a file. If you wish to create a different lesson using this same structure, please upload it again.")
        return HttpResponseRedirect("/charmming/fileupload/")


    lesson_obj.name = lesson_name
    lesson_obj.intro_text = intro_text
    lesson_obj.title = lesson_title
    lesson_obj.save() #we don't mark it as saved UNTIL it's actually done.

    for objective in objectives:
        new_obj = lesson_maker.models.Objective()
        new_obj.lesson = lesson_obj
        new_obj.obj_num = int(objective[0])
        new_obj.obj_text = objective[1]
        new_obj.save()

    for step in steps:
        try:
            new_step = step_objects.get(step=str(step[1]))
        except:
            traceback.print_exc(file=logfp)
            logfp.close()
            messages.error(request, "There is more than one step with the same number. Please report this error to your system administrator.")
            return HttpResponseRedirect("/charmming/fileupload/")
        if "running" in step[2]:
            new_step.running_text = step[0]
        elif "done" in step[2]:
            new_step.done_text = step[0]
        else: #uh oh. This is likely that someone messed with the values of the name attributes on the page.
            logfp.write(step[2])
            logfp.close()
            messages.error(request, "Invalid step. Please contact your system administrator about this message.")
            return HttpResponseRedirect("/charmming/fileupload/")
        new_step.save()

    #At this point we're done.
    lesson_obj.finished = True
    lesson_obj.save()
    logfp.close()
    messages.success(request, "Lessons Data successfully saved to the database. Please run " + charmming_config.charmming_root + "/create_lesson.py to complete the creation process.")
    return HttpResponseRedirect("/charmming/")

def lessonmaker_mockup(request):
    tdict = {}
    tdict['messages'] = get_messages(request)
    return render_to_response('html/lesson_maker_mockup.html', tdict)
