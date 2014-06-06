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
from django import forms
from django.http import HttpResponseRedirect, HttpResponse
from django.contrib import messages
from django.contrib.messages import get_messages
from django.shortcuts import render_to_response
from django.db.models import Q
from structure.models import Structure, WorkingStructure, WorkingFile, Task, CGWorkingStructure, noNscaleFound
from solvation.ionization import neutralize_tpl
from solvation.models import solvationTask
from account.views import checkPermissions
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from lessons.models import LessonProblem
import input, output, lesson1, lesson2
import re, os, copy, math
import lessonaux, charmming_config

#processes form data for minimization
def solvationformdisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)
    try:
        struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request,"Please submit a structure first.")
        return HttpResponseRedirect("/charmming/fileupload/")
#        return HttpResponse("Please submit a structure first.")
    try:
        ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "Please build a working structure before performing any calculations.")
        return HttpResponseRedirect("/charmming/buildstruct/")
#       return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    try:
        cgws = CGWorkingStructure.objects.get(workingstructure_ptr=ws.id)
        messages.error(request, "Solvation is not supported for coarse-grain models. If you wish to perform a solvation on this structure, please build a working structure with the CHARMM all-atom model.")
        return HttpResponseRedirect("/charmming/buildstruct/")
    except:
        pass

    if request.POST.has_key('solvation_structure'):
        # if there is a previous solvation structure, deactivate it
        try:
            oldtsk = solvationTask.objects.filter(workstruct=ws,active='y')[0]
            oldtsk.active = 'n'
            oldtsk.save()
        except:
            pass

        st = solvationTask()
        st.setup(ws)
        st.active = 'y'
        st.action = 'solvation'
        st.save()
      
        if ws.isBuilt != 't':
            isBuilt = False
            try:
                pTask = ws.build(st)
            except noNscaleFound, e:
                return output.returnSubmission('Minimization', error='The nScale parameterization process has not yet completed. It may take 1-2 hours.')
            pTaskID = pTask.id
        else:
            isBuilt = True
            pTaskID = int(request.POST['ptask'])

        return solvate_tpl(request,st,pTaskID)
    else:
        # get all completed tasks associated with this struct
        lesson_ok, dd_ok = checkPermissions(request)
        tasks = Task.objects.filter(workstruct=ws,status='C',active='y',modifies_coordinates=True)
#        for task in tasks:
#            if task.action == "mutation":
#                messages.error(request, "Solvation is not currently supported for mutated structures.")
#                return HttpResponseRedirect("/charmming/buildstruct/")
#
        return render_to_response('html/solvationform.html', {'ws_identifier': ws.identifier,'tasks': tasks, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})


def solvate_tpl(request,solvTask,pTaskID):
    postdata = request.POST

    workingstruct = solvTask.workstruct
    os.chdir(workingstruct.structure.location)

    solvTask.solvation_structure = postdata['solvation_structure']
    solvTask.concentration = '0.0' # needs to be set here, will be set to proper value in ionization.py
    solvTask.structure = workingstruct

    # set a, b, c, alpha, beta, gamma for run
    if postdata['solvtype'] == 'sphere_solvation':
        solvTask.spradius = postdata.spradius
    elif postdata['solvtype'] == 'set_dimensions':
        solvTask.xtl_x = postdata['set_x']
        solvTask.xtl_y = postdata['set_y']
        solvTask.xtl_z = postdata['set_z']
    elif postdata['solvtype'] == 'determine_dimensions':
        dimensions = workingstruct.dimension
        ldim = list(dimensions)
        ldim.sort()
        buffer = float(postdata['buffer'])

        if solvTask.solvation_structure in ['cubic','rhdo']:
            solvTask.xtl_x = ldim[2] + (2*buffer)
            solvTask.xtl_y = ldim[2] + (2*buffer)
            solvTask.xtl_z = ldim[2] + (2*buffer)
        elif solvTask.solvation_structure in ['tetra','hexa']:
            # the solvation input template does a coor orient, so it's safe to
            # assume that the X edge is the longest and the Y edge is the second
            # longest.
            solvTask.xtl_x = ldim[2] + (2*buffer)
            solvTask.xtl_y = ldim[2] + (2*buffer)
            solvTask.xtl_z = ldim[0] + (2*buffer)
        else:
            return output.returnSubmission('Solvation', error='Unknown crystal structure %s' % solvTask.solvation_structure)
    else:
        raise AssertionError('Unknown type')

    solvTask.save()

    if not solvTask.check_valid():
        return output.returnSubmission('Solvation', error='Invalid crystal structure definition')

    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['angles'] = "%5.2f %5.2f %5.2f" % (solvTask.angles[0],solvTask.angles[1],solvTask.angles[2])
    template_dict['topology_list'], template_dict['parameter_list'] = workingstruct.getTopparList()
    template_dict['solvation_structure'] = solvTask.solvation_structure

    pTask = Task.objects.get(id=pTaskID)
    template_dict['input_file'] = solvTask.workstruct.identifier + '-' + pTask.action
    template_dict['output_name'] = solvTask.workstruct.identifier + '-solvation'
    solvTask.parent = pTask
    solvTask.save()

    if solvTask.solvation_structure == 'sphere':
        template_dict['spradius'] = postdata['spradius']
    else:
        template_dict['xtl_x'] = float(solvTask.xtl_x)
        template_dict['xtl_y'] = float(solvTask.xtl_y)
        template_dict['xtl_z'] = float(solvTask.xtl_z)

        # use spradius as the minimum safe sphere to cut around (replaces the old
        # stuff in calc_delete.inp.
        maxl = max(template_dict['xtl_x'], template_dict['xtl_y'], template_dict['xtl_z'])
        template_dict['spradius'] = math.sqrt(3*maxl**2)/2.0

    # make sure that toppar_water_ions.str is in the TopparList
    if workingstruct.topparStream:
        if not 'toppar_water_ions.str' in workingstruct.topparStream:
            workingstruct.topparStream += ' %s/toppar/toppar_water_ions.str' % charmming_config.data_home
    else:
        workingstruct.topparStream = '%s/toppar/toppar_water_ions.str' % charmming_config.data_home

    workingstruct.save() # probably isn't necessary -- being safe
    if workingstruct.topparStream:
        template_dict['tpstream'] = workingstruct.topparStream.split()
    else:
        template_dict['tpstream'] = []

    t = get_template('%s/mytemplates/input_scripts/solvation_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    
    solvate_input_filename = workingstruct.structure.location + "/" + workingstruct.identifier + "-solvation.inp"
    inp_out = open(solvate_input_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()
    

    #change the status of the file regarding solvation
    solvTask.scripts += ',%s' % solvate_input_filename
    solvTask.save()

    doneut = postdata.has_key('salt') and postdata['salt'] != 'none'
    if doneut:
        solvTask.action = 'neutralization'
        neutralize_tpl(solvTask,postdata)
        ctask = 'Solvation'
    else:
        solvTask.start()
        solvTask.save()
        workingstruct.save() 
        ctask = 'Neutralization'

    #YP lessons stuff

    try:
        lnum=workingstruct.structure.lesson_type
        lesson_obj = eval(lnum+'.models.'+lnum.capitalize()+'()')
    except:
        lesson_obj = None


    if lesson_obj:
        lessonaux.doLessonAct(workingstruct.structure,"onSolvationSubmit",solvTask,"")
    #end YP lesson stuff
    
    return output.returnSubmission(ctask)


