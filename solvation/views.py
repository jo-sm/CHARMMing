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
from django.shortcuts import render_to_response
from structure.models import Structure, WorkingStructure, WorkingFile
from solvation.ionization import neutralize_tpl
from solvation.models import solvationParams
from account.views import isUserTrustworthy
from structure.editscripts import generateHTMLScriptEdit
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
import input, output
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
        return HttpResponse("Please submit a structure first.")
    try:
        ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
       return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    if request.POST.has_key('solvation_structure'):
        scriptlist = []
        if ws.isBuilt != 't':
            isBuilt = False
            pstruct = ws.build(scriptlist)
            pstructID = pstruct.id
        else:
            isBuilt = True
            pstructID = int(request.POST['pstruct'])

        return solvate_tpl(request,ws,isBuilt,pstructID,scriptlist)
    else:
        # get all workingFiles associated with this struct
        wfs = WorkingFile.objects.filter(structure=ws,type='crd')
        return render_to_response('html/solvationform.html', {'ws_identifier': ws.identifier,'workfiles': wfs})


def solvate_tpl(request,workingstruct,isBuilt,pstructID,scriptlist):
    postdata = request.POST
    #deals with changing the selected minimize_params
    try:
        oldparam = solvationParams.objects.filter(struct=ws, selected='y')[0]
        oldparam.selected = 'n'
        oldparam.save()
    except:
        pass

    os.chdir(workingstruct.structure.location)

    sp = solvationParams()
    sp.selected = 'y'    
    sp.pdb = file
    sp.statusHTML = "<font color=yellow>Processing</font>"

    #file.solvation_status = "<font color=yellow>Processing</font>"
    sp.solvation_structure = postdata['solvation_structure']
    sp.concentration = '0.0' # needs to be set here, will be set to proper value in ionization.py
    sp.structure = workingstruct

    # set a, b, c, alpha, beta, gamma for run
    if sp.solvation_structure == 'sphere':
        sp.spradius = postdata.sphere_radius
    else:
        sp.xtl_x = postdata['set_x']
        sp.xtl_y = postdata['set_y']
        sp.xtl_z = postdata['set_z']
    sp.save()

    if not sp.check_valid():
        return HttpResponse("Invalid crystal structure definition")

    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['angles'] = "%5.2f %5.2f %5.2f" % (sp.angles[0],sp.angles[1],sp.angles[2])
    template_dict['topology_list'] = workingstruct.getTopologyList()
    template_dict['parameter_list'] = workingstruct.getParameterList()
    template_dict['solvation_structure'] = sp.solvation_structure

    pstruct = WorkingFile.objects.filter(id=pstructID)[0]
    template_dict['input_file'] = pstruct.basename
    template_dict['output_name'] = 'solv-' + workingstruct.identifier
    sp.inpStruct = pstruct
    sp.save()

    if sp.solvation_structure == 'sphere':
        template_dict['spradius'] = postdata['spradius']
    else:
        template_dict['xtl_x'] = float(sp.xtl_x)
        template_dict['xtl_y'] = float(sp.xtl_y)
        template_dict['xtl_z'] = float(sp.xtl_z)

        # use spradius as the minimum safe sphere to cut around (replaces the old
        # stuff in calc_delete.inp.
        maxl = max(template_dict['xtl_x'], template_dict['xtl_y'], template_dict['xtl_z'])
        template_dict['spradius'] = math.sqrt(3*maxl**2)/2.0

    t = get_template('%s/mytemplates/input_scripts/solvation_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    
    user_id = workingstruct.structure.owner.id
    solvate_input_filename = workingstruct.structure.location + "/solvate-" + workingstruct.identifier + ".inp"
    inp_out = open(solvate_input_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()

    #change the status of the file regarding solvation
    scriptlist.append(solvate_input_filename)
    workingstruct.save()

    doneut = postdata.has_key('salt') and postdata['salt'] != 'none'
    if doneut:
        neutralize_tpl(workingstruct,sp,postdata,scriptlist)
    else:
        si = schedInterface()
        newJobID = si.submitJob(user_id,workingstruct.structure.location,scriptlist)
        # Lessons are borked at the moment...
        #if file.lesson_type:
        #    lessonaux.doLessonAct(file,"onSolvationSubmit",postdata)
        workingstruct.solvation_jobID = newJobID
        sstring = si.checkStatus(newJobID)
        sp.statusHTML = statsDisplay(sstring,newJobID)
        sp.save()
        workingstruct.save() 

    return HttpResponse("Done")


