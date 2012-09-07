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
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from account.views import checkPermissions
from structure.models import Structure, WorkingStructure, WorkingFile, Segment, goModel, Task
from structure.qmmm import makeQChem_tpl, makeQchem_val
from structure.aux import checkNterPatch
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from minimization.models import minimizeTask
from solvation.models import solvationTask
import charmming_config, input, output, lessonaux, lessons, lesson1, lesson2, lesson3, lesson4
import re, copy
import os, shutil
import commands, traceback

#processes form data for minimization
def minimizeformdisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)
    #chooses the file based on if it is selected or not
    try:
        struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    try:
       ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
       return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    if request.POST.has_key('sdsteps') or request.POST.has_key('abnrsteps'):
        #deals with changing the selected minimize_params
        try:
            oldtsk = minimizeTask.objects.filter(workstruct=ws,active='y')[0]
            oldtsk.active = 'n'
            oldtsk.save()
        except:
            pass

        mp = minimizeTask()
        mp.sdsteps = 0
        mp.abnrsteps = 0
        mp.tolg = 0.01
        mp.setup(ws)
        mp.active = 'y'
        mp.action = 'minimization'
        mp.save()

        if ws.isBuilt != 't':
            isBuilt = False
            pTask = ws.build(mp)
            pTaskID = pTask.id
        else:
            isBuilt = True
            pTaskID = int(request.POST['ptask'])
                
        return minimize_tpl(request,mp,pTaskID)
    else:
        # get all workingFiles associated with this struct
        tasks = Task.objects.filter(workstruct=ws,status='C',active='y').exclude(action='energy')

        # segments are also needed for QM/MM
        seg_list = []
        for seg in ws.segments.all():
            seg_list.append(seg.name)       

        lesson_ok, dd_ok = checkPermissions(request)
        return render_to_response('html/minimizeform.html', {'ws_identifier': ws.identifier,'tasks': tasks, 'seg_list': seg_list, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

def minimize_tpl(request,mp,pTaskID):
    postdata = request.POST

    #change the status of the file regarding minimization 
    sdsteps = postdata['sdsteps']
    abnr = postdata['abnrsteps']
    tolg = postdata['tolg']
    os.chdir(mp.workstruct.structure.location)
    
    # fill in the model for the minimization
    try:
        mp.sdsteps = int(sdsteps)
        mp.abnrsteps = int(abnr)
        mp.tolg = tolg
    except:
        return output.returnSubmission('Minimization',error='Form not filled out.')        

    if mp.sdsteps > charmming_config.minimize_steplimit or mp.abnrsteps > charmming_config.minimize_steplimit:
        return output.returnSubmission('Minimization',error='Minimization has a limit of %d steps.' % charmming_config.minimize_steplimit)

    if postdata.has_key('useqmmm'):
        mp.useqmmm = 'y'
        input.checkForMaliciousCode(postdata['qmsele'],postdata)
        try:
            mp.qmmmsel = postdata['qmsele']
        except:
            pass
    else:
        mp.useqmmm = 'n'


    # template dictionary passes the needed variables to the template 
    template_dict = {}
    template_dict['topology_list'], template_dict['parameter_list'], junk = mp.workstruct.getTopparList()
    if mp.workstruct.topparStream:
        template_dict['tpstream'] = mp.workstruct.topparStream.split()
    else:
        template_dict['tpstream'] = []
    template_dict['output_name'] = mp.workstruct.identifier + '-minimization'

    pTask = Task.objects.filter(id=pTaskID)[0]
    template_dict['input_file'] = mp.workstruct.identifier + '-' + pTask.action

    mp.parent = pTask
    mp.save()

    solvate_implicitly = 0
    try:
        if(postdata['solvate_implicitly']):
            solvate_implicitly = 1
    except:
        pass

    template_dict['solvate_implicitly'] = solvate_implicitly
    if postdata.has_key('fixnonh'):
        template_dict['fixnonh'] = 1
    else:
        template_dict['fixnonh'] = 0

    # handles shake
    template_dict['shake'] = request.POST.has_key('apply_shake')
    if request.POST.has_key('apply_shake'):
        template_dict['which_shake'] = postdata['which_shake']
        ### WTF is this here? --> template_dict['qmmmsel'] = mp.qmmmsel

        if postdata['which_shake'] == 'define_shake':
            template_dict['shake_line'] = postdata['shake_line']
            if postdata['shake_line'] != '':
                file.checkForMaliciousCode(postdata['shake_line'],postdata)

    template_dict['restraints'] = ''
    try:
        postdata['apply_restraints']
        template_dict['restraints'] = file.handleRestraints(request)
    except:
        pass

    # check to see if PBC needs to be used -- if so we have to set up Ewald
    dopbc = False
    if request.POST.has_key('usepbc'):
        if solvate_implicitly:
            return output.returnSubmission('Minimization', error='Conflicting solvation options given')

        # decide if the structure we're dealing with has
        # been solvated.
        solvated = False
        wfc = pTask
        while True:
            if wfc.action == 'solvation' or wfc.action == 'neutralization':
                solvated = True
                break
            if wfc.parent:
                wfc = wfc.parent
            else:
                break
        if not solvated:
            return output.returnSubmission('Minimization', error='Requested PBC on unsolvated structure')

        if mp.useqmmm == 'y':
            return output.returnSubmission('Minimization', error='You cannot use periodic boundary conditions with QM/MM')

        dopbc = True
        try:
            sp = solvationTask.objects.get(workstruct=mp.workstruct,active='y')
        except:
            return output.returnSubmission('Minimization',error='You tried to use PBC on an un-solvated system.')        
        template_dict['xtl_x'] = sp.xtl_x
        template_dict['xtl_y'] = sp.xtl_y
        template_dict['xtl_z'] = sp.xtl_z
        template_dict['xtl_angles'] = "%10.6f %10.6f %10.6f" % (sp.angles[0],sp.angles[1],sp.angles[2])
        template_dict['xtl_ucell'] = sp.solvation_structure
        template_dict['ewalddim'] = sp.calcEwaldDim()

        if template_dict['xtl_ucell'] == 'sphere':
            return output.returnSubmission('Minimization',error='You tried to use PBC with a spherivcal water shell.')

    template_dict['dopbc'] = dopbc
    template_dict['useqmmm'] = mp.useqmmm == "y"
    template_dict['sdsteps'] = sdsteps
    template_dict['abnr'] = abnr
    template_dict['tolg'] = tolg

    if dopbc:
        mp.usepbc = 't'
    else:
        mp.usepbc = 'f'
    
    if mp.useqmmm == 'y':
        qmparams = makeQchem_val(postdata,mp.qmmmsel)
        qmparams['jobtype'] = 'Force'
        template_dict = makeQChem_tpl(template_dict,qmparams,mp.workstruct)

    ## There is probably a better way to do this...
    #template_dict['headqmatom'] = 'blankme'
    #if mp.useqmmm == 'y':
    #    headstr = writeQMheader("", "SELE " + qmsel + " END")
    #    template_dict['headqmatom'] = headstr.strip() 
    mp.save()
    t = get_template('%s/mytemplates/input_scripts/minimization_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    
    user_id = mp.workstruct.structure.owner.id
    minimize_filename = mp.workstruct.structure.location + "/" + mp.workstruct.identifier + "-minimize.inp"
    inp_out = open(minimize_filename ,'w')
    inp_out.write(charmm_inp)
    inp_out.close()	
    mp.scripts += ',%s' % minimize_filename
    mp.start()
    mp.save()

    #YP lessons status update
    try:
        lnum=mp.workstruct.structure.lesson_type
        lesson_obj = eval(lnum+'.models.'+lnum.capitalize()+'()')
    except:
        lesson_obj = None

    if lesson_obj:
        lessonaux.doLessonAct(mp.workstruct.structure,"onMinimizeSubmit",mp,"")
    #YP

    mp.workstruct.save()
    return output.returnSubmission('Minimization')


