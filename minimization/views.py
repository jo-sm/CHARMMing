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
from account.views import isUserTrustworthy
from structure.models import Structure, WorkingStructure, WorkingFile, Segment, goModel 
from structure.qmmm import makeQChem, makeQChem_tpl, handleLinkAtoms, writeQMheader
from structure.aux import checkNterPatch
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from minimization.models import minimizeParams
from solvation.models import solvationParams
import charmming_config, input, output, lessonaux
import re, copy
import os, shutil
import commands

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
            oldtsk = minimizeTask.objects.filter(workstruct=workstruct, selected='y')[0]
	    oldtsk.active = 'n'
	    oldtsk.save()
        except:
            pass

        mp = minimizeTask(self)
        mp.active = 'y'
        mp.workstruct = ws

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
        tasks = structure.model.Task.objects.filter(workstruct=ws,status='C',active='y')
        return render_to_response('html/minimizeform.html', {'ws_identifier': ws.identifier,'tasks': tasks})

def minimize_tpl(request,mp,pTaskID):
    postdata = request.POST

    #change the status of the file regarding minimization 
    sdsteps = postdata['sdsteps']
    abnr = postdata['abnrsteps']
    tolg = postdata['tolg']
    os.chdir(workstruct.structure.location)
    
    # create a model for the minimization
    try:
        mp.sdsteps = int(sdsteps)
        mp.abnrsteps = int(abnr)
        mp.tolg = tolg
    except:
        # FixMe: alert user to errors
        return "Error"        

    if postdata.has_key('useqmmm'):
        mp.useqmmm = 'y'
        file.checkForMaliciousCode(postdata['qmsele'],postdata)
        try:
            mp.qmmmsel = postdata['qmsele']
        except:
            pass
    else:
        mp.useqmmm = 'n'


    # template dictionary passes the needed variables to the template 
    template_dict = {}
    template_dict['topology_list'] = workstruct.getTopologyList()
    template_dict['parameter_list'] = workstruct.getParameterList()
    template_dict['output_name'] = 'mini-' + workstruct.identifier

    pTask = structure.model.Task.objects.get(id=pTaskID)
    template_dict['input_file'] = pTask.action

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
        template_dict['qmmmsel'] = mp.qmmmsel

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
            return HttpResponse('Invalid options')

        # decide if the structure we're dealing with has
        # been solvated.
        solvated = False
        wfc = pTask
        while True:
            if wfc.action == 'solvation':
                solvated = True
                break
            if wfc.parent:
                wfc = wfc.parent
            else:
                break
        if not solvated:
            return HttpResponse('Requested PBC on unsolvated structure')

        dopbc = True
        try:
            sp = solvationParams.objects.filter(structure=workstruct,selected='y')[0]
        except:
            return HttpResponse("Err ... couldn't find solvation parameters")
        template_dict['xtl_x'] = sp.xtl_x
        template_dict['xtl_y'] = sp.xtl_y
        template_dict['xtl_z'] = sp.xtl_z
        template_dict['xtl_angles'] = "%10.6f %10.6f %10.6f" % (sp.angles[0],sp.angles[1],sp.angles[2])
        template_dict['xtl_ucell'] = sp.solvation_structure
        template_dict['ewalddim'] = sp.calcEwaldDim()

        if template_dict['xtl_ucell'] == 'sphere':
            return HttpResponse('Cannot do PBC on a sphere')

    template_dict['dopbc'] = dopbc
    template_dict['useqmmm'] = postdata.has_key("useqmmm")
    template_dict['sdsteps'] = sdsteps
    template_dict['abnr'] = abnr
    template_dict['tolg'] = tolg
    if dopbc:
        mp.usepbc = 't'
    else:
        mp.usepbc = 'f'
    
    if postdata.has_key("useqmmm"):
        # validate input
	if postdata['qmmm_exchange'] in ['HF','B','B3']:
	    exch = postdata['qmmm_exchange']
	else:
	    exch = 'HF' 
	if postdata['qmmm_correlation'] in ['None','LYP']:
	    corr = postdata['qmmm_correlation']
	else:
	    corr = 'None' 
	if postdata['qmmm_basisset'] in ['STO-3G','3-21G*','6-31G*']:
	    bs = postdata['qmmm_basisset']
	else:
	    bs = 'sto3g' 
	qmsel = postdata['qmsele']
	if qmsel == '':
	    qmsel = 'resid 1'
	if postdata['qmmm_charge'] in ['-5','-4','-3','-2','-1','0','1','2','3','4','5']:
	    charge = postdata['qmmm_charge']
	else:
	    charge = '0' 
	if postdata['qmmm_multiplicity'] in ['0','1','2','3','4','5','6','7','8','9','10']:
	    multi = postdata['qmmm_multiplicity']
	else:
	    multi = '0'
	if int(postdata['num_linkatoms']) > 0:
	    linkatoms = handleLinkAtoms(file,postdata)
	else:
	    linkatoms = None
        template_dict = makeQChem_tpl(template_dict, file, exch, corr, bs, qmsel, "Force", charge, multi, file.stripDotPDB(final_pdb_name) + ".crd", linkatoms)

    template_dict['headqmatom'] = 'blankme'
    if mp.useqmmm == 'y':
        headstr = writeQMheader("", "SELE " + qmsel + " END")
        template_dict['headqmatom'] = headstr.strip() 
    mp.statusHTML = "<font color=yellow>Processing</font>"
    mp.save()
    t = get_template('%s/mytemplates/input_scripts/minimization_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    
    user_id = workstruct.structure.owner.id
    minimize_filename = workstruct.structure.location + "/minimize-" + workstruct.identifier + ".inp"
    inp_out = open(minimize_filename ,'w')
    inp_out.write(charmm_inp)
    inp_out.close()	
    scriptlist.append(minimize_filename)
    si = schedInterface()
    newJobID = si.submitJob(user_id,workstruct.structure.location,scriptlist)

    # lessons are borked under the new order, ToDo: come back and fix this
    #if file.lesson_type:
    #    lessonaux.doLessonAct(file,"onMinimizeSubmit",postdata,final_pdb_name)

    if newJobID < 0:
       mp.statusHTML = "<font color=red>Failed</font>"
       mp.save()
    else:
       workstruct.minimization_jobID = newJobID
       workstruct.save()
       sstring = si.checkStatus(newJobID)
       mp.statusHTML = statsDisplay(sstring,newJobID)
       mp.save()


    workstruct.save()
    return HttpResponse("Done.")


