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
from structure.models import Structure, WorkingStructure, Segment, goModel 
from structure.qmmm import makeQChem, makeQChem_tpl, handleLinkAtoms, writeQMheader
from structure.editscripts import generateHTMLScriptEdit
from structure.aux import checkNterPatch
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from minimization.models import minimizeParams
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
        scriptlist = []
        if ws.isBuilt != 't':
            ws.build(scriptlist)
        return minimize_tpl(request,ws,scriptlist)
    else:
        return render_to_response('html/minimizeform.html', {'ws_identifier': ws.identifier})

def minimize_tpl(request,workstruct,scriptlist):
    postdata = request.POST
    #deals with changing the selected minimize_params
    try:
        oldparam = minimizeParams.objects.filter(pdb = file, selected = 'y')[0]
	oldparam.selected = 'n'
	oldparam.save()
    except:
        pass

    #change the status of the file regarding minimization 
    sdsteps = postdata['sdsteps']
    abnr = postdata['abnrsteps']
    tolg = postdata['tolg']
    os.chdir(workstruct.structure.location)
    
    try:
        selectedparam = minimizeParams.objects.filter(struct=workstruct,selected='y')[0]
        selectedparam.selected = 'n'
	selectedparam.save()
    except:
        #No selected minimizeparam
	pass

    # create a model for the minimization
    mp = minimizeParams(selected='y')
    try:
        mp.sdsteps = int(sdsteps)
        mp.abnrsteps = int(abnr)
        mp.tolg = tolg
    except:
        # FixMe: alert user to errors
        return "Error"        

    try:
        if postdata['usepbc']:
            mp.usepbc = 'y'
        else:
            mp.usepbc = 'n'        
    except:
        mp.usepbc = 'n'

    if postdata.has_key('useqmmm'):
        mp.useqmmm = 'y'
        file.checkForMaliciousCode(postdata['qmsele'],postdata)
        try:
            mp.qmmmsel = postdata['qmsele']
        except:
            pass
    else:
        mp.useqmmm = 'n'

    mp.save()

    # template dictionary passes the needed variables to the template 
    template_dict = {}
    template_dict['topology_list'] = workstruct.getTopologyList()
    template_dict['parameter_list'] = workstruct.getParameterList()
    template_dict['filebase'] = workstruct.identifier    

    # ack, this is totally wrong, ToDo, FIXME
    template_dict['restraints'] = '' 
    try:
        postdata['apply_restraints']
        template_dict['restraints'] = workstruct.handleRestraints(request)
    except:
        pass
    
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
    template_dict['usepbc'] = ''
    template_dict['dopbc'] = 0
    try:
        if postdata['usepbc']:
            template_dict['usepbc'] =  postdata['usepbc']
            if file.solvation_structure != 'sphere':
                dopbc = 1
                template_dict['dopbc'] = 1
            else:
                dopbc = 0
        else:
            dopbc = 0
    except:
        dopbc = 0

    template_dict['solvation_structure'] = file.solvation_structure
    template_dict['relative_boundary'] = 0
    if dopbc:
        relative_boundary = 0
        if file.solvation_structure != '' and file.crystal_x < 0:
            relative_boundary = 1

        template_dict['relative_boundary'] = relative_boundary  
        template_dict['dim_x'] = str(file.crystal_x)
        template_dict['dim_z'] = str(file.crystal_z)
        # Tim Miller test
        greaterval = max(file.crystal_x,file.crystal_z)
        template_dict['greaterval'] = str(greaterval)

        # set up images
        if file.solvation_structure == '' or solvate_implicitly:
            pass  
        else:
            # we should have a solvation file to read from
            try:
                os.stat(file.location + "new_" + file.stripDotPDB(file.filename) + ".xtl")
            except:
                # need to throw some sort of error ... for now just toss cookies
                return HttpResponse("Oops ... transfer file not found.")

        # we need to get the ewald parameter
        template_dict['file_location'] = file.location
    template_dict['useqmmm'] = postdata.has_key("useqmmm")
    template_dict['sdsteps'] = sdsteps
    template_dict['abnr'] = abnr
    template_dict['tolg'] = tolg
    
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
    mp.pdb = file
    mp.save()
    file.save()
    t = get_template('%s/mytemplates/input_scripts/minimization_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    
    user_id = file.owner.id
    minimize_filename = file.location + "/" + file.name + "/minimize.inp"
    inp_out = open(minimize_filename ,'w')
    inp_out.write(charmm_inp)
    inp_out.close()	
    scriptlist.append(minimize_filename)
    file.save()
    if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
        return generateHTMLScriptEdit(charmm_inp,scriptlist,'minimization')
    else:
        si = schedInterface()
        newJobID = si.submitJob(user_id,file.location,scriptlist)
	if file.lesson_type:
            lessonaux.doLessonAct(file,"onMinimizeSubmit",postdata,final_pdb_name)
        file.save()

        if newJobID < 0:
           mp.statusHTML = "<font color=red>Failed</font>"
           mp.save()
        else:
           file.minimization_jobID = newJobID
           sstring = si.checkStatus(newJobID)
           mp.statusHTML = statsDisplay(sstring,newJobID)
           mp.save()
           file.save()   
	return "Done."


