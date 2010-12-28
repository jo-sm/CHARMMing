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
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.template.loader import get_template
from pdbinfo.models import PDBFile, PDBFileForm 
from account.views import isUserTrustworthy
from pdbinfo.editscripts import generateHTMLScriptEdit
from pdbinfo.aux import checkNterPatch
from minimization.views import append_tpl
from dynamics.models import mdParams,ldParams,sgldParams
import rexModule
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
import django.forms
import re
import copy, shutil
import os
import lessonaux, output, charmming_config

#processes form data for ld simulations
def lddisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    PDBFile().checkRequestData(request)
    #chooses the file based on if it is selected or not
    try:
        file =  PDBFile.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    os.chdir(file.location)


    #creates a list of filenames associated with the PDB
    #The "md" option takes away all -md,-ld, and -sgld names
    #filename_list = file.getLimitedFileList("md")
    filename_list = file.getLimitedFileList('blank')
    solv_pdb = "new_" + file.stripDotPDB(file.filename) + "-solv.pdb"
    neu_pdb = "new_" + file.stripDotPDB(file.filename) + "-neutralized.pdb"
    min_pdb = "new_" + file.stripDotPDB(file.filename) + "-min.pdb"
    md_pdb = "new_" + file.stripDotPDB(file.filename) + "-md.pdb"
    ld_pdb = "new_" + file.stripDotPDB(file.filename) + "-ld.pdb"
    sgld_pdb = "new_" + file.stripDotPDB(file.filename) + "-sgld.pdb"
    het_list = file.getNonGoodHetPDBList()
    tip_list = file.getGoodHetPDBList()
    seg2_list = file.getProteinSegPDBList()
    protein_list = [x for x in seg2_list if x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb")]
    go_list = [x for x in seg2_list if x.endswith("-go.pdb") or x.endswith("-go-final.pdb")]
    bln_list = [x for x in seg2_list if x.endswith("-bln.pdb") or x.endswith("-bln-final.pdb")]
    nucleic_list = [x for x in seg2_list if x.endswith("-rna.pdb") or x.endswith("-dna.pdb") or x.endswith("-rna-final.pdb") or x.endswith("-dna-final.pdb")]
    userup_list = [x for x in seg2_list if not ( x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb") or x.endswith("-dna.pdb") or x.endswith("-dna-final.pdb") \
                     or x.endswith("-rna.pdb") or x.endswith("-rna-final.pdb") or x.endswith("het.pdb") or x.endswith("het-final.pdb") or x.endswith("-go.pdb") \
                     or x.endswith("-go-final.pdb") or x.endswith("-bln.pdb") or x.endswith("-bln-final.pdb") ) ]
    seg_list = file.segids.split(' ')
    disulfide_list = file.getPDBDisulfides()

    # Tim Miller: 03-09-2009 we need to pass a segpatch_list (containing
    # both the list of segments and their default NTER patches) to the
    # template.
    propatch_list = []
    nucpatch_list = []
    for seg in seg_list:
        if seg.endswith("-pro"):
            defpatch = checkNterPatch(file,seg)
            propatch_list.append((seg,defpatch))
        elif seg.endswith("-dna") or seg.endswith("-rna"):
            nucpatch_list.append((seg,"5TER","3TER"))

    #Each checkbox or radio button on the ldform.html page has a name equal
    #to the filename it's associated with. This for loop tries to get the post
    #data using the name as a key
    tempid = None
    for i in range(len(filename_list)):
        try:
            tempid = request.POST[filename_list[i]]
            break
        except:
            try:
                tempid = request.POST['solv_or_min']
                break
            except:
                tempid = None

    if(tempid):
        # if the name is solv/min pdb, it does not need to be appended/patched
        # otherwise reappend them
        scriptlist = []
        if tempid != min_pdb and tempid != solv_pdb and tempid != md_pdb and tempid != ld_pdb and tempid != sgld_pdb and tempid != neu_pdb:
            seg_list = append_tpl(request.POST,filename_list,file,scriptlist)
            ld_this_file = 'new_' + file.stripDotPDB(file.filename) + "-final.pdb"
            html = applyld_tpl(request,file,seg_list,ld_this_file,scriptlist)
            return HttpResponse(html)
	else:
            html = applyld_tpl(request,file,seg_list,tempid,scriptlist)
            return HttpResponse(html)

    doCustomShake = 1
    if file.ifExistsRtfPrm() < 0:
        doCustomShake = 0

    if file.natom > 30000:
        solv_time_est = "at least 4 hours"
        unsolv_time_est = "at least half an hour"
    elif file.natom > 15000:
        solv_time_est = "between 2 to 4 hours"
        unsolv_time_est = "between 15 and 30 minutes"
    elif file.natom > 10000:
        solv_time_est = "from 1 to 2 hours"
        unsolv_time_est = "between 5 and 15 minutes"
    elif file.natom > 5000:
        solv_time_est = "from 30 minutes to an hour"
        unsolv_time_est = "from 1 to 5 minutes"    
    elif file.natom > 1000:
        solv_time_est = "from 15 to 30 minutes"
        unsolv_time_est = "less than a minute" 
    else:
        solv_time_est = "less than 15 minutes"
        unsolv_time_est = "less than a minute" 
    trusted = isUserTrustworthy(request.user)
    return render_to_response('html/ldform.html', {'filename_list': filename_list,'min_pdb':min_pdb,'solv_pdb':solv_pdb, 'neu_pdb': neu_pdb, 'md_pdb': md_pdb,'ld_pdb': ld_pdb, \
                                                   'sgld_pdb': sgld_pdb,'seg_list':seg_list,'protein_list':protein_list,'tip_list':tip_list,'het_list':het_list, 'file': file, \
                                                   'docustshake': doCustomShake, 'solv_time_est': solv_time_est, 'unsolv_time_est': unsolv_time_est,'trusted':trusted, \
                                                   'disulfide_list': disulfide_list, 'propatch_list': propatch_list, 'nucpatch_list': nucpatch_list, 'nucleic_list': nucleic_list, \
                                                   'userup_list': userup_list, 'go_list': go_list, 'bln_list': bln_list})

#processes form data for md simulations
def mddisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    PDBFile().checkRequestData(request)
    #chooses the file based on if it is selected or not
    try:
        file =  PDBFile.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    os.chdir(file.location)
    temp = file.stripDotPDB(file.filename)
    #creates a list of filenames associated with the PDB
    filename_list = file.getLimitedFileList("blank")
    #filename_list = file.getFileList()
    solv_pdb = "new_" + file.stripDotPDB(file.filename) + "-solv.pdb"
    neu_pdb = "new_" + file.stripDotPDB(file.filename) + "-neutralized.pdb"
    min_pdb = "new_" + file.stripDotPDB(file.filename) + "-min.pdb"
    md_pdb = "new_" + file.stripDotPDB(file.filename) + "-md.pdb"
    ld_pdb = "new_" + file.stripDotPDB(file.filename) + "-ld.pdb"
    sgld_pdb = "new_" + file.stripDotPDB(file.filename) + "-sgld.pdb"
    seg_list = file.segids.split(' ')
    het_list = file.getNonGoodHetPDBList()
    tip_list = file.getGoodHetPDBList()
    disulfide_list = file.getPDBDisulfides()
    seg2_list = file.getProteinSegPDBList()
    protein_list = [x for x in seg2_list if x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb")]
    go_list = [x for x in seg2_list if x.endswith("-go.pdb") or x.endswith("-go-final.pdb")]
    bln_list = [x for x in seg2_list if x.endswith("-bln.pdb") or x.endswith("-bln-final.pdb")]
    nucleic_list = [x for x in seg2_list if x.endswith("-dna.pdb") or x.endswith("-dna-final.pdb") or x.endswith("-rna.pdb") or x.endswith("-rna-final.pdb")]
    userup_list = [x for x in seg2_list if not ( x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb") or x.endswith("-dna.pdb") or x.endswith("-dna-final.pdb") \
                     or x.endswith("-rna.pdb") or x.endswith("-rna-final.pdb") or x.endswith("het.pdb") or x.endswith("het-final.pdb") or x.endswith("-go.pdb") \
                     or x.endswith("-go-final.pdb") or x.endswith("-bln.pdb") or x.endswith("-bln-final.pdb") ) ]
    #filename_list = protein_list + tip_list + het_list
    done = re.compile('Done')
    atom_counts = {}

    # Tim Miller: 03-09-2009 we need to pass a segpatch_list (containing
    # both the list of segments and their default NTER patches) to the
    # template.
    propatch_list = []
    nucpatch_list = []
    for seg in seg_list:
        if seg.endswith("-pro"):
            defpatch = checkNterPatch(file,seg)
            propatch_list.append((seg,defpatch))
        elif seg.endswith("-dna") or seg.endswith("-rna"):
            nucpatch_list.append((seg,"5TER","3TER"))


    for i in range(len(filename_list)):
        #First check and see if the selected choices are segids
	#Otherwise see if it is a solvated/minimized PDB
        try:
            tempid = request.POST[filename_list[i]]
	    filename = filename_list[i]
        except:
            try:
                tempid = request.POST['solv_or_min']
		filename = request.POST['solv_or_min']
            except:
                tempid = "null"
        if(tempid!="null"):
            seg_list = file.segids.split(' ')
            try:
                if(request.POST['usepatch']):
                    file.handlePatching(request.POST)
            except:
                #If there is no patch, make sure patch_name is zero
                file.patch_name = ""
                file.save()
            scriptlist = []
            exelist    = []
            nproclist  = []
	    if(filename != min_pdb and filename != solv_pdb and filename != md_pdb and filename != ld_pdb and filename != sgld_pdb and filename != neu_pdb):
	        seg_list = append_tpl(request.POST,filename_list,file,scriptlist)
		md_this_file = 'new_' + file.stripDotPDB(file.filename) + "-final.pdb"
                html = applymd_tpl(request,file,seg_list,md_this_file,scriptlist)
                return HttpResponse(html)
	    else:
                html = applymd_tpl(request,file,seg_list,filename,scriptlist)
                return HttpResponse(html)

    atomCounts = {}
    for seg in seg_list:
        atomCounts[seg] = file.countAtomsInSeg(seg)
    atomCounts['solv'] = file.countSolvatedAtoms()

    doCustomShake = 1
    if file.ifExistsRtfPrm() < 0:
        doCustomShake = 0

    if file.natom > 30000:
        solv_time_est = "at least 4 hours"
        unsolv_time_est = "at least half an hour"
    elif file.natom > 15000:
        solv_time_est = "between 2 to 4 hours"
        unsolv_time_est = "between 15 and 30 minutes"
    elif file.natom > 10000:
        solv_time_est = "from 1 to 2 hours"
        unsolv_time_est = "between 5 and 15 minutes"
    elif file.natom > 5000:
        solv_time_est = "from 30 minutes to an hour"
        unsolv_time_est = "from 1 to 5 minutes"    
    elif file.natom > 1000:
        solv_time_est = "from 15 to 30 minutes"
        unsolv_time_est = "less than a minute" 
    else:
        solv_time_est = "less than 15 minutes"
        unsolv_time_est = "less than a minute" 
    trusted = isUserTrustworthy(request.user)

    # decide whether or not this run is restartable (check for non-empty restart file)
    try:
        os.stat(file.location + "new_" + file.stripDotPDB(file.filename) + "-md.res")
        canrestart = True
    except:
        canrestart = False

    return render_to_response('html/mdform.html', {'filename_list': filename_list,'seg_list':seg_list,'solv_pdb':solv_pdb,'neu_pdb':neu_pdb,'min_pdb':min_pdb,'md_pdb':md_pdb, \
                                                   'ld_pdb':ld_pdb,'sgld_pdb':sgld_pdb,'het_list':het_list,'tip_list':tip_list,'protein_list':protein_list, 'file': file, \
                                                   'docustshake': doCustomShake, 'solv_time_est': solv_time_est, 'unsolv_time_est': unsolv_time_est, 'atomCounts': atomCounts, \
                                                   'trusted':trusted, 'disulfide_list': disulfide_list, 'propatch_list': propatch_list, 'nucpatch_list': nucpatch_list, \
                                                   'nucleic_list': nucleic_list, 'userup_list': userup_list, 'go_list': go_list, 'bln_list': bln_list, 'canrestart': canrestart})

def applyld_tpl(request,file,seg_list,min_pdb,scriptlist):
    postdata = request.POST
    rtf_prm_dict = file.getRtfPrmPath()
    fbeta = postdata['fbeta']   
    nstep = postdata['nstep']
    if postdata.has_key('usesgld'):
        usesgld = postdata['usesgld']
	ld_suffix = '-sgld'
        try:
            oldparam = sgldParams.objects.filter(pdb=file, selected='y')[0]
            oldparam.selected = 'n'
            oldparam.save()
        except:
            pass
	ldp = sgldParams(selected='y',pdb=file)
    else:
        usesgld = None
	ld_suffix = '-ld'
        try:
            oldparam = ldParams.objects.filter(pdb = file, selected = 'y')[0]
            oldparam.selected = 'n'
            oldparam.save()
        except:
            pass
	ldp = ldParams(selected='y',pdb=file)
    ldp.fbeta = str(fbeta)
    ldp.nstep = str(nstep)
    
    try:
        make_movie = postdata['make_movie']
        if usesgld:
           file.sgld_movie_req = True
           ldp.sgld_movie_req = True
        else:        
           file.ld_movie_req = True
           ldp.ld_movie_req = True
    except:
        make_movie = None
        if usesgld:
           file.sgld_movie_req = False
        else:
           file.ld_movie_req = False
    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'] = file.getTopologyList()
    template_dict['parameter_list'] = file.getParameterList()
    template_dict['filebase'] = file.stripDotPDB(file.filename)
    template_dict['input_file'] = file.stripDotPDB(min_pdb)
    template_dict['useqmmm'] = ''
    template_dict['qmmmsel'] = ''
    template_dict['headqmatom'] = 'blankme'
    template_dict['restraints'] = ''
    template_dict['seg_list'] = seg_list
    try:
        postdata['apply_restraints']
        template_dict['restraints'] = file.handleRestraints(request)
    except:
        pass

    #If the user wants to solvate implicitly the scpism line is needed
    #84 will be the scpism number in this program
    template_dict['solvate_implicitly'] = False
    try:
        template_dict['solvate_implicitly'] = postdata['solvate_implicitly']
        if(postdata['solvate_implicitly']):
	    ldp.scipism = True
    except:
        pass

    template_dict['fbeta'] = fbeta
    template_dict['ld_suffix'] = ld_suffix
    template_dict['shake'] = request.POST.has_key('apply_shake')
    if request.POST.has_key('apply_shake'):
        template_dict['which_shake'] = postdata['which_shake']
        if postdata['which_shake'] == 'define_shake':
            template_dict['shake_line'] = postdata['shake_line']
            if postdata['shake_line'] != '':
                file.checkForMaliciousCode(postdata['shake_line'],postdata)
    template_dict['rtfprm'] = False
    if file.ifExistsRtfPrm() < 0:
        template_dict['rtfprm'] = True
    template_dict['nstep'] = nstep
    template_dict['usesgld'] = usesgld
    if(usesgld):
        tsgavg = '0.0'
        tempsg = '0.0'

        tsgavg = postdata['tsgavg']
	ldp.tsgavg = str(tsgavg)
	tempsg = postdata['tempsg']
	ldp.tempsg = str(tempsg)

        template_dict['tsgavg'] = tsgavg
        template_dict['tempsg'] = tempsg
    template_dict['output_name'] = "new_" + file.stripDotPDB(file.filename) + ld_suffix
    user_id = file.owner.id
    os.chdir(file.location)
    ld_filename = "charmm-" + file.stripNew(file.stripDotPDB(file.filename)) +\
                  ld_suffix + ".inp"
    t = get_template('%s/mytemplates/input_scripts/applyld_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))

    inp_out = open(file.location + ld_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()  
    #change the status of the file regarding minimization
    if usesgld:
        #file.sgld_status = "<font color=yellow>Processing</font>"
        ldp.statusHTML = "<font color=yellow>Processing</font>"
    else:
        #file.ld_status = "<font color=yellow>Processing</font>"
        ldp.statusHTML = "<font color=yellow>Processing</font>"
    ldp.save()
    si = schedInterface()
    scriptlist.append(file.location + ld_filename)

    if file.lesson_type:
        if usesgld:
            lessonaux.doLessonAct(file,"onSGLDSubmit",postdata,None)
        else:
            lessonaux.doLessonAct(file,"onLDSubmit",postdata,None)

    if make_movie:
        if usesgld:
            print "Make sgld dyna movie"
            return makeJmolMovie(file,postdata,min_pdb,scriptlist,'sgld')
        else:
            return makeJmolMovie(file,postdata,min_pdb,scriptlist,'ld')
    else:
        si = schedInterface()
        if usesgld:
	    #For user editable scripts
            if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
                return generateHTMLScriptEdit(charmm_inp,scriptlist,'sgld')
	    else:
                newJobID = si.submitJob(user_id,file.location,scriptlist,{},{})
                file.sgld_jobID = newJobID
                sstring = si.checkStatus(newJobID)
                file.sgld_movie_status = None
                
                ldp.statusHTML = statsDisplay(sstring,newJobID)
                ldp.sgld_movie_status = None
                ldp.save()
        else:
	    #for user editable scripts
            if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
                return generateHTMLScriptEdit(charmm_inp,scriptlist,'ld')
	    else:
                newJobID = si.submitJob(user_id,file.location,scriptlist,{},{})
                file.ld_jobID = newJobID
                sstring = si.checkStatus(newJobID)
                file.ld_movie_status = None

                ldp.statusHTML = statsDisplay(sstring,newJobID)
                ldp.ld_movie_status = None
                ldp.save()
    file.save() 
    return "Done."

#Generates MD script and runs it
def applymd_tpl(request,file,seg_list,min_pdb,scriptlist):
    postdata = request.POST

    #deals with changing the selected minimize_params
    seqno = 1
    try:
        oldparam = mdParams.objects.filter(pdb=file, selected='y')[0]
        oldparam.selected = 'n'
        seqno = oldparam.sequence + 1
        oldparam.save()
    except:
        pass

    mdp = mdParams(selected='y',pdb=file)
    mdp.sequence = seqno
    try:
       make_movie = postdata['make_movie']
       file.md_movie_req = True
       mdp.md_movie_req = True
    except:
       make_movie = False
       file.md_movie_req = False
       mdp.md_movie_req = False
    try:
       solvate_implicitly = postdata['solvate_implicitly']
       mdp.scipism = True
    except:
       solvate_implicitly = 0
       mdp.scipism = False
    
    #check for replica exchange
    try:
        postdata['use_replica_exchange']
        useReplicaExchange = True
    except:
        useReplicaExchange = False

    # template dictionary passes the needed variables to the templat
    template_dict = {}     
    template_dict['topology_list'] = file.getTopologyList()
    template_dict['parameter_list'] = file.getParameterList()
    template_dict['input_file'] = file.stripDotPDB(min_pdb)
    template_dict['restraints'] = ''
    template_dict['solvate_implicitly'] = solvate_implicitly
    template_dict['useqmmm'] = ''
    template_dict['qmmmsel'] = ''
    template_dict['headqmatom'] = 'blankme'
    template_dict['choose_me'] = 'y'

    orig_rst = file.location + "new_" + file.stripDotPDB(file.filename) + "-md.res"
    new_rst  = file.location + "new_" + file.stripDotPDB(file.filename) + "-md-old.res"
    if postdata.has_key('dorestart'):
        template_dict['strtword'] = "restart"
        try:
            shutil.copy(orig_rst,new_rst)
        except:
            pass
    else:
        template_dict['strtword'] = "start"

    try:
        postdata['apply_restraints']
        template_dict['restraints'] = file.handleRestraints(request)
    except:
        pass

    # if we have a solvated structure that with crystal_x < 0 or we're
    # doing vacuum or scpism with dynamics, then we need to figure out the structure
    # diameter
    relative_boundary = 0
    if file.solvation_structure != '' and file.crystal_x < 0:
        relative_boundary = 1

    template_dict['filebase'] = file.stripDotPDB(file.filename)
    template_dict['relative_boundary'] = relative_boundary
    template_dict['dim_x'] = str(file.crystal_x) 
    template_dict['dim_z'] = str(file.crystal_z) 
    template_dict['solvation_structure'] = file.solvation_structure

    # if we solvate implicitly or in vacuum, we don't wantto do PBC/PME
    # 20 beyond the structure.
    dopbc = 1
    if file.solvation_structure == '' or solvate_implicitly:  
        dopbc = 0
    elif file.solvation_structure == 'sphere':
        # we can't build a crystal for sphere so we just leave as is and do
        # vacuum boundary conditions
        dopbc = 0

    # we need to set up images, but only if we're doing PBC/PME
    template_dict['dopbc'] = dopbc 
    template_dict['file_location'] = file.location
    if dopbc:
        if file.solvation_structure == '' or solvate_implicitly:
            pass  
        else:
            # we should have an image file from solvation to use...
            try:
                os.stat(file.location + "new_" + file.stripDotPDB(file.filename) + ".xtl")
            except:
                # need to throw some sort of error ... for now just toss cookies 
                return HttpResponse("Oops ... transfer file not found.")

    template_dict['shake'] = request.POST.has_key('apply_shake')
    if request.POST.has_key('apply_shake'):
 	template_dict['which_shake'] = postdata['which_shake']
        if postdata['which_shake'] == 'define_shake':
            template_dict['shake_line'] = postdata['shake_line']
            if postdata['shake_line'] != '':
                file.checkForMaliciousCode(postdata['shake_line'],postdata)

    #md with either be useheat or useequi
    md = postdata['md'] 
    template_dict['md'] = md

    if request.POST.has_key('gen_avgstruct'):
        template_dict['genavg_struct'] = True
    else:    
        template_dict['genavg_struct'] = False

    if md == 'useheat':
        mdp.type = 'heat'
        nstep = postdata['nstepht']
        
        template_dict['nstep'] = nstep
        #If the user does replica exchange, the steps should be divided by 100 to work with the rex scripts
        if(useReplicaExchange):
            nstep = float(nstep)/100
        
        mdp.nstep = str(nstep)
        ndcd = int(int(nstep)/10) # the dynamics trajectories we make have 10 frames
        template_dict['ndcd'] = str(ndcd)
        firstt = postdata['firstt']
        mdp.firstt = str(firstt)
        finalt = postdata['finalt']
        mdp.finalt = str(finalt)
        template_dict['firstt'] = firstt
        template_dict['finalt'] = finalt
        teminc = postdata['teminc']
        mdp.teminc = str(teminc)
        tbath = postdata['tbath']
        mdp.tbath = str(tbath)
        ihtfrq = postdata['ihtfrq']
        mdp.ihtfrq = str(ihtfrq)
        template_dict['ihtfrq'] = ihtfrq
        template_dict['teminc'] = teminc
        template_dict['tbath'] = tbath

        mdp.temp = str(tbath)
        mdp.save()

    else:
        # I hate Django...
        mdp.firstt = '0.0'
        mdp.finalt = '0.0'
        mdp.teminc = '0.0'
        mdp.ihtfrq = '0.0'
        # Done hating Django...
        mdp.type = 'equi'
        nstep = postdata['nstepeq']
       
       #If the user does replica exchange, the steps should be divided by 100 to work with the rex scripts
        if(useReplicaExchange):
            nstep = float(nstep)/100
       
        mdp.nstep = str(nstep)
        template_dict['nstep'] = nstep
        ndcd = int(int(nstep)/10) # the dynamics trajectories we make have 10 frames
        template_dict['ndcd'] = str(ndcd)
        temp = postdata['temp']
        mdp.temp = str(temp)
        mdp.tbath = mdp.temp
        template_dict['temp'] = temp
        mdp.save()

    template_dict['output_name'] = "new_"  + file.stripDotPDB(file.filename) + "-md"
    mdp.save()
    user_id = file.owner.id
    os.chdir(file.location)
    t = get_template('%s/mytemplates/input_scripts/applymd_template.inp' % charmming_config.charmming_root)

    #In case of replica exchange, two output scripts need to be written, one with the start keyword
    #thep second as a restart script
    if(useReplicaExchange == True):

        #Create and fill rexParams object
        rexmdl = rexModule.createRexObjectFromPostData(request,mdp)
        #Create configuration file to send to replica exchange scripts
        rexModule.createConfigFileFromPostData(request,file,mdp)
        #Create the rex wrapper c-shell script to call all of the rex scripts
        rexModule.createRexWrapperCShellScript(file,rexmdl)
        #Clean old files associated with run
        rexModule.removeOldRexRunFiles(file)
        
        #write the non-restart input script
        template_dict['restart'] = 0
        charmm_inp = output.tidyInp(t.render(Context(template_dict)))
        md_filename = "rexstart.inp"
        inp_out = open(file.location + md_filename,'w')
        inp_out.write(charmm_inp)
        inp_out.close()  
       
        #write the restart input script
        template_dict['restart'] = 1
        charmm_inp = output.tidyInp(t.render(Context(template_dict)))
        md_filename = "rex.inp"
        inp_out = open(file.location + md_filename,'w')
        inp_out.write(charmm_inp)
        inp_out.close() 

    #Write configuration file for using the cover script to do replica exchange
    else:
        #else, write the normal start file and continue on
        template_dict['restart'] = 0
        md_filename = "charmm-" + file.stripDotPDB(file.filename) + "-md.inp"
        charmm_inp = output.tidyInp(t.render(Context(template_dict)))
        inp_out = open(file.location + md_filename,'w')
        inp_out.write(charmm_inp)
        inp_out.close()  

    #change the status of the file regarding minimization 
    if file.lesson_type:
        # we need final_pdb_name
        lessonaux.doLessonAct(file,"onMDSubmit",postdata,min_pdb)

    #file.md_status = "<font color=yellow>Processing</font>"
    #If using replica exchange, a custom script must be used to wrap it 
    #having "customscript" in the script list will invoke the backend job submission to use a different script
    #rather than a charmm executable
    exedict = {}
    nprocdict = {}
    if(useReplicaExchange):
        scriptlist.append("customscript")
        exedict["customscript"] = file.location + file.stripDotPDB(file.filename) + "-rexwrapper.csh"
        nprocdict["customscript"] = rexmdl.nbath
    else:
        scriptlist.append(file.location + md_filename)

    if make_movie:
       print "applymd ... makeJmolMovie"
       return makeJmolMovie(file,postdata,min_pdb,scriptlist,'md')
    else:
        print "applymd ... direct submit"
        if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
            return generateHTMLScriptEdit(charmm_inp,scriptlist,'md')
        else:
            si = schedInterface()
            newJobID = si.submitJob(user_id,file.location,scriptlist,exedict,nprocdict)
            file.md_jobID = newJobID
            sstring = si.checkStatus(newJobID)

            mdp.statusHTML = statsDisplay(sstring,newJobID)
            mdp.md_movie_status = None

    file.save()
    return "Done."

#pre: Requires a file object, and the name of the psf/crd file to read in
#Jmol requires the PDB to be outputted a certain a certain way for a movie to be displayed
#it does not take DCD files
def makeJmolMovie(file,postdata,psf_filename,scriptlist,type):
    file.md_movie_status = ""
    charmm_inp = """* Movie making
*

open read unit 2 card name """ + file.stripDotPDB(psf_filename) + """.psf
read psf card unit 2
close unit 2

open read unit 10 file name new_""" + file.stripDotPDB(file.filename) + """-""" + type + """.dcd
traj firstu 10 skip 10

set i 1
label loop
 traj read
 open unit 1 card write name new_""" + file.stripDotPDB(file.filename) + """-movie@i.pdb
 write coor pdb unit 1
 **Cords
 *
 incr i by 1
 if i lt 11 goto loop

stop"""
    #Run the script through charmm, this is not done under job queue at the moment
    #because the PDBs must be appended together after the above script has been run.
    #Once the DAG and query stuff has been implemented, this is how the following code should
    #be changed to
    #1. Create a new field in the PDBFile object called movie_status
    #2. In the status method, check to see if movie_status is done
    #3. If it is done, call the method combinePDBsForMovie(...): right below the following code.
    if(type == 'md'):
        movie_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-mdmovie.inp'
    elif(type == 'ld'):
        movie_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-ldmovie.inp'
    elif(type == 'sgld'):
        movie_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-sgldmovie.inp'
    movie_handle = open(file.location + movie_filename,'w')
    movie_handle.write(charmm_inp)
    movie_handle.close()
    user_id = file.owner.id
    scriptlist.append(file.location + movie_filename)

    si = schedInterface()
    if(type == 'md'):
        if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
            file.save()
            return generateHTMLScriptEdit(charmm_inp,scriptlist,'md')
	else:
            newJobID = si.submitJob(user_id,file.location,scriptlist,{},{})
            file.md_jobID = newJobID
            sstring = si.checkStatus(newJobID)

            mdp = mdParams.objects.filter(pdb=file,selected='y')[0]
            mdp.statusHTML = statsDisplay(sstring,newJobID)
            mdp.md_movie_status = None
            mdp.save()
    elif(type == 'ld'):
        if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
            file.save()
            return generateHTMLScriptEdit(charmm_inp,scriptlist,'ld')
	else:
            newJobID = si.submitJob(user_id,file.location,scriptlist,{},{})
            file.ld_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            file.ld_movie_status = None
            
            ldp = ldParams.objects.filter(pdb=file,selected='y')[0]
            ldp.statusHTML = statsDisplay(sstring,newJobID)
            ldp.ld_movie_status = None
            ldp.save()
    elif(type == 'sgld'):
        if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
            file.save()
            return generateHTMLScriptEdit(charmm_inp,scriptlist,'sgld')
	else:
            newJobID = si.submitJob(user_id,file.location,scriptlist,{},{})
            file.sgld_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            file.sgld_movie_status = None
            
            ldp = sgldParams.objects.filter(pdb=file,selected='y')[0]
            ldp.statusHTML = statsDisplay(sstring,newJobID)
            ldp.save()
    file.save()
#    combinePDBsForMovie(file)
 
#pre: Requires a PDBFile object
#Combines all the smaller PDBs make in the above method into one large PDB that
#jmol understands
#type is md,ld,or sgld
def combinePDBsForMovie(file,type):
    ter = re.compile('TER')
    remark = re.compile('REMARK')
    #movie_handle will be the final movie created
    #mini movie is the PDBs which each time step sepearated into a new PDB
    if(type == 'md'): 
        movie_handle = open(file.location + 'new_' + file.stripDotPDB(file.filename) + "-md-mainmovie.pdb",'a')
    elif(type == 'ld'): 
        movie_handle = open(file.location + 'new_' + file.stripDotPDB(file.filename) + "-ld-mainmovie.pdb",'a')
    elif(type == 'sgld'): 
        movie_handle = open(file.location + 'new_' + file.stripDotPDB(file.filename) + "-sgld-mainmovie.pdb",'a')
    for i in range(10):
        i = i+1
	minimovie_handle = open(file.location +  "new_" + file.stripDotPDB(file.filename) + "-movie" + `i` + ".pdb",'r')
	movie_handle.write('MODEL ' + `i` + "\n")
	for line in minimovie_handle:
	    if(not remark.search(line) and not ter.search(line)):
	        movie_handle.write(line)
	movie_handle.write('ENDMDL\n')
	minimovie_handle.close()
	#os.remove(file.location +  "new_" + file.stripDotPDB(file.filename) + "-movie" + `i` + ".pdb")
	if(type == 'md'):
            mdp = mdParams.objects.filter(pdb=file,selected='y')[0]
	    mdp.md_movie_status = "<font color=33CC00>Done</font>"
            mdp.save()
	elif(type == 'ld'):
            ldp = ldParams.objects.filter(pdb=file,selected='y')[0]
	    ldp.ld_movie_status = "<font color=33CC00>Done</font>"
            ldp.save()
	elif(type == 'sgld'):
            sgldp = sgldParams.objects.filter(pdb=file,selected='y')[0]
	    sgldp.sgld_movie_status = "<font color=33CC00>Done</font>"
            sgldp.save()
	file.save()
    return "Done."
	    
