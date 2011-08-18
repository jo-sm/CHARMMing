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
from structure.models import Structure, WorkingStructure, WorkingFile, Task
from account.views import isUserTrustworthy
from structure.aux import checkNterPatch
from dynamics.models import mdTask, ldTask, sgldTask
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from solvation.models import solvationTask
import django.forms
import copy, shutil, os, re
import lessonaux, input, output, charmming_config

#processes form data for ld simulations
def lddisplay(request):
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

    if request.POST.has_key('nstep'):
        # form has been submitted, check and see if there is an old ldTask, and
        # if so, deactivate it.
        ldt = None

        if request.POST.has_key('usesgld'):
            try:
                oldtsk = sgldTask.objects.filter(workstruct=workstruct,active='y')[0]
                oldtsk.active = 'n'
                oldtsk.save()
            except:
                pass
            ldt = sgldTask()
            ldt.setup(ws)
            ldt.action = 'sgld'
            ldt.save()
        else:
            try:
                oldtsk = ldTask.objects.filter(workstruct=workstruct,active='y')[0]
                oldtsk.active = 'n'
                oldtsk.save()
            except:
                pass
            ldt = ldTask()
            ldt.setup(ws)
            ldt.action = 'ld'
            ldt.save()

        ldt.nstep = int(request.POST['nstep'])
        ldt.active = 'y'
        ldt.save()

        if ws.isBuilt != 't':
            isBuilt = False
            pTask = ws.build(ldt)
            pTaskID = pTask.id
        else:
            isBuilt = True
            pTaskID = int(request.POST['ptask'])

        return applyld_tpl(request,ldt,pTaskID)
  
    else:
        # get all workingFiles associated with this struct
        tasks = Task.objects.filter(workstruct=ws,status='C',active='y')
        return render_to_response('html/ldform.html', {'ws_identifier': ws.identifier,'tasks': tasks})

#processes form data for md simulations
def mddisplay(request):
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


    if request.POST.has_key('mdtype'):
        try:
            oldtsk = mdTask.objects.filter(workstruct=workstruct,active='y')[0]
            oldtsk.active = 'n'
            oldtsk.save()
        except:
             pass

        mdt = mdTask()
        mdt.setup(ws)
        mdt.active = 'y'
        mdt.action = 'md'
        mdt.save()

        if ws.isBuilt != 't':
            isBuilt = False
            pTask = ws.build(scriptlist)
            pTaskID = pTask.id
        else:
            isBuilt = True
            pTaskID = int(request.POST['ptask'])

        return applymd_tpl(request,mdt,pTaskID)

    else:
        # decide whether or not this run is restartable (check for non-empty restart file)
        try:
            os.stat(ws.structure.location + "/md-" + ws.identifier + ".res")
            canrestart = True
        except:
            canrestart = False

        # get all workingFiles associated with this struct
        tasks = Task.objects.filter(workstruct=ws,status='C',active='y')
        return render_to_response('html/mdform.html', {'ws_identifier': ws.identifier,'tasks': tasks, 'canrestart': canrestart})

def applyld_tpl(request,ldt,pTaskID):
    postdata = request.POST

    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'] = ldt.workstruct.getTopologyList()
    template_dict['parameter_list'] = ldt.workstruct.getParameterList()
    template_dict['fbeta'] = postdata['fbeta']
    template_dict['nstep'] = postdata['nstep']
    template_dict['usesgld'] = postdata.has_key('usesgld')
    template_dict['identifier'] = ldt.workstruct.identifier

    ldt.fbeta = template_dict['fbeta']
    ldt.nstep = template_dict['nstep']

    pTask = Task.objects.get(id=pTaskID)
    ldt.parent = pTask
    template_dict['input_file'] = ldt.workstruct.identifier + '-' + pTask.action
    template_dict['useqmmm'] = ''
    template_dict['qmmmsel'] = ''
    ldt.save()

    # deal with whether or not we want to go to Hollywood (i.e. make a movie)
    # NB: this is not working at present; flesh this code in!!!
    try:
        make_movie = postdata['make_movie']
    except:
        make_movie = None


    #If the user wants to solvate implicitly the scpism line is needed
    #84 will be the scpism number in this program
    template_dict['solvate_implicitly'] = False
    try:
        template_dict['solvate_implicitly'] = postdata['solvate_implicitly']
        if(postdata['solvate_implicitly']):
	    ldt.scpism = True
    except:
        pass

    if ldt.action == 'sgld':
        tsgavg = '0.0'
        tempsg = '0.0'

        tsgavg = postdata['tsgavg']
	ldt.tsgavg = str(tsgavg)
	tempsg = postdata['tempsg']
	ldt.tempsg = str(tempsg)

        template_dict['tsgavg'] = tsgavg
        template_dict['tempsg'] = tempsg

    # determine whether periodic boundary conditions should be applied, and if
    # so pass the necessary crystal and image parameters
    dopbc = False
    if request.POST.has_key('usepbc'):
        if request.POST.has_key('solvate_inplicitly'):
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
            sp = solvationTask.objects.filter(structure=workstruct,active='y')[0]
        except:
            return HttpResponse("Err ... couldn't find solvation parameters")
        template_dict['xtl_x'] = sp.xtl_x
        template_dict['xtl_y'] = sp.xtl_y
        template_dict['xtl_z'] = sp.xtl_z
        template_dict['xtl_angles'] = "%10.6f %10.6f %10.6f" % (sp.angles[0],sp.angles[1],sp.angles[2])
        template_dict['xtl_ucell'] = sp.solvation_structure
        template_dict['ewalddim'] = sp.calcEwaldDim()

    template_dict['dopbc'] = dopbc

    template_dict['output_name'] = ldt.workstruct.identifier + '-' + ldt.action
    user_id = ldt.workstruct.structure.owner.id
    ld_filename = ldt.workstruct.structure.location + '/' + ldt.workstruct.identifier + "-" + ldt.action + ".inp"
    t = get_template('%s/mytemplates/input_scripts/applyld_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))

    inp_out = open(ld_filename, 'w')
    inp_out.write(charmm_inp)
    inp_out.close()  
    ldt.scripts += ',%s' % ld_filename

    # change the status of the file regarding LD
    ldt.save()


    # lessons are borked at the moment
    #if file.lesson_type:
    #    if usesgld:
    #        lessonaux.doLessonAct(file,"onSGLDSubmit",postdata,None)
    #    else:
    #        lessonaux.doLessonAct(file,"onLDSubmit",postdata,None)

    if make_movie:
        # NB, fix me up...
        if usesgld:
            return makeJmolMovie(file,postdata,min_pdb,scriptlist,'sgld')
        else:
            return makeJmolMovie(file,postdata,min_pdb,scriptlist,'ld')
    else:
        ldt.start()
        ldt.save()

    ldt.workstruct.save()
    return HttpResponse("Done.")

#Generates MD script and runs it
def applymd_tpl(request,mdt,pTaskID):
    postdata = request.POST

    mdt.make_movie = postdata.has_key('make_movie')
    mdt.scpism = postdata.has_key('solvate_implicitly')
    
    # template dictionary passes the needed variables to the templat
    template_dict = {}     
    template_dict['topology_list'] = mdt.workstruct.getTopologyList()
    template_dict['parameter_list'] = mdt.workstruct.getParameterList()

    template_dict['output_name'] = mdt.workstruct.identifier + '-md'

    pTask = Task.objects.get(id=pTaskID)
    template_dict['input_file'] = mdt.workstruct.identifier + '-' + pTask.action
    template_dict['solvate_implicitly'] = mdt.scpism
    mdp.inpStruct = pstruct

    orig_rst = mdt.workstruct.structure.location + "/" + pstruct.basename + "-md.res"
    new_rst  = mdt.workstruct.structure.location + "/" + pstruct.basename + "-md-old.res"
    if postdata.has_key('dorestart'):
        template_dict['restartname'] = new_rst
        template_dict['strtword'] = "restart"
        try:
            shutil.copy(orig_rst,new_rst)
        except:
            pass
    else:
        template_dict['strtword'] = "start"

    # determine whether periodic boundary conditions should be applied, and if
    # so pass the necessary crystal and image parameters
    dopbc = False
    if request.POST.has_key('usepbc'):
        if request.POST.has_key('solvate_inplicitly'):
            return HttpResponse('Invalid options')
        
        # decide if the structure we're dealing with has
        # been solvated.
        solvated = False
        wfc = pTask
        while True:
            if wfc.action == 'solvate':
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
            sp = solvationTask.objects.filter(structure=workstruct,active='y')[0]
        except:
            return HttpResponse("Err ... couldn't find solvation parameters")
        template_dict['xtl_x'] = sp.xtl_x
        template_dict['xtl_y'] = sp.xtl_y
        template_dict['xtl_z'] = sp.xtl_z
        template_dict['xtl_angles'] = "%10.6f %10.6f %10.6f" % (sp.angles[0],sp.angles[1],sp.angles[2])
        template_dict['xtl_ucell'] = sp.solvation_structure
        template_dict['ewalddim'] = sp.calcEwaldDim()

    template_dict['dopbc'] = dopbc 

    # The MD type can have three values: useheat, usenve, or
    # usenvt (todo add NPT)
    template_dict['mdtype'] = postdata['mdtype']
    template_dict['genavg_structure'] = postdata.has_key('gen_avgstruct')

    if template_dict['mdtype'] == 'useheat':
        mdt.ensemble = 'heat'
        mdt.firstt = postdata['firstt']
        mdt.finalt = postdata['finalt']
        template_dict['nstep'] = postdata['nstepht']
        template_dict['firstt'] = postdata['firstt']
        template_dict['finalt'] = postdata['finalt']
        template_dict['tbath'] = postdata['tbath']
        template_dict['teminc'] = postdata['teminc']
        template_dict['ihtfrq'] = postdata['ihtfrq']
    elif template_dict['mdtype'] == 'usenve':
        mdt.ensemble = 'nve'
        template_dict['nstep'] = postdata['nstepnve']
        template_dict['ieqfrq'] = postdata['ieqfrq']
    elif template_dict['mdtype'] == 'usenvt':
        mdt.ensemble = 'nvt'
        mdt.temp = postdata['hoovertemp']
        template_dict['nstep'] = postdata['nstepnvt']
        template_dict['hoovertemp'] = postdata['hoovertemp']
    elif template_dict['mdtype'] == 'usenpt':
        mdt.ensemble = 'npt'
        mdp.temp = postdata['hoovertemp']
        template_dict['nstep'] = postdata['nstepnpt']
        template_dict['hoovertemp'] = postdata['hoovertemp']
        template_dict['pgamma'] = postdata['pgamma']

    if dopbc or mdt.ensemble == 'nvt' or mdt.ensemble == 'npt':
        template_dict['cpt'] = 'cpt'
    else:
        template_dict['cpt'] = ''

    user_id = mdt.workstruct.structure.owner.id
    os.chdir(mdt.workstruct.structure.location)
    t = get_template('%s/mytemplates/input_scripts/applymd_template.inp' % charmming_config.charmming_root)

    # write out the file and let it go...
    md_filename = workstruct.structure.location + "/" + workstruct.identifier + "-md.inp"
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    inp_out = open(md_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()  

    mdt.scripts += ',%s' % md_filename
    mdt.save()

    # lessons are still borked
    #if file.lesson_type:
    #    lessonaux.doLessonAct(file,"onMDSubmit",postdata,min_pdb)

    make_movie = postdata.has_key('make_movie')
    if make_movie:
        return makeJmolMovie(workstruct,postdata,scriptlist,'md')
    else:
        mdt.start()

    mdt.save()
    mdt.workstruct.save()
    return HttpResponse("Done.")

#pre: Requires a file object, and the name of the psf/crd file to read in
#Jmol requires the PDB to be outputted a certain a certain way for a movie to be displayed
#it does not take DCD files
def makeJmolMovie(workstruct,postdata,scriptlist,type):
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
    #1. Create a new field in the Structure object called movie_status
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
	else:
            newJobID = si.submitJob(user_id,file.location,scriptlist,{},{})
            file.md_jobID = newJobID
            sstring = si.checkStatus(newJobID)

            mdp = mdTask.objects.filter(pdb=file,selected='y')[0]
            mdp.statusHTML = statsDisplay(sstring,newJobID)
            mdp.md_movie_status = None
            mdp.save()
    elif(type == 'ld'):
        if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
            file.save()
	else:
            newJobID = si.submitJob(user_id,file.location,scriptlist,{},{})
            file.ld_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            file.ld_movie_status = None
            
            ldp = ldTask.objects.filter(pdb=file,selected='y')[0]
            ldp.statusHTML = statsDisplay(sstring,newJobID)
            ldp.ld_movie_status = None
            ldp.save()
    elif(type == 'sgld'):
        if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
            file.save()
	else:
            newJobID = si.submitJob(user_id,file.location,scriptlist,{},{})
            file.sgld_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            file.sgld_movie_status = None
            
            ldp = sgldTask.objects.filter(pdb=file,selected='y')[0]
            ldp.statusHTML = statsDisplay(sstring,newJobID)
            ldp.save()
    file.save()
#    combinePDBsForMovie(file)
 
#pre: Requires a Structure object
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
            mdp = mdTask.objects.filter(pdb=file,selected='y')[0]
	    mdp.md_movie_status = "<font color=33CC00>Done</font>"
            mdp.save()
	elif(type == 'ld'):
            ldp = ldTask.objects.filter(pdb=file,selected='y')[0]
	    ldp.ld_movie_status = "<font color=33CC00>Done</font>"
            ldp.save()
	elif(type == 'sgld'):
            sgldp = sgldTask.objects.filter(pdb=file,selected='y')[0]
	    sgldp.sgld_movie_status = "<font color=33CC00>Done</font>"
            sgldp.save()
	file.save()
    return "Done."
	    
