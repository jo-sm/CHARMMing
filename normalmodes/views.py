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
from structure.models import Structure, noNscaleFound
from django.contrib import messages
from structure.qmmm import makeQChem_tpl, makeQchem_val, writeQMheader
from normalmodes.aux import getNormalModeMovieNum
from normalmodes.models import nmodeTask
from account.views import checkPermissions
from structure.aux import checkNterPatch
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from structure.models import Structure, WorkingStructure, WorkingFile, Task, CGWorkingStructure
from atomselection_aux import getAtomSelections, saveAtomSelections
import selection.models, structure.mscale
from django.contrib.messages import get_messages
import input, output
import re, copy, os, shutil
import lessonaux, charmming_config
import cPickle
from lesson_config import *

def normalmodesformdisplay(request):
    """
    Display the form for normal modes calculation
    """

    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)
    try:
        struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please submit a structure first.")
        return HttpResponseRedirect("/charmming/fileupload/")
#        return HttpResponse("Please submit a structure first.")
    try:
        ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "Please build your structure before performing any calculations.")
        return HttpResponseRedirect("/charmming/buildstruct/")

    cg_model = False
    try:
        cgws = CGWorkingStructure.objects.get(workingstructure_ptr=ws.id)
        cg_model = True
#        messages.error(request, "Normal mode analysis is not supported for coarse-grain models. If you wish to perform normal mode analysis on this structure, please build a working structure with the CHARMM all-atom model.")
#        return HttpResponseRedirect("/charmming/buildstruct/")
    except:
        pass
#       return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")
    tdict = {}
    # check to see if we have prior normal modes to display to the user
    tdict['cg_model'] = cg_model
    nmode_lines = []
    nmode_fname = ws.structure.location + '/' + ws.identifier + '-nmode.out'
    try:
        os.stat(nmode_fname)
    except:
        pass
    else:
        nmodefp = open(nmode_fname, 'r')
        getnxt = False
        ntoget = 10
        for line in nmodefp:
            line = line.strip()
            if line.startswith('VIBRATION MODE'):
                freq = float(line.split()[4])
                if freq > 0.00001:
                    nmode_lines.append(line)
                    getnxt = True
            elif getnxt:
                nmode_lines.append(line)
                getnxt = False
                ntoget -= 1

            if ntoget == 0:
                break
        nmodefp.close()

    if request.POST.has_key('num_normalmodes'):
        # if there is a previous solvation structure, deactivate it
        try:
            oldtsk = nmodeTask.objects.filter(workstruct=ws,active='y')[0]
            oldtsk.active = 'n'
            oldtsk.save()
        except:
            pass

        nt = nmodeTask()
        nt.setup(ws)
        nt.active = 'y'
        nt.action = 'nmode'
        nt.modifies_coordinates = False
        nt.save()

        if ws.isBuilt != 't':
            isBuilt = False
            try:
                pTask = ws.build(nt)
            except noNscaleFound, e:
                return output.returnSubmission('Minimization', error='The nScale parameterization process has not yet completed. It may take 1-2 hours.')
            pTaskID = pTask.id
        else:
            isBuilt = True
            pTaskID = int(request.POST['ptask'])
            pTask = Task.objects.get(id=pTaskID)
            ptask_path = "%s/%s-%s.psf"%(struct.location,ws.identifier,pTask.action)
            try:
                os.stat(ptask_path)
            except: #probably a qchem thing, so copy the build PSF
                shutil.copyfile("%s/%s-build.psf"%(struct.location,ws.identifier),ptask_path) #TODO: warn user?
                shutil.copyfile("%s/%s-build.crd"%(struct.location,ws.identifier),ptask_path.replace(".psf",".crd")) #TODO: warn user?
                #crd doesn't change for qchem stuff...

        if request.POST.has_key('useqmmm'):
            nt.useqmmm = 'y'
            nt.save()
            saveAtomSelections(request, ws, pTask)

        return applynma_tpl(request,ws,pTaskID,nt)
    else:
        sl = []
        sl.extend(structure.models.Segment.objects.filter(structure=struct,is_working='n'))
        sl = sorted(map(lambda x: x.name,sl))
        #Somehow this works differently than ener or minimization, and displays "Segment object" instead of the segment name...I Don't know what's going on here.
        #Actually, I don't understand how seg_list even works in the others, since the loop goes:
            #for item in seg_list:
            #    {{item}}
            #Which would just print Segment object...how it works outside of here is a bit of a mystery./
        tdict['seg_list'] = sl
        tasks = Task.objects.filter(workstruct=ws,status='C',active='y',modifies_coordinates=True)
        tdict = getAtomSelections(tdict,ws)
        lesson_ok, dd_ok = checkPermissions(request)
        tdict['messages'] = get_messages(request)
        tdict['ws_identifier'] = ws.identifier
        tdict['tasks'] = tasks
        tdict['nmode_lines'] = nmode_lines
        return render_to_response('html/normalmodesform.html', tdict)


def applynma_tpl(request,workstruct,pTaskID,nmTask):
    """
    Generates the Normal modes script and runs it
    """

    # try to decide how many atoms this structure has
    pfile = workstruct.localpickle if workstruct.localpickle else workstruct.structure.pickle
    pfp = open(pfile, 'r')
    pdb = cPickle.load(pfp)
    pfp.close()

    # get the parent task
    pTask = Task.objects.get(id=pTaskID)

    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['input_file'] = workstruct.identifier + '-' + pTask.action
    template_dict['topology_list'], template_dict['parameter_list'] = workstruct.getTopparList()
    template_dict['nma'] = request.POST['nma']
    if workstruct.topparStream:
        template_dict['tpstream'] = workstruct.topparStream.split()
    else:
        template_dict['tpstream'] = []


    # save model
    if template_dict['nma'] == 'useenm':
        nmTask.type = 2
        nmTask.rcut = float(request.POST['rcut'])
        nmTask.kshort = float(request.POST['kshort'])
        nmTask.klong = float(request.POST['klong'])
        template_dict['rcut'] = nmTask.rcut
        template_dict['kshort'] = nmTask.kshort
        template_dict['klong'] = nmTask.klong
    else:
        nmTask.type = 1
    nmTask.nmodes = int(request.POST['num_normalmodes'])

    if 'append_' + workstruct.identifier in pdb.keys():
        if template_dict['nma'] == 'usevibran':
            if len(pdb['append_' + workstruct.identifier]) > charmming_config.max_nma_atoms:
                messages.error(request, "You may only do NMA on structures with less than %d atoms" % charmming_config.max_nma_atoms)
                return HttpResponseRedirect("/charmming/normalmodes/")
#                return output.returnSubmission('Normal modes', error='You may only do NMA on structures with less than %d atoms' % charmming_config.max_nma_atoms)
    else:
        # this must be the case where the appending has not happened yet
        # for now let's allow this but we need to figure out some way
        # of counting the atoms of a pre-appended struct
        pass


    template_dict['identifier'] = workstruct.identifier
    template_dict['numnormalmodes'] = nmTask.nmodes
    template_dict['useqmmm'] = request.POST.has_key("useqmmm")

    # take care of any QM/MM that the user might want.
    if template_dict['nma'] == 'usevibran' and request.POST.has_key("useqmmm"):
        #input.checkForMaliciousCode(request.POST['qmsele'],request.POST) #It's being passed escaped...
        #We pre-validate. This is obsolete.
#        if "script" in request.POST['qmsele'] or "&lt;" in request.post['qmsele']: #Script injection huh?
#Obsolete due to ONIOM and other multi-scale improvements.
#            messages.error(request, "Bad QM selection.")
#            return HttpResponseRedirect('/charmming/about/')
        if request.POST['model_selected'] not in ["oniom","qmmm"]:
            messages.error(request, "Invalid Multi-Scale modeling type.")
            return HttpResponseRedirect("/charmming/about/")
        else:
            modelType = request.POST['model_selected']
        template_dict['modelType'] = modelType
        nmTask.modelType = modelType
        if modelType == "qmmm":
            atomselection = selection.models.AtomSelection.objects.get(workstruct=workstruct)
            qmparams = makeQchem_val("qmmm",atomselection)
            qmparams['jobtype'] = 'freq'
            template_dict = makeQChem_tpl(template_dict,qmparams,False,nmTask)
        elif modelType == "oniom":
            template_dict = structure.mscale.make_mscale(template_dict, request, modelType, nmTask)
        nmTask.save()

    # If need be, print out trajectories for the modes the user requested.
    template_dict['gen_trj'] = request.POST.has_key("gen_trj") 
    template_dict['num_trjs'] = request.POST.has_key("num_trjs") 
    if request.POST.has_key("gen_trj") and request.POST.has_key("num_trjs"):
        # movies are liable to be broken for the time being
        nmTask.nma_movie_req = True
        try:
            ntrjs = int(request.POST['num_trjs'])
        except:
            ntrjs = 5
        if ntrjs > 20:
            ntrjs = 20
        if ntrjs > nmTask.nmodes:
            ntrjs = nmTask.nmodes
        #if request.POST.has_key("useqmmm"):
        #    headstr = writeQMheader("", "SELE " + qmsel + " END")
        #else:
        #    headstr = "* Not using QM/MM\n"
        template_dict['ntrjs'] = str(ntrjs)
        nmTask.num_trjs = ntrjs
        #template_dict['headstr'] = headstr
    else:
        nmTask.nma_movie_req = False

    t = get_template('%s/mytemplates/input_scripts/applynma_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))    

    nma_filename = workstruct.structure.location + '/' + workstruct.identifier + "-nmode.inp"
    if nmTask.useqmmm == 'y' and modelType == "oniom":
        nmTask.add_script('charmm-mscale',nma_filename,charmming_config.default_mscale_nprocs)
    else:
        nmTask.add_script('charmm',nma_filename,charmming_config.default_charmm_nprocs)
#    nmTask.scripts += ',%s' % nma_filename
    nmTask.save()

    inp_out = open(nma_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()

    #change the status of the file regarding minimization 
    if request.POST.has_key("gen_trj") and request.POST.has_key("num_trjs"):
        nmTask.make_nma_movie = True
        nmTask.save()
        return makeNmaMovie_tpl(workstruct,request.POST,int(request.POST['num_trjs']),template_dict['nma'],nmTask,pTask) #I have no idea if pTask is stored in nmTask, so let's be safe
    else:
        nmTask.start()

    ## There are no NMA lessons at the moment ... when there are, this needs to be redone correctly
    #Lesson_maker implements NMA lesson construction. It has been redone.
    if workstruct.structure.lesson_type:
        lessonaux.doLessonAct(workstruct.structure,"onNMASubmit",nmTask,None)

    workstruct.save()
    return output.returnSubmission('Normal modes')


#makes NMA Movie
def makeNmaMovie_tpl(workstruct,postdata,num_trjs,typeoption,nmm,ptask):
    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'], template_dict['parameter_list'] = workstruct.getTopparList()
    template_dict['orig_coords'] = workstruct.identifier + "-" + ptask.action #Only inp/out files don't match the task.action scheme. Careful with dsfdocking though.
    template_dict['filebase'] = workstruct.identifier
    template_dict['typeoption'] = typeoption
    template_dict['num_trjs'] = str(num_trjs)
    if workstruct.topparStream:
        template_dict['tpstream'] = workstruct.topparStream.split()
    else:
        template_dict['tpstream'] = []


    #Run the script through charmm, this is not done under job queue at the moment
    #because the PDBs must be appended together after the above script has been run.
    #Once the DAG and query stuff has been implemented, this is how the following code should
    #be changed to
    #1. Create a new field in the Structure object called movie_status
    #2. In the status method, check to see if movie_status is done
    #3. If it is done, call the method combinePDBsForMovie(...): right below the following code.

    t = get_template('%s/mytemplates/input_scripts/makeNmaMovie_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    
    movie_filename = workstruct.structure.location + "/" + workstruct.identifier + '-nmamovie.inp'
    movie_handle = open(movie_filename,'w')
    movie_handle.write(charmm_inp)
    movie_handle.close()
    if nmm.useqmmm == 'y' and nmm.modelType == "oniom":
        nmm.add_script("charmm-mscale",movie_filename,charmming_config.default_mscale_nprocs)
    else:
        nmm.add_script('charmm',movie_filename,charmming_config.default_charmm_nprocs)
#    nmm.scripts += ',' + movie_filename

    #nmm.nma_movie_status = "<span style='color:33CC00'>Done</span>"
    nmm.start()
    nmm.save()
    return output.returnSubmission('Normal modes')

#pre: Requires a PDBFile object
#Combines all the smaller PDBs make in the above method into one large PDB that
#jmol understands
#type is nma
#Ignore this for now...there's a better one attached to the nmTask
def combineNmaPDBsForMovie(file):
    nmm = nmodeTask.objects.filter(pdb=file, selected = 'y')[0]
    ter = re.compile('TER')
    remark = re.compile('REMARK')
    #movie_handle will be the final movie created
    #mini movie is the PDBs which each time step sepearated into a new PDB
    trj_num = 0
    continueit = True
    while(continueit):
        try:
            pdbno = trj_num + 1
            os.stat(file.location + 'new_' + file.stripDotPDB(file.filename) + '-mtraj_' + str(pdbno) + '.trj')
            trj_num += 1 
        except:
            continueit = False
    currtrjnum = 1
    while currtrjnum < (trj_num+1):
        movie_handle = open(file.location + 'new_' + file.stripDotPDB(file.filename) + "-nma-mainmovie-" + str(currtrjnum) + ".pdb",'a')
        for i in range(12):
            i = i+1
            minimovie_handle = open(file.location +  "new_" + file.stripDotPDB(file.filename) + "-nmapremovie" + `i` + "-" + str(currtrjnum) + ".pdb",'r')
            movie_handle.write('MODEL ' + `i` + "\n")
            for line in minimovie_handle:
                if(not remark.search(line) and not ter.search(line)):
                    movie_handle.write(line)
            movie_handle.write('ENDMDL\n')
            minimovie_handle.close()
            os.remove(file.location +  "new_" + file.stripDotPDB(file.filename) + "-nmapremovie" + `i` + "-" + str(currtrjnum) + ".pdb")
        currtrjnum += 1

    nmm.nma_movie_status = "<span class='done'>Done</span>" #where do we even use this...?
    nmm.make_nma_movie = False
    nmm.save()
    file.save()
