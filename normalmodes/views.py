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
from structure.models import Structure 
from structure.qmmm import makeQChem_tpl, writeQMheader
from normalmodes.aux import getNormalModeMovieNum
from normalmodes.models import nmodeTask
from account.views import checkPermissions
from structure.aux import checkNterPatch
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from structure.models import Structure, WorkingStructure, WorkingFile, Task
import input, output
import re, copy, os, shutil
import lessonaux, charmming_config
import cPickle

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
        return HttpResponse("Please submit a structure first.")
    try:
        ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
       return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    # check to see if we have prior normal modes to display to the user
    nmode_lines = []
    nmode_fname = ws.structure.location + '/' + ws.identifier + '-nmodes.out'
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
        nt.save()

        if ws.isBuilt != 't':
            isBuilt = False
            pTask = ws.build(nt)
            pTaskID = pTask.id
        else:
            isBuilt = True
            pTaskID = int(request.POST['ptask'])

        return applynma_tpl(request,ws,pTaskID,nt)
    else:
        tasks = Task.objects.filter(workstruct=ws,status='C',active='y').exclude(action='energy')

        lesson_ok, dd_ok = checkPermissions(request)
        return render_to_response('html/normalmodesform.html', {'ws_identifier': ws.identifier,'tasks': tasks, 'nmode_lines': nmode_lines, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})


def applynma_tpl(request,workstruct,pTaskID,nmTask):
    """
    Generates the Normal modes script and runs it
    """

    # try to decide how many atoms this structure has
    pfp = open(workstruct.structure.pickle, 'r')
    pdb = cPickle.load(pfp)
    pfp.close()

    # get the parent task
    pTask = Task.objects.get(id=pTaskID)

    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['input_file'] = workstruct.identifier + '-' + pTask.action
    template_dict['topology_list'], template_dict['parameter_list'], junk = workstruct.getTopparList()
    template_dict['nma'] = request.POST['nma']
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
                return output.returnSubmission('Normal modes', error='You may only do NMA on structures with less than %d atoms' % charmming_config.max_nma_atoms)
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
        input.checkForMaliciousCode(request.POST['qmsele'],request.POST)
        qmparams = makeQchem_val(request.POST,request.POST['qmsele'])
        qmparams['jobtype'] = 'freq'
        template_dict = makeQChem_tpl(template_dict,qmparams,nmTask.workstruct)

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
        #template_dict['headstr'] = headstr
    else:
        nmTask.nma_movie_req = False

    t = get_template('%s/mytemplates/input_scripts/applynma_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))    

    nma_filename = workstruct.identifier + "-nmodes.inp"
    nmTask.scripts += ',%s' % nma_filename
    nmTask.save()

    inp_out = open(workstruct.structure.location + '/' + nma_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()

    #change the status of the file regarding minimization 
    if request.POST.has_key("gen_trj") and request.POST.has_key("num_trjs"):
        return makeNmaMovie_tpl(workstruct,request.POST,int(request.POST['num_trjs']),template_dict['nma'],nmTask)
    else:
        nmTask.start()
        
    workstruct.save()
    return output.returnSubmission('Normal modes')


#makes NMA Movie
def makeNmaMovie_tpl(workstruct,postdata,num_trjs,typeoption,nmm):
    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'], template_dict['parameter_list'], junk = workstruct.getTopparList()
    template_dict['filebase'] = workstruct.identifier
    template_dict['typeoption'] = typeoption
    template_dict['num_trjs'] = str(num_trjs)

    #Run the script through charmm, this is not done under job queue at the moment
    #because the PDBs must be appended together after the above script has been run.
    #Once the DAG and query stuff has been implemented, this is how the following code should
    #be changed to
    #1. Create a new field in the Structure object called movie_status
    #2. In the status method, check to see if movie_status is done
    #3. If it is done, call the method combinePDBsForMovie(...): right below the following code.

    t = get_template('%s/mytemplates/input_scripts/makeNmaMovie_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    
    movie_filename = workstruct.identifier + '-nmamovie.inp'
    movie_handle = open(workstruct.structure.location + movie_filename,'w')
    movie_handle.write(charmm_inp)
    movie_handle.close()
    nmm.scripts += ',' + movie_filename

    ## There are no NMA lessons at the moment ... when there are, this needs to be redone correctly
    ##if file.lesson_type:
    ##    lessonaux.doLessonAct(file,"onNMASubmit",postdata,None)
    #nmm.nma_movie_status = "<font color=33CC00>Done</font>"
    nmm.start()
    nmm.save()
    return output.returnSubmission('Normal modes')

#pre: Requires a PDBFile object
#Combines all the smaller PDBs make in the above method into one large PDB that
#jmol understands
#type is nma
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

    nmm.nma_movie_status = "<font color=33CC00>Done</font>"
    nmm.make_nma_movie = False
    nmm.save()
    file.save()
