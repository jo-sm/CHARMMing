
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
from django.contrib import messages as messages
from django.template import *
from django.template import RequestContext
from account.views import checkPermissions
from lessons.models import LessonProblem
from mutation.models import mutateTask
import output, lesson1, lesson2, lesson3, lesson4, lesson5, lesson6, lessonaux
import structure.models, input
from lesson_config import *
import os, copy, json, mimetypes, string, re
import pychm.io, charmming_config
from pychm.const.bio import aaAlphabet
import cPickle

def mutestructure(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    input.checkRequestData(request)
    postdata = request.POST

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please upload a structure first.") #There should only ever be one message but we'll generalize.
        return HttpResponseRedirect('/charmming/fileupload/') #I have no idea why this would happen but you never know.

    try:
        ws = structure.models.WorkingStructure.objects.get(structure=struct,selected='y')
    except:
        messages.error(request, "Please build a working structure before performing any calculations.")
        return HttpResponseRedirect('/charmming/buildstruct')

    try:
        os.stat(workstruct.structure.location + "/" + postdata['MutFile'] + "-mutation.pdb")
        messages.error(request, "There already exists a mutation under that name. Please delete this mutation file or choose a different name.")
        return HttpResponseRedirect('/charmming/mutation') #This is unlikely to happen. I dont' want to restructure the whole page because of this one error,
    #an error which would only happen if the user REFUSES to use our automatic generator and picks a name that already exists, on their own.
    except: #Basically if the file DOESN'T exist, keep going
        pass

    MutResi = postdata['MutResi']
    MutResName = postdata['MutResName'] #The script doesn't actually care about this one, it's just for user reference
    NewResn = postdata['NewResn']
    MutSegi = postdata['MutSegi']
    MutFile = postdata['MutFile']

    if len(MutFile) > 20:
        messages.error(request, "Working structure identifier is too long. Please make sure your identifier is under 20 characters.")
        return HttpResponseRedirect("/charmming/mutation/")

    modstruct = structure.models.WorkingStructure()
    modstruct.structure = ws.structure
    modstruct.doblncharge = copy.deepcopy(ws.doblncharge) #This may be a primitive but I don't care, let's be safe
    modstruct.isBuilt = copy.deepcopy(ws.isBuilt)
    modstruct.modelName = copy.deepcopy(ws.modelName)
    modstruct.qmRegion = copy.deepcopy(ws.qmRegion) #What is this?
    modstruct.finalTopology = copy.deepcopy(ws.finalTopology)
    modstruct.finalParameter = copy.deepcopy(ws.finalParameter)
    modstruct.topparStream = copy.deepcopy(ws.topparStream)
    modstruct.localpickle = ws.structure.location + "/" + MutFile + "-pickle.dat"
    modstruct.lesson = None #No association - I don't think we should keep track of this unless there's a lesson for mutation?
    modstruct.extraStreams = copy.deepcopy(ws.extraStreams)
    modstruct.identifier = MutFile #Can we do this? I need to do it to make sure the filename percolates up correctly
    modstruct.save()

    # horrific hack alert --
    # Copy all of the working segment identifiers from the parent workingstructure into this one. The reason that
    # this works is because the mutate task does this actually mutation AND builds a new PSF, so the PSFs of the
    # working segments are never again used in a CHARMM script. This is definitely not a clean way of doing things,
    # but it ought to work...
    for seg in ws.segments.all():
        modstruct.segments.add(seg)

    modstruct.save()

    ws.selected = 'n'
    ws.save()
    modstruct.selected = 'y'
    modstruct.save()

    if modstruct.localpickle:
        woof = modstruct.localpickle
    else:
        woof = struct.pickle #No specific pickle file, possibly not mutated before

    mt = mutateTask()
    mt.setup(modstruct)
    mt.lpickle = woof
    mt.MutResi = MutResi
    mt.MutResName = MutResName
    mt.NewResn = NewResn
    mt.MutSegi = MutSegi
    mt.MutFile = MutFile + "-mutation" #THis way all other tasks don't break...
    mt.active = 'y'
    mt.action = 'mutation' #I assume these are inherited from Task...
    mt.parent = None #Important; mutate Tasks have no parent; equivalent to build Tasks
    mt.save()
    #Since there's no parent task we can skip over it...
    return mutate_task_process(request,mt,ws) #PDB is passed so we can work with it...

def mutate_task_process(request,mt,oldws):
    postdata = request.POST
    logfp = open("/tmp/log_mutate_process.txt","w")
    os.chdir(mt.workstruct.structure.location)
    logfp.write("localpickle = " + str(mt.lpickle) + "\n")
    logfp.write("MutResi = " + str(mt.MutResi) + "\n")
    logfp.write("MutSegi = " + str(mt.MutSegi) + "\n")
    logfp.write("MutFile = " + str(mt.MutFile) + "\n")
    template_dict = {}
    template_dict['topology_list'], template_dict['parameter_list'] = mt.workstruct.getTopparList()
    if mt.workstruct.topparStream:
        template_dict['tpstream'] = mt.workstruct.topparStream.split()
    else:
        template_dict['tpstream'] = []
    template_dict['output_name'] = mt.workstruct.identifier + '-mutate'
    template_dict['MutResi'] = mt.MutResi
    template_dict['NewResn'] = mt.NewResn
    template_dict['MutSegi'] = mt.MutSegi
    template_dict['MutFile'] = mt.MutFile
    template_dict['BaseStruct'] = mt.workstruct.structure.name
    logfp.close()
#    pTask = Task.objects.filter(id=pTaskID)[0]
    action = ""
    ws = mt.workstruct
    if ws.isBuilt != 't': #In theory this should work, but it doesn't. Not sure why.
        isBuilt = False #??? It serves no purpose in minimization either...
        pTask = oldws.build(mt)
        pTaskID = pTask.id
    else:
        isBuilt = True
        pTaskID = int(request.POST['ptask']) #THis comes from the previous form...
    pTask = Task.objects.filter(id=pTaskID)[0]
    action = pTask.action
    template_dict['input_file'] = oldws.identifier + "-" + action
    #So now there's a reference to the parent task which makes the rest of charmming not break

    mt.parent = None
    mt.save()

    #THis way we ge tthe correct parent and don't break everything...
    template_dict['restraints'] = ''
    try:
        postdata['apply_restraints']
        template_dict['restraints'] = file.handleRestraints(request)
    except:
        pass

    mt.save()
    t = get_template('%s/mytemplates/input_scripts/mutation_template.inp' %(charmming_config.charmming_root))
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))

    user_id = ws.structure.owner.id
    mutate_filename = ws.structure.location + "/" + ws.identifier + "-mutation.inp"
    inp_out = open(mutate_filename, 'w')
    inp_out.write(charmm_inp)
    inp_out.close()
    mt.scripts += ',%s' % mutate_filename
    mt.start()
    mt.save()
    mt.workstruct.save()
    return output.returnSubmission('Mutation') #Will this break everything?




#Point-mutation page load. Actual mutation happens elsewhere.
def selectstructure(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    logfp = open("/tmp/sortedtest.txt","w")
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0] #Gets the currently-selected structure, which you can switch later
    except:
        messages.error(request, "Please build a structure first.") #There should only ever be one message but we'll generalize.
        return HttpResponseRedirect('/charmming/fileupload/')

    proposedname = struct.name
    existingWorkStructs = structure.models.WorkingStructure.objects.filter(structure=struct)
    if existingWorkStructs:
        usedNameList = [ews.identifier for ews in existingWorkStructs]

        while proposedname in usedNameList:
           m = re.search("([^-\s]+)-mt([0-9]+)$", proposedname)
           if m:
               basename = m.group(1)
               numbah   = int(m.group(2))
               proposedname = basename + "-mt" + str(numbah+1)
           else:
               proposedname += "-mt1"
    else:
        messages.error(request, "Please build a working structure before performing a mutation.", extra_tags="arf")
        return HttpResponseRedirect('/charmming/buildstruct/')
    try: #Just in case...
        ws = structure.models.WorkingStructure.objects.get(structure=struct,selected='y')
    except:
        messages.error(request, "Please build a working structure before performing a mutation.", extra_tags="arf")
        return HttpResponseRedirect('/charmming/buildstruct/')
    if ws.isBuilt != "t":
        messages.error(request, "Please perform a calculation on this structure before performing a mutation.")
        return HttpResponseRedirect("/charmming/energy/")

    try:
        cgws = structure.models.CGWorkingStructure.objects.get(workingstructure_ptr=ws.id)
        messages.error(request, "Mutation is not supported for coarse-grain models. If you wish to perform a point-mutation on this structure, please build a working structure with the CHARMM all-atom model.")
        return HttpResponseRedirect("/charmming/buildstruct/")
    except:
        pass

    segments = ws.segments.all()
    proCheck = False
    for segment in segments:
        if "good" in segment.name:
            messages.error(request, "GOOD hets are not currently supported by the mutation procedure. Please build a working structure with only PRO-type segments.")
            return HttpResponseRedirect('/charmming/buildstruct')
        if "pro" in segment.name:
            proCheck = True
    #If there's no "pro" segment, we shouldn't be mutating at all.
    if not proCheck:
        messages.error(request, "Your molecule does not have 'PRO' segments. CHARMMing only support mutating on 'PRO' segments at this time.")
        return HttpResponseRedirect('/charmming/buildstruct/')
    #This way we need less JS black magic
    tdict = {}
    tdict['proposedname'] = proposedname
    tdict['structname'] = struct.name
    filepath = struct.location.replace(charmming_config.user_home, '') + '/'
    filename = ws.identifier #This makes it so we can change the coordinates on-the-fly
#    ws = []
#    if len(existingWorkStructs) > 0:
#        tdict['haveworkingstruct'] = True
#        ws.extend(existingWorkStructs)
#        ws.sort()
#    else:
#        tdict['haveworkingstruct'] = False
#    tdict['built_list'] = ws
#Unsure if we need these since you can't swap working structs in mutation module ~VS
    # get list of all of the models that we can use
#    st_models = []
#    fp = open(struct.pickle, 'r')
#    pdb = cPickle.load(fp)
#    fp.close()
#    tdict['model_list'] = pdb.keys()
#Do we need the model? This seems like something we might actually want...
    amino_list = sorted(list(set(pychm.const.bio.aaAlphabet.keys())))
    tdict['amino_list']= amino_list
#    sl = []
#    sl.extend(structure.models.Segment.objects.filter(structure=struct,is_working='n'))
#    sl.sort()
# This seems excessive since in the iteration later down (iter_seg) we get all the segments anyway...
    segnames = []
    seglist = []
    if ws.localpickle:
        pdb = open(ws.localpickle)
    else:
        pdb = open(struct.pickle)
    #This way you can do multi-level mutations without issue
    pdbfile = cPickle.load(pdb)
    taco = pdbfile.iter_models().next() #we don't care about which model we load, so long as we load one
    pdb.close()
    tasks = Task.objects.filter(workstruct=ws,status='C',active='y',modifies_coordinates=True).exclude(action="solvation") #Things that were solvated on, you run at your own risk
    isMutated = False
    for task in tasks:
        if task.action == 'mutation':
            isMutated = True
            break
    tdict['filename'] = filename
    tdict['isMutated'] = isMutated
    #This way it gets replaced if there's mutation, or otherwise it'll just stay the same
    try:
        segments = ws.segments.all()
    except Exception as ex:
        return HttpResponse(str(ex)) #Something went very very wrong.
    segmentnames = map(lambda x:x.name,segments)
    #Pickle loading! Time to make this code even slower. We will use model00 since I honestly don't think the connectivity/segmentation of the molecule itself changes
    #with minimization and the like, it just changes conformation, and we don't need coords to do mutation (we can fetch those from elsewhere) 
    #segnames holds the segment names, seglist holds a list of dictionaries matching
    #the amino acids in each segment to their resIDs.
    #However we need to make SURE that only the segments from the current ws are loaded, or things can get quite bad.
    segnames = [] #We leave this in or otherwise the webpage gets the indices wrong
    chain_terminators = [] #Stores when each chain terminates (atom number); this way we don't have problems
    for seg in taco.iter_seg():
        if seg.segid in segmentnames and "pro" in seg.segid: #Careful...all proteins or non-purely-hetatm segments MUST have PRO in them!
            segdict = dict(zip(amino_list, [ [] for i in range (len(amino_list)) ]))
            atom = seg.iter_res().next().iter_atom().next().atomNum0 #Gets the (original) atom number to figure out the chain termination
            chain_terminators.append(atom)
            for res in seg.iter_res():
                try:
                    segdict[res.resName].append(res.resid)
                except KeyError:
                    pass
            segdict = dict (( (k,v) for k, v in segdict.iteritems() if v )) #Generator function
            seglist.append(segdict)
            segnames.append(seg.segid)
    dictlist = json.dumps(seglist) #Holds JSON strings of each dictionary so that JS can use them for shiny stuff
    tdict['segnames'] = segnames
    tdict['dictlist'] = dictlist
    tdict['chain_terminators'] = json.dumps(chain_terminators)
    tdict['tasks'] = tasks
    tdict['messages'] = messages.get_messages(request)
    tdict['ws_identifier'] = ws.identifier
    tdict['filepath'] = filepath
    tdict['super_user'] = request.user.is_superuser
    if (ws.isBuilt != 't'):
        messages.error(request, "Please perform a calculation on the whole atom set before performing a mutation.")
        return HttpResponseRedirect("/charmming/energy/")
#        tdict['no_coords'] = True

    return render_to_response('html/mutationselect.html', tdict)



