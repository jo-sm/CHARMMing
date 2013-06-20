
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
from django.template import *
from account.views import checkPermissions
from lessons.models import LessonProblem
from mutation.models import mutateTask
import output, lesson1, lesson2, lesson3, lesson4, lessonaux
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
        return HttpResponse("No structure selected.") #I have no idea why this would happen but you never know.


    try:
        ws = structure.models.WorkingStructure.objects.get(structure=struct,selected='y')
    except:
        return HttpResponse("You must first build a working structure before applying mutations.")

    try:
        os.stat(workstruct.structure.location + "/" + postdata['MutFile'] + "-mutation.pdb")
        return HttpResponse("Structure with that name already exists.")
    except: #Basically if the file DOESN'T exist, keep going
        pass

    MutResi = postdata['MutResi']
    MutResName = postdata['MutResName'] #The script doesn't actually care about this one, it's just for user reference
    NewResn = postdata['NewResn']
    MutSegi = postdata['MutSegi']
    MutFile = postdata['MutFile']

    modstruct = copy.deepcopy(ws)
    modstruct.id = None
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

    if modstruct.localpickle:
        woof = modstruct.localpickle
    else:
        woof = struct.pickle #No specific pickle file, possibly not mutated before

    mt = mutateTask()
    mt.setup(modstruct)
    mt.localpickle = woof
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
    logfp.write("localpickle = " + str(mt.localpickle) + "\n")
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
    if ws.isBuilt != 't':
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
    mutate_filename = ws.structure.location + "/" + ws.identifier + "-mutate.inp"
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
        return HttpResponse("No structures present to mutate.")

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
        return HttpResponse("You must first build a working structure.")
    try:
        ws = structure.models.WorkingStructure.objects.get(structure=struct,selected='y')
    except:
        return HttpResponse("You must first build a working structure.")
    #This way we need less JS black magic
    tdict = {}
    tdict['proposedname'] = proposedname
    tdict['structname'] = struct.name
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
    taco = pdbfile['model00']
    pdb.close()
    tasks = Task.objects.filter(workstruct=ws,status='C',active='y').exclude(action='energy')
    try:
        segments = ws.segments.all()
    except Exception as ex:
        return HttpResponse(str(ex))
    segmentnames = map(lambda x:x.name,segments)
    #Pickle loading! Time to make this code even slower. We will use model00 since I honestly don't think the connectivity/segmentation of the molecule itself changes
    #with minimization and the like, it just changes conformation, and we don't need coords to do mutation (we can fetch those from elsewhere) 
    #segnames holds the segment names, seglist holds a list of dictionaries matching
    #the amino acids in each segment to their resIDs.
    #However we need to make SURE that only the segments from the current ws are loaded, or things can get quite bad.
    segnames = [] #We leave this in or otherwise the webpage gets the indices wrong
    for seg in taco.iter_seg():
        if seg.segid in segmentnames and "pro" in seg.segid: #Careful...all proteins or non-purely-hetatm segments MUST have PRO in them!
            segdict = dict(zip(amino_list, [ [] for i in range (len(amino_list)) ]))
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
    tdict['tasks'] = tasks
#    tdict['disulfide_list'] = struct.getDisulfideList() #NO idea if we even need this, or proto_list, it's just more cPickle open reads and those make everything slow
#    tdict['proto_list'] = []
    tdict['super_user'] = request.user.is_superuser
#    for seg in sl:
#        tdict['proto_list'].extend(seg.getProtonizableResidues()) #Is this necessary?

    tdict['lesson_ok'], tdict['dd_ok'] = checkPermissions(request)
    return render_to_response('html/mutationselect.html', tdict)


