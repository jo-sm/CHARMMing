
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
from django.contrib import messages
from django.contrib.messages import get_messages
from django.shortcuts import render_to_response
from account.models import *
from structure.models import Task
from django.contrib.auth.models import User
from django.template import *
from account.views import checkPermissions
import structure.models
from lesson_config import *
import pychm.io, charmming_config
from pychm.io.pdb import PDBFile
from selection.models import AtomSelection, OniomSelection
from selection.models import LonePair
from atomselection_aux import saveAtomSelections, validateInputs
import cPickle, json
import traceback

def selectstructure(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    if request.POST:
        postdata=request.POST
        source = postdata['source']
        highest_qm_layer = 1 #TODO: This should be changed for MM/MM.
        if postdata.has_key('model_selected'): #If it doesn't have it going in from the energy/minim/whatever page...
            if postdata['model_selected'] not in ['qmmm','oniom']:
                messages.error(request, "Please use our graphical interface to perform atom selections. If you are seeing this message in error, please report it.")
                return HttpResponseRedirect("/charmming/"+dest)
            else:
                modelType = postdata['model_selected']
            try:
                layers = int(postdata['oniom_layers'])
            except:
                messages.error(request, "Invalid number of layers.")
                return HttpResponseRedirect("/charmming/"+dest)
            if modelType == "qmmm":
                layers = 1
            else:
                highest_qm_layer = layers - 1
                for i in range(2,layers):
                    if postdata.has_key('mm_box_layer_'+str(i)):
                        highest_qm_layer = i - 1
        elif postdata.has_key('modelType'): #Then we're in the selection page going back.
            if postdata['modelType'] not in ['qmmm','oniom']:
                messages.error(request, "Invalid model type. Please use our graphical interface to perform atom selections. If you are seeing this message in error, please report it.")
                return HttpResponseRedirect("/charmming/"+dest)
            else:
                modelType = postdata['modelType']
            try:
                layers = int(postdata['layers'])
            except:
                messages.error(request, "Invalid number of layers.")
                return HttpResponseRedirect("/charmming/"+dest)
        #If your modelType isn't oniom, don't worry about the layers.
        dest = source + "/" #Yes this is a bad naming convention. However, it works.
        try:
            struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
        except:
            messages.error(request, "No structure selected. Please upload a structure.")
            return HttpResponseRedirect("/charmming/fileupload/")
#            return HttpResponse("No structure selected.")

        try:
            ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
        except:
            logfp = open('/tmp/selectbarf.txt','w')
            traceback.print_exc(file=logfp)
            logfp.close()

            messages.error(request, "Please build a working structure.")
            return HttpResponseRedirect("/charmming/buildstruct/")
#           return HttpResponse("Please build a working structure.")

        try:
            task_id = postdata['task_id']
        except:
            messages.error(request, "No tasks present. Please perform a calculation on your working structure before doing QM/MM selection.")
            return HttpResponseRedirect("/charmming/energy/")

        try:
            task = Task.objects.filter(workstruct=ws,id=task_id)[0] #This might break with build tasks?
        except:
            messages.error(request, "Task " + task_id + " does not exist. Please report this bug.") #This means something went REALLY wrong in the database
            return HttpResponseRedirect("/charmming/"+dest)
#            return HttpResponse("Task " + task_id + " does not exist.")
        if ws.isBuilt == 'f':
            messages.error(equest, "Your working structure has not been built. Please perform a calculation on it before doing QM/MM selection.")
            return HttpResponseRedirect("/charmming/energy/")

        tdict = {}
        #Remember to validate input BEFORE going in to the rest.
        filepath = struct.location.replace(charmming_config.user_home,'') + "/" + ws.identifier + "-" + task.action + ".pdb"
        #for example 1yjp-ws1-build.pdb
        tdict['modelType'] = modelType
        tdict['layers'] = layers
        tdict['highest_qm_layer'] = highest_qm_layer
        tdict['filepath'] = filepath
        tdict['source'] = source
        tdict['task_id'] = task_id
        tdict['starting_layer'] = layers - 1 #This is not an issue with QM/MM since we override.
        logfp = open("/tmp/oniom.txt","w")
        logfp.write(str(tdict))
        logfp.close()
        #The following will all be stored in hidden fields and passed back and forth. We need to put checks on the receiving end (i.e. either the Python or the JS for qmmm_params)
        #such that no script injection can occur
        #It's safe to assume all these fields exist since they must by default
        if not (postdata.has_key('atomselection') or postdata.has_key('atomselection_layer_1')): #Make sure we're coming in from the right page.
            if layers > 1: #layers >1 => ONIOM model. This MAY require changing later if implementing more models.
                layer_values = [] #Hack because django doesn't take variables within variables, like {{exchange_layer{{layer}}}}
                for i in range(1,layers): #Since you don't care about the final layer...
                    current_layer = []
                    current_layer.append(i)
                    current_layer.append(postdata['qmmm_exchange_layer_'+str(i)])
                    current_layer.append(postdata['qmmm_charge_layer_'+str(i)])
                    current_layer.append(postdata['qmmm_correlation_layer_'+str(i)])
                    current_layer.append(postdata['qmmm_basisset_layer_'+str(i)])
                    current_layer.append(postdata['qmmm_multiplicity_layer_'+str(i)])
                    layer_values.append(current_layer) #Therefore layer 1 is stored in 0, but contains the number, so we can template easily.
                    #Assume user does not input a text selection, but also make sure these values don't change.
                    #TODO: Make sure that when the form is submitted for any QM/MM calc that these values are checked and updated on the AtomSelection/OniomSelection.
                tdict['layer_values'] = layer_values
            else:
                    tdict['exchange'] = postdata['qmmm_exchange']
                    tdict['charge'] = postdata['qmmm_charge']
                    tdict['correlation'] = postdata['qmmm_correlation']
                    tdict['basis_set'] = postdata['qmmm_basisset']
                    tdict['multiplicity'] = postdata['qmmm_multiplicity']
        else:
            validate_inputs_result = validateInputs(tdict,postdata,layers)
            if validate_inputs_result != "Passed":
                messages.error(request,validate_inputs_result)
                return HttpResponseRedirect('/charmming/'+dest)
            error_message = saveAtomSelections(request,ws,task)
            try:
                arf = int(error_message)
            except: #If it's not an int, we have a problem
                messages.error(request, error_message)
                return HttpResponseRedirect('/charmming/'+dest)
            return HttpResponseRedirect('/charmming/'+ dest)
            #Do stuff with the database and return redirect to source
        #TODO: Consider rewriting the above behemoth into a generic method like saveAtomSelections in etc.atomselection_aux.py
        #This code is applied when there is NOT an atomselection active.
        if ws.localpickle:
            pdb = open(ws.localpickle)
            pdbfile = cPickle.load(pdb)
            pdb.close()
        else:
            pdbfile = PDBFile(struct.location + "/" + ws.identifier + "-" + task.action + ".pdb")
        taco = pdbfile['model00']
        #Important note: this does NOT work the same as in mutation. Since the segment iteration
        #Goes by alphabet rather than by atom number order, this becomes a mess. So we
        #Get the segment objects in the list, then sort by their "first atom" number. HOpefully
        #The fact that we make a new generator object means we don't repeat ourselves and thus
        #have issues with single-residue things.
        segmentlist = [] #Holds the segment names
        chain_terminators = [] #Holds the "chain terminator". We pair them up in different lists because we need to preserve order
        for seg in taco.iter_seg(): #segments are in alphabetical order...
            atom = seg.iter_res().next().iter_atom().next().atomNum0 #Gets the original atom number of the first atom in the first residue of this segment
            segmentlist.append(seg)
            chain_terminators.append(atom)
        #The next line makes it so that the segments are sorted by the atom number of their first residue.
        #While we do get this same thing while iterating, the segments are in alphabetical order rather than in order by atomic number.
        #Therefore, we need to sort them by the number of that atom once they're done, but we don't really "know" the order until the end anyway.
        #Thus, the sort below.
        segmentlist = sorted(segmentlist,key=lambda x:x.iter_res().next().iter_atom().next().atomNum0) #This is to get the thing properly sorted...
        segmentlist = map(lambda x:x.segid,segmentlist) #This is to get it back to names...
        chain_terminators = sorted(chain_terminators)
        tdict['chain_terminators'] = json.dumps(chain_terminators)
        tdict['segmentlist'] = json.dumps(segmentlist)
        return render_to_response("html/selection.html", tdict)
    else:
        messages.error(request, "Atom selection must be accessed via a calculation page. No such source has been found.")
        return HttpResponseRedirect("/charmming/about/")
#        return HttpResponse("No source.")
