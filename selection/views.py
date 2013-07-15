
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
from selection.models import AtomSelection
from selection.models import LonePair
import cPickle, json

def selectstructure(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    if request.POST:
        postdata=request.POST
        source = postdata['source']
        dest = source + "/" #Yes this is a bad naming convention. However, it works.
        try:
            struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
        except:
            messages.error(request, "No structure selected. Please upload a structure.")
            return HttpResponseRedirect("/charmming/fileupload/")
#            return HttpResponse("No structure selected.")

        try:
            ws = structure.models.WorkingStructure.objects.get(structure=struct,selected='y')
        except:
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
        filepath = struct.location.replace(charmming_config.user_home,'') + "/" + ws.identifier + "-" + task.action + ".pdb"
        #for example 1yjp-ws1-build.pdb
        tdict['filepath'] = filepath
        tdict['source'] = source
        tdict['task_id'] = task_id
        if postdata.has_key('atomselection'):
            try:
                foo = postdata['atomselection']
                specialchars = set("#$/;\n\\_+=[]{}()&^%")
                if not(foo.startswith("bynum")) or (len(specialchars.intersection(foo)) > 0): #i.e. no special chars
                    messages.error(request, "Invalid selection string. Please use the automated selection system.")
                    tdict['messages'] = get_messages(request)
                    return render_to_response("html/selection.html", tdict)
            except: #I'm not sure what would break here, but let's do it anyway.
                messages.error(request, "You have selected too many atoms, or performed an invalid selection. Please try again.")
                tdict['messages'] = get_messages(request)
                return render_to_respose("html/selection.html", tdict)
            try:
               tdict['num_linkatoms'] = int(postdata['linkatom_num']) #This will fail if the user put in anything other than a number because of how int works.
            except:
                messages.error(request, "Number of link atoms is not an integer.")
                tdict['messages'] = get_messages(request)
                return render_to_response("html/selection.html", tdict)
            qmhosts = []
            mmhosts = []
            #TODO: Switch to the new DB structure
            for key in postdata.keys():
                #Thankfully we don't have to modify anything here to add in the atom number support
                if key.startswith("qmhost"):
                    if len(specialchars.intersection(request.POST[key])) == 0:
                        foo = postdata[key].split("\t")
                        qmhosts.append((str(key[-1]),foo[0],foo[1],foo[2]))
                        #The format above goes as follows:
                        #(divid,resid,atomtype,segid). 
                        #As describe in selection.models, we only need divid once to make the correct orderings.
                    else:
                        messages.error(request, "Invalid QM atom selection. Please use the automated selection system.")
                        tdict['messages'] = get_messages(request)
                        return render_to_response('html/selection.html', tdict)
                elif key.startswith("mmhost"):
                    if len(specialchars.intersection(request.POST[key])) == 0:
                        foo = postdata[key].split("\t")
                        mmhosts.append((str(key[-1]),foo[0],foo[1],foo[2]))
                    else:
                        messages.error(request, "Invalid MM atom selection. Please use the automated selection system.")
                        tdict['messages']= get_messages(request)
                        return render_to_response('html/selection.html', tdict)

            try: #First see if there's already an atomselection present
                oldselect = AtomSelection.objects.filter(workstruct=ws)[0]
                oldselect.task = task
                oldselect.workstruct = ws #Just in case...
                oldselect.selectionstring = postdata['atomselection'] #replace it
                for i in range(0, len(qmhosts)):
                    lonepair = LonePair()
                    lonepair.selection = oldselect
                    lonepair.divid = qmhosts[i][0]
                    lonepair.qmresid = qmhosts[i][1]
                    lonepair.qmatomname = qmhosts[i][2]
                    lonepair.qmsegid = qmhosts[i][3]
                    lonepair.mmresid = mmhosts[i][1]
                    lonepair.mmatomname = mmhosts[i][2]
                    lonepair.mmsegid = mmhosts[i][3]
                    lonepair.save()
                oldselect.linkatom_num = postdata['linkatom_num']
                oldselect.save()
            except: #There's two exceptions possible, so first we try to make a new atomselect
                atomselect = AtomSelection()
                atomselect.task = task
                atomselect.workstruct = ws
                try:
                    atomselect.selectionstring = postdata['atomselection']
                    atomselect.linkatom_num = postdata['linkatom_num']
                    atomselect.save()
                except Exception as ex:
                    return HttpResponse("Woof" + str(ex) + "\n" + str(lonepair.selection) + "\t" + str(lonepair.divid))
                try:
                    for i in range(0, len(qmhosts)):
                        lonepair = LonePair()
                        lonepair.selection = atomselect
                        lonepair.divid = qmhosts[i][0]
                        lonepair.qmresid = qmhosts[i][1]
                        lonepair.qmatomname = qmhosts[i][2]
                        lonepair.qmsegid = qmhosts[i][3]
                        lonepair.mmresid = mmhosts[i][1]
                        lonepair.mmatomname = mmhosts[i][2]
                        lonepair.mmsegid = mmhosts[i][3]
                        lonepair.save()
                    atomselect.save()
                except Exception as ex: #THis is if you have so many atoms you break the postdata request size or the maximum char field size
                    return HttpResponse("Woof" + str(ex) + "\n" + str(lonepair.selection) + "\t" + str(lonepair.divid))
#                    messages.error(request, "You have selected too many atoms, or performed an invalid selection. Please try again.")
#                    tdict['messages'] = get_messages(request)
#                    return render_to_response("html/selection.html", tdict)
#                    return HttpResponse("You have selected too many atoms.")
            return HttpResponseRedirect('/charmming/'+ dest)
            #Do stuff with the database and return redirect to source
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
