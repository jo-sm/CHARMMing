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

from django.template.loader import get_template
from django.contrib import messages
from django.contrib.messages import get_messages
from django.shortcuts import render_to_response
from structure.models import Structure, WorkingStructure, WorkingFile, Task
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from django.template import *
from pychm.lib.mol import Mol
import output, charmming_config, lessonaux, input, lessons, lesson1, lesson2, lesson3, lesson4
from account.views import checkPermissions
from pychm.scripts.getprop import getProp
#import structure.models
#import Structure
import os, re, time, cPickle

def doRMSD(postdata,location,id,picklefile,stlist):

    pfp = open(picklefile, 'r')
    pdb = cPickle.load(pfp)
    pfp.close()

    # get all of the Mol objects corresponding to the structures that we want RMSDs of
    mols = {}
    for t in stlist:
        if t.action == 'build':
            mols[t.action] = pdb['append_' + id]
        elif t.action == 'minimization':
            mols[t.action] = pdb['mini_' + id]
        elif t.action == 'solvation':
            mols[t.action] = pdb['neut_' + id]
        elif t.action == 'neutralization':
            mols[t.action] = pdb['neut_' + id]
        elif t.action == 'md':
            mols[t.action] = pdb['md_' + id]
        elif t.action == 'ld':
            mols[t.action] = pdb['ld_' + id]
        else:
            return output.returnSubmission('RMS Calculation', error='RMSD does not handle %s' % t.action)

    # Now, strip out crap like solvent and orient them
    # we do this by only grabbing pro, dna, and rna residues...
    smols = {}
    for k in mols.keys():
        tmp1 = mols[k].find(segtype='pro')
        tmp2 = mols[k].find(segtype='dna')
        tmp3 = mols[k].find(segtype='rna')
        smols[k] = Mol(sorted(tmp1 + tmp2 + tmp3))

    # Last step: if we want only the backbone, pull those atoms out of the first struct,
    # reorient it, and get the transform matrix to apply to all of the rest of the structs.
    if postdata['comp_allatom'] == '0':
        fmols = {}
        for k in smols.keys():
            finmol = Mol()
            if postdata.has_key('bb_n'):
                tmp = smols[k].find(atomtype=' n  ')
                finmol = Mol(sorted(finmol + tmp))
            if postdata.has_key('bb_ca'):
                tmp = smols[k].find(atomtype=' ca ')
                finmol = Mol(sorted(finmol + tmp))
            if postdata.has_key('bb_c'):
                tmp = smols[k].find(atomtype=' c  ')
                finmol = Mol(sorted(finmol + tmp))
            if postdata.has_key('bb_o'):
                tmp = smols[k].find(atomtype=' o  ')
                finmol = Mol(sorted(finmol + tmp))
            if postdata.has_key('bb_hn'):
                tmp = smols[k].find(atomtype=' hn ')
                finmol = Mol(sorted(finmol + tmp))
            if postdata.has_key('bb_ha'):
                tmp = smols[k].find(atomtype=' ha ')
                finmol = Mol(sorted(finmol + tmp))
            fmols[k] = Mol(sorted(finmol))
    else:
        fmols = smols

    # now go ahead and calculate the RMSDs:

    rmsdlines = ['<table class="table_center" border="1">\n']

    fp = open(location + '/rmsd-' + id + '.html', 'w')
    structlist = fmols.keys()
    structlist.sort()

    tabwidth = float(100.0/(len(structlist)+1))
    headline = '<tr><td style="width:%f%%"></td>' % tabwidth
    for k in structlist:
        headline += '<td style="width:%f%%">%s</td>' % (tabwidth,k)
    headline += '</tr>\n'
    rmsdlines.append(headline)

    for i in range(len(structlist)):
        thisline = '<tr><td>%s</td>' % structlist[i]
        for j in range(len(structlist)):
            if j < i:
                thisline += '<td></td>'
            elif j == i:
                thisline += '<td>0.00</td>'
            else:
                key1 = structlist[i]
                key2 = structlist[j]
                rmsd = fmols[key1].get_rmsd(fmols[key2], orient=True)
                thisline += '<td>%.2f</td>' % rmsd
        thisline += '</tr>\n'
        rmsdlines.append(thisline)
    rmsdlines.append('</table>\n')

    for line in rmsdlines:
        fp.write(line)
    fp.close()

    return render_to_response('html/rmsddisplay.html', {'rmsdlines': rmsdlines})

def rmsformdisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)

    #chooses the file based on if it is selected or not
    try:
         struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please upload a structure first.")
        return HttpResponseRedirect("/charmming/fileupload/")
#         return output.returnSubmission("RMS calculation", error="Please submit a structure first.")
    try:
         ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "Please build your structure first.")
        return HttpResponseRedirect("/charmming/buildstruct/")
#        return output.returnSubmission("RMS calculation", error="Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    feedback = ''
    tasks = Task.objects.filter(workstruct=ws,status='C',active='y',modifies_coordinates=True)
    tdict = {}
    tdict['ws_identifier'] = ws.identifier
    tdict['tasks'] = tasks
    rmsd_list = []
    for t in tasks:
        if request.POST.has_key('dotask_%d' % t.id):
            rmsd_list.append(t)
    if len(rmsd_list) > 0:
        if len(rmsd_list) == 1:
            messages.error(request, "More than 1 structure must be selected.")
            lesson_ok, dd_ok = checkPermissions(request)
            tdict['lesson_ok'] = lesson_ok
            tdict['dd_ok'] = dd_ok
            tdict['messages'] = messages.get_messages(request)
            return render_to_response('html/rmsdform.html', tdict)
#            return output.returnSubmission('RMS Calculation', error='More than 1 structure must be selected.')
        else:
            picklefile = struct.pickle
            if ws.localpickle != None:
                picklefile = ws.localpickle
            return doRMSD(request.POST,ws.structure.location,ws.identifier,picklefile,rmsd_list)
    else:
        prior_matrix = ''
        try:
            os.stat(struct.location + '/rmsd-' + ws.identifier + '.html')
        except:
            pass
        else:
            fp = open(struct.location + '/rmsd-' + ws.identifier + '.html')
            for line in fp:
                prior_matrix += line
            fp.close()
            tdict['prior_matrix'] = prior_matrix
    lesson_ok, dd_ok = checkPermissions(request)
    tdict['lesson_ok'] = lesson_ok
    tdict['dd_ok'] = dd_ok
    return render_to_response('html/rmsdform.html', tdict)

import dynamics.models
import traceback

def domdprop(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)
    #chooses the file based on if it is selected or not
    try: 
         struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
         return output.returnSubmission("MD properties", error="Please submit a structure first.")
    try:
         ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        return output.returnSubmission("MD properties", error="Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")
    try:
        mdp = dynamics.models.mdTask.objects.filter(workstruct=ws, status='C', active='y')[0]
    except:
        return output.returnSubmission('MD properties', error='Could not find completed MD task.')         


    fname = "%s/%s-md.out" % (ws.structure.location,ws.identifier)
#    cmdline = "%s/prop.py %s %s-mdproperties.dat" % (charmming_config.data_home,fname,ws.identifier)
    os.chdir(ws.structure.location)
    nprop = 0
    props = []
    if request.POST.has_key('getener'):
        props.append("averener")
#        cmdline += " averener"
    if request.POST.has_key('getvolu'):
        props.append("avervolu")
#        cmdline += " avervolu"
    if request.POST.has_key('getpressi'):
        props.append("averpressi")
#        cmdline += " averpressi"
    if request.POST.has_key('gettemp'):
        props.append("avertemp")
#        cmdline += " avertemp"
    if request.POST.has_key('gettote'):
        props.append("avertote")
#        cmdline += " avertote"
    if request.POST.has_key('gettotk'):
        props.append("avertotk")
#        cmdline += " avertotk"
    props.append("avertime")
    inpfile = open(fname)
    taco = getProp(inpfile,*props)
    inpfile.close()

    outLabels = []
    outData = []
    outLabels.append('avertime')
    outData.append(map(lambda x: '%10.5f' % x,taco['avertime']))
    for key in taco.keys():
        if key.endswith("time"):
            continue
        outLabels.append(key)
        outData.append(map(lambda x: '%10.5f' % x,taco[key]))
    outData = map(None, *outData)

    rawOutput = []
    rawOutput.append('    '.join(outLabels))
    for line in outData:
        rawOutput.append('    '.join(line))
    rawOutput = '\n'.join(rawOutput)

    outFile = open(ws.identifier + "-mdproperties.dat", "w")
    outFile.write(rawOutput)
    outFile.close()

#    logfp = open('/tmp/mdprop.txt', 'w')
#    logfp.write('cmdline = %s\n' % cmdline)
#    logfp.close()
#    os.system(cmdline)

    # create a working file, so that this thing appears on the download page
    wf = WorkingFile()
    wf.task = mdp
    wf.path = '%s/%s-mdproperties.dat' % (ws.structure.location,ws.identifier)
    wf.canonPath = wf.path
    wf.type = 'prp'
    wf.description = 'MD property analysis data'
    wf.save()
    tdict = {}
    lesson_ok, dd_ok = checkPermissions(request)
    tdict['lesson_ok'] = lesson_ok
    tdict['dd_ok'] = dd_ok
#    tdict['messages'] = get_messages(request)
#    return render_to_response("html/didprop.html", {'msg': "Properties generated, please download %s-mdproperties.dat from the download files page." % ws.identifier})
    tdict['outLabels'] = outLabels
    tdict['outData'] = outData
    return render_to_response("html/didprop.html", tdict)

def getmdprop(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)

    #chooses the file based on if it is selected or not
    try:
         struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
         return output.returnSubmission("MD properties", error="Please submit a structure first.")
    try:
         ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        return output.returnSubmission("MD properties", error="Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    try:
        mdp = dynamics.models.mdTask.objects.filter(workstruct=ws, status='C', active='y')[0]
    except:
        logfp = open('/tmp/getmdprop.txt', 'w')
        traceback.print_exc(file=logfp)
        logfp.close()
        return output.returnSubmission("MD properties", error="You must have completed molecular dynamics to get properties")

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/mdanalysis.html', {'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

