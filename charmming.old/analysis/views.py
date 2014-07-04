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
from django.http import HttpResponseRedirect
from structure.models import Structure, WorkingStructure, CGWorkingStructure, WorkingFile, Task
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from django.template import *
from pychm.lib.mol import Mol
from pychm.tools import expandPath, flatten
from pychm.io.pdb import get_molFromPDB
from pychm.future.io.charmm.dcd import open_dcd
import numpy as np
import output, charmming_config, lessonaux, input, lessons, lesson1, lesson2, lesson3, lesson4, lesson5
from account.views import checkPermissions
from pychm.scripts.getprop import getProp
#import structure.models
#import Structure
import os, re, time, cPickle
import pychm
import traceback

def doRMSD(user,postdata,location,id,picklefile,stlist):

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
        elif t.action == 'sgld':
            mols[t.action] = pdb['sgld_' + id]
        else:
            return output.returnSubmission('RMS Calculation', error='RMSD does not handle %s' % t.action)

    # Now, strip out crap like solvent and orient them
    # we do this by only grabbing pro, dna, and rna residues...
    smols = {}
    for k in mols.keys():
        tmp1 = mols[k].find(segtype='pro')
        tmp2 = mols[k].find(segtype='dna')
        tmp3 = mols[k].find(segtype='rna')
        tmp4 = mols[k].find(segtype='go')
        tmp5 = mols[k].find(segtype='bln')
        tmp6 = mols[k].find(segtype='bad')
        smols[k] = Mol(sorted(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6))

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

    # lesson time...
    try:
         struct = Structure.objects.filter(owner=user,selected='y')[0]
    except:
        return HttpResponseRedirect("/charmming/fileupload/")
    try:
         ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        return HttpResponseRedirect("/charmming/buildstruct/")

    try:
        lnum=ws.structure.lesson_type
        lesson_obj = eval(lnum+'.models.'+lnum.capitalize()+'()')
    except:
        lesson_obj = None

    if lesson_obj:
        lessonaux.doLessonAct(ws.structure,"onRMSDSubmit",postdata)
    # end lesson time

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
            return doRMSD(request.user,request.POST,ws.structure.location,ws.identifier,picklefile,rmsd_list)
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

## aux routines for NatQ
def get_native_contacts(mymol, AACUT=4.5):
    # build AA.resid -> CG.atomid map
    mymap = [ None, ]
    ngly = 0
    for i, res in ( (res.resid, res) for res in mymol.iter_res()):
        if res.resName == 'gly':
            ngly += 1
            mymap.append(None)
        cur_index = i*2 + ngly
        mymap.append(cur_index)
    # nuke backbone atoms
    atoms = ( atom for atom in mymol if not atom.is_backbone() )
    # nuke hydrogens
    atoms = ( atom for atom in atoms if atom.element != 'h' )
    # rebuild Mol obj
    mymol = Mol(atoms)
    # build list of residue pairs in native contact
    native_contacts = []
    iterres = [ (res.resid, res) for res in mymol.iter_res() ]
    for i, res_i in iterres:
        for j, res_j in iterres:
            if i >= j-1:
                continue # dont consider adjacent residues
            try:
                for atom_i in res_i:
                    for atom_j in res_j:
                        r_ij = atom_i.calc_length(atom_j)
                        if r_ij <= AACUT:
                            native_contacts.append((i, j))
                            raise AssertionError # break a double loop
            except AssertionError:
                pass
    # apply AA.resid -> CG.atomid map to native_contacts
    def fun(inp):
        return (mymap[inp[0]], mymap[inp[1]])
    #
    return map(fun, native_contacts)

def extract_xyz(mydcd, native_contacts):

    logfp = open('/tmp/extract_xyz.txt','w')

    # dump all data
    with open_dcd(mydcd) as fp:
        data = fp.get_massive_dump()
    # figure out which cg atoms we need
    logfp.write("len data = %d %d %d\n\n\n" % (len(data['x']),len(data['y']),len(data['z'])))
    cg_ids = sorted(set(flatten(native_contacts)))
    logfp.write('cg_ids = ' + str(cg_ids) + '\n')
    logfp.flush()
    # extract the xyz data from the dump, for our coordinates of interest
    x = {}
    y = {}
    z = {}
    for id in cg_ids:
        try:
            logfp.write('id = %s\n' % id)
            x[id] = data['x'][:,id-1]
            logfp.write('got x\n')
            y[id] = data['y'][:,id-1]
            logfp.write('got y\n')
            z[id] = data['z'][:,id-1]
            logfp.write('got z\n')
        except:
            logfp.write('got exception ... skipping\n')
            logfp.flush()

    logfp.close()
    return (x, y, z)

def calc_q(native_contacts, x, y, z, CGCUT = 8.0):
    # look at cut**2, instead of cut, sqrt's are expensive
    CUT = CGCUT * CGCUT
    # build scratch space for r**2
    dim0 = len(native_contacts)
    dim1 = x[native_contacts[0][0]].shape[0]
    data = np.zeros((dim0, dim1))
    # populate scratch with r**2 values for all native contacts
    for n, (i, j) in enumerate(native_contacts):
        dx = x[i]-x[j]
        dy = y[i]-y[j]
        dz = z[i]-z[j]
        rr = (dx*dx+dy*dy+dz*dz)
        data[n, :] = rr
    # build scratch for Q
    tmp = np.zeros(data.shape)
    # use fancy numpy mapping to check if R**2 is .gt. or .lt. CUT**2
    tmp[data <= CUT] = 1
    tmp[data > CUT] = 0
    # return Q(T)
    return tmp.mean(axis=0)

## main routine for native contacts...
def natqdisplay(request):
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

    try:
        cgws = CGWorkingStructure.objects.get(workingstructure_ptr=ws.id)
    except:
        logfp = open('/tmp/cdws.txt', 'w')
        traceback.print_exc(file=logfp)
        logfp.close()

        messages.error(request, "Native contact analysis is only available for CG models")
        return HttpResponseRedirect("/charmming/buildstruct/")
    if cgws.cg_type != 'go':
        messages.error(request, "Native contact analysis is only available for KT Go models")
        return HttpResponseRedirect("/charmming/buildstruct/")

    pdbfile = struct.location + '/cgbase-' + str(cgws.id) + '.pdb'
    dcdfile = struct.location + '/' + ws.identifier + '-md.dcd'
    try:
        os.stat(pdbfile)
    except:
        messages.error(request,"Please build your structure and run dynamics first")
        return HttpResponseRedirect("/charmming/buildstruct/")
        
    baddcd = False
    try:
        os.stat(dcdfile)
    except:
        baddcd = True

    if baddcd:
        dcdfile = struct.location + '/' + ws.identifier + '-ld.dcd'
        try:
            os.stat(dcdfile)
        except:
            messages.error(request,"Please run dynamics first")
            return HttpResponseRedirect("/charmming/buildstruct/")


    natqlines = []

    # --- Frank's main routine here ---
    my_aa_mol = get_molFromPDB(pdbfile)
    my_aa_mol.parse()

    # scan the AA Mol object, to build a list of residues that are in native
    # contact in the crystal structure.
    my_native_contacts = get_native_contacts(my_aa_mol)
    # read the xyz data from the DCD file, throw away the coordinates that
    # arent part of a native contact pair
    x, y, z = extract_xyz(dcdfile,my_native_contacts)
    q = calc_q(my_native_contacts, x, y, z)

    # we only have 1000 frames in a CHARMMing trajectory file by default.
    # One frame every 100 steps
    mystep = 100
    for i in range(len(q)):
        natqlines.append((mystep,'%5.3f' % q[i]))
        mystep += 100

    logfp = open('/tmp/l5thing.txt','w')
    logfp.write('at end\n')
    try:
        lnum=ws.structure.lesson_type
        lesson_obj = eval(lnum+'.models.'+lnum.capitalize()+'()')
    except:
        lesson_obj = None

    logfp.write('lesson_obj = %s\n' % str(lesson_obj))
    logfp.close()

    if lesson_obj:
        lessonaux.doLessonAct(ws.structure,"onNATQSubmit",request.POST)

    return render_to_response('html/natqdisplay.html', {'natqlines': natqlines})

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

