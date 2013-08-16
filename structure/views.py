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
# Create your views here.
from httplib import HTTPConnection
from django import forms
from django.db.models import Q
from django.template import RequestContext
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.contrib import messages
from django.contrib.messages import get_messages #These are for fancy redirects...
from minimization.models import minimizeTask
from dynamics.models import mdTask, ldTask, sgldTask
from solvation.models import solvationTask
from account.models import *
from normalmodes.views import combineNmaPDBsForMovie
from normalmodes.aux import getNormalModeMovieNum
from normalmodes.models import nmodeTask
import dd_infrastructure.models
from apbs.models import redoxTask
from structure.models import Task, energyTask
from django.contrib.auth.models import User
from django.core.mail import mail_admins
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from account.views import checkPermissions
from structure.aux import checkNterPatch
from lessons.models import LessonProblem
from structure.qmmm import makeQChem_tpl, makeQchem_val
from atomselection_aux import getAtomSelections, saveAtomSelections
import output, lesson1, lesson2, lesson3, lesson4, lessonaux
import structure.models, input
import selection.models
import structure.mscale
import subprocess
# import all created lessons by importing lesson_config
# also there is a dictionary called 'file_type' in lesson_config.py specififying the file type of files uploaded by the lessons  
from lesson_config import *
import os, sys, re, copy, datetime, time, stat, json, openbabel
import mimetypes, string, random, glob, traceback
import pychm.io, charmming_config, minimization
import cPickle


# problem during upload
"""
def uploadError(request,problem):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    prob_string = "An unknown error occurred during file upload."
    if problem == 'overflow':
       prob_string = "The file you're trying to upload has more than 50,000 atoms."
    elif problem == 'parse_error':
       prob_string = "There was a problem trying to parse your structure file."
    foo = open(charmming_config.charmming_root+"/mytemplates/html/problem.html")
    html_request = foo.read()
    foo.close()
    html_request.replace("prob_string",prob_string)
    messages.error(request, html_request)
    return HttpResponseRedirect("/charmming/fileupload/")
#    return render_to_response('html/problem.html',{'prob_string': prob_string})
The above does not appear to work.
We will simply redirect - this has no point anymore.
"""

# Deletes user specified file from database and files relating to it
def deleteFile(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    #make sure request data isn't malicious
    input.checkRequestData(request)
    if request.POST.has_key('filename'):
        delete_filename = request.POST['filename']
    else:
        messages.error(request, "Bad")
        return HttpResponseRedirect("/charmming/buildstruct/")

    deleteAll = 0
    remove_extension = re.compile('\.\S*')
    #If user wants to remove all files, visit /deletefile/all_files/

    if(delete_filename == "all_files"):
        deleteAll = 1
        total_files = structure.models.Structure.objects.filter(owner=request.user)
    else:
        file = structure.models.Structure.objects.filter(owner=request.user,name=delete_filename)[0]
        total_files = [file]

    logfp = open('/tmp/delfile.txt','w')
    logfp.write('In delete\n')
    logfp.flush()

    for s in total_files[::-1]:
        os.chdir(charmming_config.user_home + '/' + request.user.username)
        subprocess.call(["rm","-rf",s.name])
#        os.system("rm -rf " + s.name)

        if s.selected == 'y' and deleteAll == 0:
            s.selected = ''
            allFileList = structure.models.Structure.objects.filter(owner=request.user)
            if len(allFileList) >= 2:
                # We need to change the selection to another file, arbitrarilty we choose
                # the first one.
                for candidateSwitch in allFileList:
                     if candidateSwitch.name != delete_filename:
                         candidateSwitch.selected = 'y'
                         candidateSwitch.save()
                         break

        try:
            s.delete()
        except:
            traceback.print_exc(file=logfp)

    logfp.close()
    return HttpResponse('Done')


#Lets user view file in browser
def viewFiles(request,filename, mimetype = None):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    username = request.user.username

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "No structure has been uploaded. Please upload a structure before downloading its files.")
        return HttpResponseRedirect('/charmming/fileupload/')
    try:
#        logfp = open('/tmp/vf.txt', 'w')
#        logfp.write("%s/%s\n" % (struct.location,filename))
#        logfp.close()
        os.stat("%s/%s" % (struct.location,filename))
    except:
        messages.error(request, "The file " + filename + " doesn't seem to exist. Maybe you deleted it?")
        return HttpResponseRedirect('/charmming/viewprocessfiles/')

    fname = "%s/%s" % (struct.location,filename)
    fp = open(fname, 'r')
    fcontent = fp.read()
    fp.close()

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/viewfile.html', {'fcontent': fcontent, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

#Wraps all files associated with the PDB into a tar and deletes it
def downloadFilesPage(request,mimetype=None):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return render_to_response('html/nopdbuploaded.html')
    try:
        ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "Please build a working structure first.")
        return HttpResponseRedirect("/charmming/buildstruct/")
    try:
        cgws = structure.models.CGWorkingStructure.objects.filter(workingstructure_ptr=ws.id)
    except structure.models.CGWorkingStructure.MultipleObjectsReturned:
        messages.error(request, "More than one CGWorkingStructure is associated to this working structure. Please report this bug, then rebuild your structure.")
        return HttpResponseRedirect("/charmming/fileupload/")
    except:
        pass


    wsname = ws.identifier
    filelst = []
    headers = []

    # files from segment construction
    headers.append('Segment setup')
    if not cgws:
        for wseg in ws.segments.all():
            filelst.append(('Segment setup', wseg.builtPSF, 'PSF of segment %s' % wseg.name))
            filelst.append(('Segment setup', wseg.builtPSF.replace('.psf','.pdb'), 'PDB of segment %s' % wseg.name))
            filelst.append(('Segment setup', wseg.builtCRD, 'CRD of segment %s' % wseg.name))

            # there would be an input and output too
            filelst.append(('Segment setup', 'build-%s.inp' % wseg.name, 'Input file to build segment %s' % wseg.name))
            filelst.append(('Segment setup', 'build-%s.out' % wseg.name, 'Output file from building segment %s' % wseg.name))
    else:
        filelst.append(('Segment setup', ws.structure.name + "-go.prm", "PRM setup for " + ws.identifier))
        filelst.append(('Segment setup', ws.structure.name + "-go.pdb", "PDB setup for " + ws.identifier))
        filelst.append(('Segment setup', ws.structure.name + "-go.rtf", "RTF setup for " + ws.identifier))


    # now get all Tasks associated with the structure and get their workingFiles
    # hack alert: we only store CRDs now but we can extrapolate the PSFs, INPs, and OUTs
    # ToDo: store this info in the workingfile
    tasks = structure.models.Task.objects.filter(workstruct=ws,active='y')
    for task in tasks:
        headers.append(task.action)
        # get a list of all of the files associated with this task
        files = structure.models.WorkingFile.objects.filter(task=task)
        for file in files:
            basename = file.canonPath.split('/')[-1]
            filelst.append((task.action, basename, file.description))
    return render_to_response('html/downloadfiles.html', {'filelst': filelst, \
                                                          'headers': headers, \
                                                          'wsname': wsname})

#Wraps all files associated with the PDB into a tar and deletes it
def downloadTarFile(request,mimetype=None):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    tar_filename = struct.name + '.tar.gz'

    username = request.user.username
    os.chdir(charmming_config.user_home + '/' + username)

    try:
        os.unlink(tar_filename)
    except:
        pass

    # remove condor logs and error files
    #Now with 100% more subprocess
    subprocess.call(["rm","-f",(struct.name + '/*.{err,log}')])
#    os.system('rm -f ' + struct.name + '/*.{err,log}')
    subprocess.call(["cp", (charmming_config.data_home + '/solvation/water.crd'), struct.name])
    subprocess.call(["cp", (charmming_config.data_home + '/calcewald.pl'), struct.name])
    subprocess.call(["cp", (charmming_config.data_home + '/savegv.py'), struct.name])
    subprocess.call(["cp", (charmming_config.data_home + '/savechrg.py'), struct.name])
#    os.system('cp ' + charmming_config.data_home + '/solvation/water.crd ' + struct.name)
#    os.system('cp ' + charmming_config.data_home + '/calcewald.pl ' + struct.name)
#    os.system('cp ' + charmming_config.data_home + '/savegv.py ' + struct.name)
#    os.system('cp ' + charmming_config.data_home + '/savechrg.py ' + struct.name)
    subprocess.call(["cp",(charmming_config.data_home + '/toppar/top_all36_prot.rtf'),(charmming_config.data_home + '/toppar/par_all36_prot.prm'),struct.name])
    subprocess.call(["tar","-czf",tar_filename,struct.name])
#    os.system('cp ' + charmming_config.data_home + '/toppar/top_all36_prot.rtf ' + charmming_config.data_home + '/toppar/par_all36_prot.prm ' + struct.name)
#    os.system('tar -czf ' + tar_filename + ' ' + struct.name)
    statinfo = os.stat(tar_filename)
    if mimetype is None:
        mimetype,encoding = mimetypes.guess_type("%s/%s/%s" % (charmming_config.user_home,username,tar_filename))
    response = HttpResponse(mimetype=mimetype)
    response['Content-Disposition'] = 'attachment; filename=%s' % tar_filename
    response['Content-Length'] = statinfo.st_size
    response.write(file(charmming_config.user_home + "/" + username + "/" + tar_filename, "rb").read())
    return response

def downloadGenericFiles(request,filename, mimetype = None):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    username = request.user.username

    fnarr = filename.split("/")
    filename = fnarr[-1]

    path = "/usr/local/charmming"

    try:
        os.stat(path + "/" + filename)
    except:
        messages.error(request, "Hmmm ... the file " + filename + " doesn't seem to exist. Please file a bug report.")
        return HttpResponseRedirect('/charmming/viewprocessfiles/')

    if mimetype is None:
        mimetype,encoding = mimetypes.guess_type(path + "/" + filename)
    response = HttpResponse(mimetype=mimetype)
    response['Content-Disposition'] = 'attachment; filename=%s' %filename
    response.write(file(path + "/" + filename, "rb").read())
    return response
"""
##TODO: Find out why this method exists. Nothing calls it. Why is it here?
def viewDynamicOutputContainer(request,jobtype):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    file =  Structure.objects.filter(owner=request.user,selected='y')[0]
    input.checkRequestData(request)
    input.checkForMaliciousCode(jobtype,request)
    if jobtype in ['minimization']:
        charmm_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-min.out'
        try:
            os.stat(file.location + charmm_filename)
        except:
            return HttpResponse('Please wait while the calculation starts.')
    elif jobtype in ['solvation']:
        charmm_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-solv.out'
        try:
            os.stat(file.location + charmm_filename)
        except:
            return HttpResponse('Please wait while the calculation starts.')
        try:
            os.stat(file.location + 'charmm-' + file.stripDotPDB(file.filename) + '-neutralize.out')
            charmm_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-neutralize.out'
        except:
            pass
    elif jobtype in ['nma']:
        charmm_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-nmodes.out'
        try:
            os.stat(file.location + charmm_filename)
        except:
            return HttpResponse('Please wait while the calculation starts.')
    elif jobtype in ['md']:
        charmm_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-md.out'
        try:
            os.stat(file.location + charmm_filename)
        except:
            return HttpResponse('Please wait while the calculation starts.')
    elif jobtype in ['ld']:
        charmm_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-ld.out'
        try:
            os.stat(file.location + charmm_filename)
        except:
            return HttpResponse('Please wait while the calculation starts.')
    elif jobtype in ['sgld']:
        charmm_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-sgld.out'
        try:
            os.stat(file.location + charmm_filename)
        except:
            return HttpResponse('Please wait while the calculation starts.')
    else:
        return HttpResponse('Dynamic output is not available for this job type!')
    return render_to_response('html/viewdynamicoutput.html',{'jobtype':jobtype,'charmm_filename':charmm_filename})


def viewDynamicOutput(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    file =  Structure.objects.filter(owner=request.user,selected='y')[0]
    input.checkRequestData(request)
    inputname = request.POST['inputname']
    return downloadProcessFiles(request,inputname)
"""

#Lets user download the actual file
def downloadProcessFiles(request,filename, mimetype = None):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    username = request.user.username

    #if the filename contains a slash, it is probably not located in the pdb_uploads/ folder
    #an example would be filename == /solvation/water.crd which is located in /usr/local/charmming
    slash = re.compile("/")
    if (slash.search(filename)):
        return downloadGenericFiles(request,filename)
    if mimetype is None:
        mimetype,encoding = mimetypes.guess_type("%s/%s/%s" % (charmming_config.user_home,username,filename))

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "No structure uploaded. Please upload a structure.")
        return HttpResponseRedirect('/charmming/fileupload/')
#        return HttpResponse('No structure uploaded')

    try:
        sval = os.stat("%s/%s" % (struct.location,filename))
    except:
        messages.error(request, "Oops...the file " + filename + " no longer exists.")
        return HttpResponseRedirect("/charmming/viewprocessfiles/")
#        return HttpResponse('Oops ... that file no longer exists.')
    response = HttpResponse(mimetype=mimetype)
    response['Content-Disposition'] = 'attachment; filename=%s' % filename
    response['Content-length'] = sval.st_size
    response.write(file("%s/%s" % (struct.location,filename), "rb").read())
    return response

#Loads template that contains all the inp/out files for download
def viewProcessFiles(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "No structure uploaded. Please upload a structure.")
        return HttpResponseRedirect('/charmming/fileupload/')
#        return HttpResponse("Please select a structure first.")
    try:
        ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "Please build a working structure first.")
        return HttpResponseRedirect('/charmming/buildstruct/')
#        return HttpResponse('No structure uploaded')

    header_list = []
    file_list   = []
    tasks = structure.models.Task.objects.filter(workstruct=ws)
    for task in tasks:
        if task.action not in header_list: #No more than one task per ws...this works fine with mutation since that creates a new WS
            header_list.append(task.action)

        # get all input and output files associated with the task
        wfiles = structure.models.WorkingFile.objects.filter( Q(task=task), \
                   Q(type='inp') | Q(type='out'))

        for wf in wfiles:
            bn = wf.canonPath.split('/')[-1]
            ds = wf.description
            file_list.append((task.action,bn,ds))

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/viewprocessfiles.html', {'headers': header_list, 'files': file_list, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok, 'messages':get_messages(request)})

def visualize(request,filename):
    """
    Allows the user to visualize their structure via various visualization methods.
    """

    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    username = request.user.username

    try:
        struct =  structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "No structure uploaded. Please upload a structure.")
        return HttpResponseRedirect('/charmming/fileupload/')
#        return HttpResponse('No structure uploaded')
    try:
        ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "Please build a working structure first.")
        return HttpResponseRedirect('/charmming/buildstruct/')
#        return HttpResponse('No structure uploaded')
    cgws = None
    try:
        cgws = structure.models.CGWorkingStructure.objects.get(workingstructure_ptr=ws.id)
    except structure.models.CGWorkingStructure.MultipleObjectsReturned:
        messages.error(request, "More than one CGWorkingStructure is associated to this working structure. Please report this bug, then rebuild your structure.")
        return HttpResponseRedirect('/charmming/fileupload')
    except:
        pass

    if not ws.isBuilt == "t": #Redirect to calculation page
        messages.error(request, "Please perform a calculation on your structure before visualizing.")
        return HttpResponseRedirect("/charmming/energy")
    filelst = []

    # files from segment construction
    if not cgws: #CG models don't have segments.
        for wseg in ws.segments.all():
            filelst.append((wseg.builtCRD.replace('.crd','.pdb'), 'Segment %s' % wseg.name))


    # now get all tasks associated with the structure
    tasks = structure.models.Task.objects.filter(workstruct=ws,active='y')
    for task in tasks:
        # check for PDB files associated with the task
        wfiles = structure.models.WorkingFile.objects.filter(task=task,type='pdb')
        for wfile in wfiles:
            s = wfile.canonPath.split('/')[-1]
            filelst.append((s, wfile.description))

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/visualize.html', {'filelst': filelst,'lesson_ok':lesson_ok,'dd_ok':dd_ok})

"""
ChemDoodle is very problematic and incomplete. Even years of development won't fix the fact
that they're released under a paid license for full functionality and
we can't use that.
#Lets user view PDB through ChemDoodle
def chemdoodle(request,filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "No structure uploaded. Please upload a structure.")
        return HttpResponseRedirect('/charmming/fileupload/')
#        return HttpResponse('No structure')

    orgfname = struct.location + '/' + filename
    newfname = struct.location + '/tmp-' + filename
    woof = open(struct.location + "/pdbpickle.dat") #THis will get modified once chemdoodle actually works...remember the localpickle!
    pdbfile = cPickle.load(woof)
    helix_info = ""
    woof.close()
    try:
        for line in pdbfile.get_metaData()['helix']: #Print out the helices/sheets
            helix_info = helix_info + "HELIX" + (" "*((67-len(line))+4)) +  line.upper() + "\n"
    except KeyError:
        pass
    pdb = pychm.io.pdb.PDBFile(orgfname)
#    woof = open(struct.location + "/" + "pdbpickle.dat")
#    pdb = cPickle.load(woof)
    pdb[0].write(newfname,outformat="pdborg",ter=True,end=False)
    woof.write("\n")
    woof = open(newfname,"a+")
    woof.write(helix_info)
    woof.write("END\n")
    woof.close()
    mycontent = ''
    fp = open(newfname, 'r')
    for line in fp:
        line = line.strip() + '\\n' #Otherwise the JS part of this will get literal newlines in the template rather than \n characters
        if line.startswith('REMARK') or line.startswith('Remark') or line.startswith('remark'):
            continue
        mycontent += line
    fp.close()
    os.unlink(newfname)
    return render_to_response('html/chemdoodle.html', {'content': mycontent})
"""
#Lets user view PDB through jsmol - same functions and scripting
def jmol(request,filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "No structure uploaded. Please upload a structure.")
        return HttpResponseRedirect('/charmming/fileupload/')
#        return HttpResponse('No structure')
    cg_model = False
    try:
        ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "Please build a working structure first.")
        return HttpResponseRedirect('/charmming/fileupload/')
    try:
        cgws = structure.models.CGWorkingStructure.objects.get(workingstructure_ptr=ws.id)
        cg_model = True
    except Exception as ex:
        pass

    if len(filename) < 4:
        messages.error(request, "Bad filename. Please verify that the file you are trying to open exists before accessing this page.")
        return HttpResponseRedirect("/charmming/visualize/")
    filename = struct.location.replace(charmming_config.user_home,'') + '/' + filename
    filename_end = filename.rsplit('/',1)[1] #i.e., "a-pro-5.pdb", etc.
    lesson_ok, dd_ok = checkPermissions(request)
    super_user = request.user.is_superuser
    structname = filename_end.split('.',1)[0] #a-pro-5
    try:
        isProtein = not("-good-" in filename_end or "-bad-" in filename_end or filename_end[3] == "-")
    except:
        messages.error(request, "Bad filename. Please verify that the file you are trying to open exists before accessing this page.")
        return HttpResponseRedirect("/charmming/visualize/")
    #Assume everything that isn't a good/bad hetatm segment to be a protein for the purposes of rendering, that way we deal with mutations
    #Also assume anything that has 3 letters followed by a dash at the start to be a ligand file.
    return render_to_response('html/jmol.html', {'structname':structname, 'filepath': filename, 'segid':'NA', 'resid':'NA','lesson_ok':lesson_ok,'dd_ok':dd_ok, 'isProtein':isProtein, 'super_user':super_user, 'cg_model':cg_model})

#I haven't the foggiest idea where this function is called, but I assume it works fine with JSmol. ~VS
def jmolHL(request,filename,segid,resid):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "No structure uploaded. Please upload a structure.")
        return HttpResponseRedirect('/charmming/fileupload/')
#        return HttpResponse('No structure')
    filename = struct.location.replace(charmming_config.user_home,'') + '/' + filename

    return render_to_response('html/jmol.html', {'filepath': filename, 'segid': segid, 'resid': resid })

#Lets user view PDB through GLmol
def glmol(request,filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "No structure uploaded. Please upload a structure.")
        return HttpResponseRedirect('/charmming/fileupload/')
#        return HttpResponse('No structure')
    try:
        ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "Please build a working structure first.")
        return HttpResponseRedirect('/charmming/fileupload/')
#        return HttpResponse('No working structure selected.')
    cg_model = False
    try:
        cgws = structure.models.CGWorkingStructure.objects.get(workingstructure_ptr=ws.id)
        cg_model = True
    except structure.models.CGWorkingStructure.MultipleObjectsReturned:
        messages.error(request, "More than one CGWorkingStructure is associated to this working structure. Please report this bug, then rebuild your structure.")
        return HttpResponseRedirect('/charmming/fileupload')
    except: #Everything else (0, query error) is probably fine, so we keep going
        pass

    if len(filename) < 4:
        messages.error(request, "Bad filename. Please verify that the file you are trying to open exists before accessing this page.")
        return HttpResponseRedirect("/charmming/visualize/")

    filepath = struct.location.replace(charmming_config.user_home,'') #ws is always in same dir as struct...
    filename = filepath + "/" + filename
    filename_end = filename.rsplit('/',1)[1]

    helix_info = "" #Holds helix/sheet info for GLmol to read off the textbox
    #We make a check for localpickle so we don't load the wrong thing...
    if ws.localpickle: #ws is only called so that we can get its localpickle object.
        datafile = open(ws.localpickle)
    else:
        datafile = open(struct.pickle)
    pdbfile = cPickle.load(datafile) #Get the PDB file...
    datafile.close()
    try:
        for line in pdbfile.get_metaData()['helix']: #Print out the helices/sheets
            helix_info = helix_info + "HELIX" + (" "*((67-len(line))+4)) +  line.upper() + "\\n" 
    except KeyError:
        pass
    try:
        for line in pdbfile.get_metaData()['sheet']:
            if len(line) < 60: #i.e. if it's a "starting" record
                helix_info = helix_info + "SHEET" + (" "*((31-len(line))+4)) + line.upper() + "\\n"
            else:
                helix_info = helix_info + "SHEET" + (" "*((60-len(line))+4)) + line.upper() + "\\n"
    except KeyError:
        pass
    #There really isn't anything to be done if there's no helix/sheet info present because while it might mean that we need to run STRIDE it also might just mean there isn't any (like in 1yjp). Running STRIDE is cheap but it would probably be best to do it at the pychm/PDBFile I/O level rather than up here to guarantee everything is stored in the pickle's metadata
    #Since our PDBPickle doesn't know about the spacing and just stores the numbers, we need to "shift" the spaces to the left such that a 1-character sheet/helix number takes up 4 spaces after the HELIX/SHEET identifier, a 2-character number takes up 3 spaces, etc. 
    #The length of a "basic" HELIX line (1-character number) is 67, so we subtract from 67, then add 4 to get the correct number of spaces
    #Similarly the length of a basic SHEET line is 60.
    #Now we have the helix information for GLmol stored. 
    # WARNING! We always take model00. Does this ever matter? I have no idea!
    lesson_ok, dd_ok = checkPermissions(request)
    try:
        isProtein = not("-good-" in filename_end or "-bad-" in filename_end or filename_end[3] == "-")
        #TODO: Modify this! It's not a good way to tell and we really need a better way. Maybe attach an isProtein to the file objects...?
    except:
        messages.error(request, "Bad filename. Please verify that the file you are trying to open exists before accessing this page.")
        return HttpResponseRedirect("/charmming/visualize/")
    return render_to_response('html/glmol.html', {'cg_model': cg_model, 'filepath': filename, 'segid':'NA', 'resid':'NA','lesson_ok':lesson_ok,'dd_ok':dd_ok, 'isProtein':isProtein, 'helix_info':helix_info})


#Allows the user to see what processes their PDBs are undergoing
def viewstatus(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        file = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return render_to_response('html/statusreport.html', {'structure': None })

    try:
        workingStruct = structure.models.WorkingStructure.objects.filter(structure=file,selected='y')[0]
        haveWorkingStruct = True
    except:
        haveWorkingStruct = False
        workingStruct = None

    statuses = []
    if workingStruct:
        if not workingStruct.locked:
            workingStruct.updateActionStatus()
        t = structure.models.Task.objects.filter(workstruct=workingStruct,active='y')
        #TODO: Modify this for multi-task support. However we don't care until tasks don't just overwrite themselves.
        foo_tasks = set() #holds task actions so we don't repeat tasks
        for task in t:
            cankill = 0
            if task.status == 'R':
               cankill = 1
            if task.action not in foo_tasks:
                foo_tasks = foo_tasks.union(set([task.action]))
                statuses.append((task.action,task.get_status_display(),cankill,task.id))

    return render_to_response('html/statusreport.html', {'structure': file, 'ws':workingStruct, 'haveWorkingStruct': True, 'statuses': statuses})

def viewtaskoutput(request, taskid):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    tdict = {}
    try:
        tid = int(taskid)
        task = Task.objects.get(id=tid) #This is unique. Hopefully.
        if task.action == "dsfdocking":
            dockjob = dd_infrastructure.models.jobs.objects.get(job_scheduler_id=task.jobID)
        output_location = task.workstruct.structure.location
        os.chdir(output_location)
        out_action = task.action
    except Exception as ex:
        tdict['intro_string'] = "Task not found. Please report this issue." + str(ex) #Since we have a jQuery frame already, let's just not worry about making one
        return render_to_response('html/view_task_output.html',tdict)
    if task.action == "dsfdocking":
        #This is where things get weird since the file paths are completely different
        fp = open(charmming_config.user_home + "/dd/jobs/" + request.user.username + "/dd_job_" + str(dockjob.id) + "/run.out")
    elif task.action == "neutralization": #We need some extra checks here because this is a two-part task
        try:
            path_to_solv = task.workstruct.structure.location + "/" + task.workstruct.identifier + "-neutralization.out"
            os.stat(path_to_solv) #TODO: MOdify this for multi-names!
            #If the file exists, we read that, if not, we read solvation.out.
            fp = open(path_to_solv)
        except:
            fp = open(path_to_solv.replace("neutraliz","solv")) #Minimal form.
    elif task.action == "redox": #There's a lot of scripts going on and no real "main" one doing execution.
    #If you can figure out a decent system for this stuff, please tell me. ~VS
        tdict['intro_string'] = "Redox information cannot be displayed. Please use the progress tracker instead."
        return render_to_response('html/view_task_output.html',tdict)
    else:
        fp = open(task.workstruct.identifier + "-" + out_action + ".out") #e.g. 1yjp-energy.out
        #actions are standardized now. Whoever reads this in the future, make SURE your task ACTION matches the file name. 
        #TODO: Adapt this code for mutli-file task recognition.
    fcontents = fp.read()
    fp.close()
    top_string = "Structure: " + task.workstruct.structure.name + ", Working Structure: " + task.workstruct.identifier + ", Action: " + task.action.capitalize()
    tdict['file_contents'] = fcontents
    tdict['intro_string'] = top_string
#    except:
#        return HttpResponse("No such task.") #This works because it'd basically replace the status page with the error. Doing messages.error just messes things up here.
    return render_to_response('html/view_task_output.html', tdict)

# kill a task
def killTask(request,taskid):
#    logfp = open('/tmp/killTask.txt', 'w')
#    logfp.write('In killTask\n')

    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)

#    logfp.write('Got task ID = %s\n' % taskid)
    try:
        tid = int(taskid)
        t = Task.objects.get(id=tid)
    except:
#        logfp.write('No such task.\n')
#        logfp.close()
        return HttpResponse('No such task')

#    logfp.write('Found task, status = %s!\n' % t.status)
    if t.workstruct.structure.owner != request.user:
#        logfp.write('No permission\n')
#        logfp.close()
        return HttpResponse('No permission')
    if t.status != 'R':
#        logfp.write('Task not running!\n')
#        logfp.close()
        return HttpResponse('Job not running')

#    logfp.write('Calling task kill method.\n')
#    logfp.close()
    t.kill()
    return HttpResponse("OK")


#Used to report an error. Used in calling problem.html
def reportError(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)
    try:
        if(request.POST['errordescription'] != ''):
            emailmsg = 'Bug Reported!\n'
            emailmsg += 'User was: ' + request.user.username + '\n'
            emailmsg += 'Subject of error is: ' + request.POST['subject'] + '\n'
            emailmsg += 'Error details given by user are: \n' + request.POST['errordescription']
            mail_admins('Bug Reported',emailmsg,fail_silently=False)
            return HttpResponse("Bug was successfully reported. Thank you for reporting this bug.")
        ##TODO: Add non-error message. This may be fixed with Django versions higher than 1.2.
        else:
            messages.error(request, "You must type something into the form...")
            return HttpResponseRedirect("/charmming/about/")
#            return HttpResponse("You must type something into the form...\n")
    except:
        messages.error(request, "Bug failed to be sent properly. Please check all fields or E-mail btamiller@gmail.com")
        return HttpResponseRedirect("/charmming/about/")
#        return HttpResponse("Bug failed to be sent properly. Please check all fields or E-mail btamiller@gmail.com.")



#Lets user see the list of PDBs they have uploaded
#THIS IS NOT FOR VIEWING JMOL/CHEMAXON STUFF
def viewpdbs(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    user_pdbs = structure.models.Structure.objects.filter(owner=request.user)

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/viewpdbs.html', {'user_pdbs': user_pdbs, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

#This is for changing the currently selected PDB
#on the SELECT/EDIT PDBs page
#VS note: This appears to be a holdover from the oldcharmming times, but it still is used in the change PDBs page,
#as an AJAX request, wiuth a template purely to hold something in place.
#There should be a more efficient way of doign this, but this seems pretty much the bleeding edge of what django can do.
def switchpdbs(request,switch_id):
#    logfp = open('/tmp/switch.txt', 'w')
#    logfp.write('In switchpdbs\n')
#    logfp.close()
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        oldfile = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
        oldfile.selected = ''
        oldfile.save()
    except:
        pass
    newfile = structure.models.Structure.objects.filter(owner=request.user,name=switch_id)[0] 

    newfile.selected = 'y'
    newfile.save()

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/switchpdb.html',{'oldfile':oldfile,'newfile':newfile,'lesson_ok':lesson_ok,'dd_ok':dd_ok})

#This calculates the energy
def energyform(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)

    #chooses the file based on if it is selected or not
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please submit a structure.")
        return HttpResponseRedirect("/charmming/fileupload/")
#        return HttpResponse("Please submit a structure first.")
    try:
         ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        messages.error(request, "Please build a working structure before performing any calculations.")
        return HttpResponseRedirect("/charmming/buildstruct/")
#        return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before attempting an energy calculation")

    cg_model = False
    try:
        cgws = structure.models.CGWorkingStructure.objects.get(workingstructure_ptr=ws.id)
        cg_model = True
    except:
        pass

    tdict = {}
    tdict['ws_identifier'] = ws.identifier
    tdict['cg_model'] = cg_model

    os.chdir(struct.location)

    energy_lines = []
    try:
        energyfp = open('%s-energy.txt' % ws.identifier,'r')
        energy_lines = energyfp.readlines()
        energyfp.close()
    except:
        pass

    tdict['energy_lines'] = energy_lines
    sl = []
    sl.extend(structure.models.Segment.objects.filter(structure=struct,is_working='n'))
    sl = sorted(map(lambda x: x.name, sl))
    tdict['seg_list'] = sl

    if request.POST.has_key('form_submit'):
        try:
            oldtsk = energyTask.objects.filter(workstruct=ws,active='y')[0]
            oldtsk.active = 'n'
            oldtsk.save()
        except:
            pass

        et = energyTask()
        et.setup(ws)
        et.active = 'y'
        et.action = 'energy'
        et.save()
        if ws.isBuilt != 't':
            pTask = ws.build(et)
            pTaskID = pTask.id
        else:
            pTaskID = int(request.POST['ptask'])
            pTask = Task.objects.get(id=pTaskID)
        if request.POST.has_key('useqmmm'):
            saveAtomSelections(request, ws, pTask)
        return calcEnergy_tpl(request,ws,pTaskID,et)

    else:
        # get all workingFiles associated with this struct
        tasks = Task.objects.filter(workstruct=ws,status='C',active='y',modifies_coordinates=True)
        tdict['tasks'] = tasks

        tdict = getAtomSelections(tdict,ws) #Gets the atom selections by putting new keys in the dictionary, then returns it and puts it into tdict
        lesson_ok, dd_ok = checkPermissions(request)
        tdict['messages'] = get_messages(request)
        tdict['lesson_ok'] = lesson_ok
        tdict['dd_ok'] = dd_ok
        return render_to_response('html/energyform.html', tdict)

#Why do we need workstruct as an arg if eobj has it attached to itself?
#TODO: Optimize this method!
def calcEnergy_tpl(request,workstruct,pTaskID,eobj):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    logfp = open("/tmp/energy-fail.txt","w")
    logfp.write(str(eobj.modifies_coordinates))

    pTask = Task.objects.get(id=pTaskID)

    postdata = request.POST
    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'], template_dict['parameter_list'] = workstruct.getTopparList()
    template_dict['output_name'] = workstruct.identifier + '-ener'
    template_dict['input_file'] = workstruct.identifier + '-' + pTask.action
    if workstruct.topparStream:
        template_dict['tpstream'] = workstruct.topparStream.split()
    else:
        template_dict['tpstream'] = []


    eobj.finale = 0.0

    #If the user wants to solvate implicitly the scpism line is needed
    #84 will be the scpism number in this program
    solvate_implicitly = 0
    try:
        if(postdata['solvate_implicitly']):
            solvate_implicitly = 1
    except:
        pass
    template_dict['solvate_implicitly'] = solvate_implicitly

    # check to see if PBC needs to be used -- if so we have to set up Ewald
    if request.POST.has_key('usepbc'):
        if solvate_implicitly:
            return output.returnSubmission('energy', error='Invalid options')

        # decide if the structure we're dealing with has
        # been solvated.
        solvated = False
        wfc = pTask
        while True:
            if wfc.action == 'solvatom' or wfc.action == 'neutralization':
                solvated = True
                break
            if wfc.parent:
                wfc = wfc.parent
            else:
                break
        if not solvated:
            return output.returnSubmission('energy', error='Requested PBC on unsolvated structure')

        dopbc = True
        eobj.usepbc = 'y'
        try:
            sp = solvationTask.objects.get(workstruct=workstruct,active='y')
        except:
            return output.returnSubmission("Energy", error="Couldn't find solvation parameters")
        template_dict['xtl_x'] = sp.xtl_x
        template_dict['xtl_y'] = sp.xtl_y
        template_dict['xtl_z'] = sp.xtl_z
        template_dict['xtl_angles'] = "%10.6f %10.6f %10.6f" % (sp.angles[0],sp.angles[1],sp.angles[2])
        template_dict['xtl_ucell'] = sp.solvation_structure
        template_dict['ewalddim'] = sp.calcEwaldDim()

        if template_dict['xtl_ucell'] == 'sphere':
            messages.error(request, "Cannot do PBC on a sphere.")
            return HttpResponseRedirect("/charmming/solvation/")
#            return HttpResponse('Cannot do PBC on a sphere')
    else:
        dopbc = False
        eobj.usepbc = 'n'
    template_dict['dopbc'] = dopbc
    modelType = False


    if postdata.has_key('useqmmm'):
        eobj.useqmmm = 'y'
        #input.checkForMaliciousCode(postdata['qmsele'],postdata)
        #We pre-validate...
        if postdata['model_selected'] not in ["oniom","qmmm"]:
            messages.error(request,"Invalid model for Multi-Scale modeling.")
            return HttpResponseRedirect("/charmming/about/")
        else:
            modelType = postdata['model_selected']
        try: #This is obsolete for ONIOM...
            eobj.qmmmsel = postdata['qmsele']
        except:
            pass

        if dopbc:
            messages.error(request, "Cannot combine QM/MM and periodic boundary conditions.")
            return HttpResponseRedirect("/charmming/solvation/")
#            return HttpResponse('Cannot combine QM/MM and periodic boundary conditions')
    else:
        eobj.useqmmm = 'n'
    if eobj.useqmmm == 'y': #We branch to include ONIOM here.
        template_dict['modelType'] = modelType
        template_dict['useqmmm'] = True
        eobj.modelType = modelType
        if modelType == "qmmm":
            atomselection = selection.models.AtomSelection.objects.get(workstruct=eobj.workstruct) #TODO: Update for names. ALso don't worry about this yet - only one record should ever be returned anyway.
            qmparams = makeQchem_val("qmmm",atomselection)
            qmparams['jobtype'] = 'Force'
            template_dict = makeQChem_tpl(template_dict,qmparams,False,eobj)
        elif modelType == "oniom":
            try:
                template_dict = structure.mscale.make_mscale(template_dict, request, modelType, eobj)
            except:
                logfp.write('Failed to generate template\n')
                logfp.write(traceback.format_exc())
                logfp.flush()
                logfp.close()
    else:
        template_dict['useqmmm'] = False

    eobj.save()
    # lessons are still borked
    #if file.lesson_type:
    #    lessonaux.doLessonAct(file,"onEnergySubmit",request.POST)
    t = get_template('%s/mytemplates/input_scripts/calcEnergy_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))

    energy_input = workstruct.structure.location + "/" + workstruct.identifier + "-energy.inp"
    energy_output = workstruct.identifier + "-energy.out"
    inp_out = open(energy_input,'w')
    inp_out.write(charmm_inp)
    inp_out.close()
    eobj.scripts += ',%s' % energy_input
    eobj.modifies_coordinates = False #We set it here, because it doesn't stick in the finish method, and it doesn't stick in the method that calls this one
    #For some very strange and unknown reason, the field is modified when this method (the one we are in) is called.
    #We force-set it to false just in case.
    eobj.save()

    # go ahead and submit
    if eobj.useqmmm == 'y' and modelType == "oniom":
        eobj.start(mscale_job=True)
    else:
        eobj.start()

    # loop until the energy calculation is finished
    #TODO: Eliminate this? It seems a bit silly since we do .query() with updateActionStatus anyway...
    while True:
        sstring = eobj.query()
        if "complete" in sstring:
            break
        if "failed" in sstring:
            return HttpResponse('Energy calculation failed. Please check output from ' +pTask.action + '.')

    #YP lessons status update
    try:
        lnum=eobj.workstruct.structure.lesson_type
        lesson_obj = eval(lnum+'.models.'+lnum.capitalize()+'()')
    except:
        lesson_obj = None

    if lesson_obj:
        lessonaux.doLessonAct(eobj.workstruct.structure,"onEnergySubmit",eobj,"")
    #YP

    workstruct.save()

    # for now, I am not going to create a workstructure for the energy since it's effectively a
    # no-op; talk to Lee about whether this should be done. It probably should be, at one point
    # VS note: We should do this UNTIL we make build an actual part of the build working structure page. At that point
    # there's no reason to carry it around for any reason other than having it as a pTask, which shouldn't matter.
    return parseEnergy(workstruct,energy_output,eobj)

#This I believe is used for the Energy render AJAX request. It's weird and a holdover from oldcharmming.
def parseEnergy(workstruct,output_filename,enerobj=None):
    time.sleep(2)

    outfp = open(workstruct.structure.location + '/' + output_filename,'r')
    ener = re.compile('ENER ENR')
    dash = re.compile('----------')
    charmm = re.compile('CHARMM>')
    initiate = 0
    #the second extern marks the end of the energy output
    dash_occurance = 0
    energy_lines = ''

    for line in outfp:
        line = line.lstrip()
        if line.startswith('ENER>'):
            enerobj.finale = float(string.splitfields(line)[2])
            enerobj.save()
        if(dash_occurance == 2 or (initiate > 0 and line.strip().startswith('CHARMM>'))):
            if(line.startswith('CHARMM>')):
                break
            energy_lines += line
            break
        if(ener.search(line) or initiate > 0):
            if(initiate == 0): initiate = 1
            energy_lines += line
        if(dash.search(line) and initiate > 0):
            dash_occurance += 1
    outfp.close()
    writefp = open(workstruct.structure.location + '/' + workstruct.identifier + '-energy.txt','w')
    writefp.write(energy_lines)
    writefp.close()

    # lessons are still borked
    #if file.lesson_type:
    #    lessonaux.doLessonAct(file=file,function="onEnergyDone",finale=enerobj.finale)

    return render_to_response('html/displayenergy.html',{'linelist':energy_lines})

def getSegs(Molecule,Struct,auto_append):
    #auto_append is never used so I will use it for the custom ligand build...
    #True means normal/PDB.org build, False means custom build
    Struct.save()
    #Why do we save before doing anything?
    logfp = open('/tmp/getsegs.txt', 'w')
    logfp.write('In getSegs\n')
    sl = []
    sl.extend(structure.models.Segment.objects.filter(structure=Struct))
    for seg in Molecule.iter_seg():
        logfp.write('Found segment %s\n' % seg.segid)

        reslist = seg.iter_res()
        firstres = reslist.next().resName
        try:
            secondres = reslist.next()
            resName = "N/A"
        except StopIteration: #End of iterator, returns exception instead of null
            resName = firstres

        newSeg = structure.models.Segment()
        newSeg.structure = Struct
        newSeg.is_working = 'n'
        newSeg.name = seg.segid
        newSeg.type = seg.segType
        newSeg.is_custom = not(auto_append) #False if auto_append is True, and vice-versa.
        newSeg.resName = resName #This is for buildstruct page and other misc uses
        logfp.write("auto_append = " + str(newSeg.is_custom) + "\n")

        foo = map(lambda x: x.name,sl)
        bar = map(lambda x: x.structure,sl)
        logfp.write(str(foo) + "\n" + str(bar)) 
        #the following checks for whether the segment you're "getting" is already in the structure, or whether it's not. THis is used for custom ligands
        #such that you can add ligands to an already-existing structure without duplicating its segments
        #If the if test fails, the save is skipped - newSeg gets stuff associated to it but doesn't get committed to the database.
        if newSeg.name in map(lambda x: x.name,sl) and newSeg.structure in map(lambda x: x.structure,sl):
            logfp.write("Segment " + newSeg.name + " skipped\n")#If a segment with that name already exists in that structure, skip it and keep going.
            continue
        logfp.write('new seg object in charmming created\n')

        if seg.segType in ['pro']:
            newSeg.rtf_list = charmming_config.data_home + '/toppar/' + charmming_config.default_pro_top
            newSeg.prm_list = charmming_config.data_home + '/toppar/' + charmming_config.default_pro_prm
        elif seg.segType in ['rna','dna']:
            newSeg.rtf_list = charmming_config.data_home + '/toppar/' + charmming_config.default_na_top
            newSeg.prm_list = charmming_config.data_home + '/toppar/' + charmming_config.default_na_prm
        elif seg.segType == 'good':
            newSeg.rtf_list = ''
            newSeg.prm_list = ''
            newSeg.stream_list = charmming_config.data_home + '/toppar/toppar_water_ions.str'

        # get list of potential REDOX only segments. These segments will only
        # be used for REDOX calculations and otherwise aren't part of the final
        # constructed model. NB, we're going to assume that the first model is
        # canonical for these.
        if seg.segType == 'bad':
            ok = True
            for res in seg.iter_res():
                if res.resName not in ['fs4','sf4']:
                   ok = False
                   break
            if ok:
                newSeg.fes4 = True

        # set default patching type
        newSeg.set_default_patches(firstres)        
        newSeg.save()
#Should we save Struct again? We don't change anything here.
    logfp.write('All done\n')
    logfp.close()

def newupload(request, template="html/fileupload.html"):
    """
    Handle a new file being uploaded.
    """
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    input.checkRequestData(request)
    nmodels = 1
    makeGoModel = False
    makeBLNModel = False

    username = request.user.username
    u = User.objects.get(username=username)
    location = charmming_config.user_home + '/' + request.user.username + '/'

    all_messages = get_messages(request) #Gets any error messages from redirects.
    if request.POST.has_key('modeltype') and request.POST['modeltype'] == 'gomodel':
        makeGoModel = True
    if request.POST.has_key('modeltype') and request.POST['modeltype'] == 'blnmodel':
        makeBLNModel = True
    
    # check and see if this post is meant to be part of a lesson.
    #NOTE: Because Django is having issues with class inheritance, each file
    #Will have a lesson type (lesson1, lesson2, etc.) and a primary key id.
    #It's not a fun hack but until Django is changed it will have to be like this
    #YP
    #Upon Django upgrade we should see if this is possible 
    #VS
    if request.POST.has_key('lesson') and request.POST['lesson'] !=  'nolesson':
        lnum = request.POST['lesson']
        lessontype = lnum
        lesson_obj = eval(lnum+'.models.'+lnum.capitalize()+'()')

    else:
        lessonid = None
        lessontype = None
        lesson_obj = None

    if lesson_obj:
        lesson_obj.user = request.user
        lesson_obj.save()
        lessonid = lesson_obj.id
    #YP end of checking for lesson   

    try:
        request.FILES['pdbupload'].name
        file_uploaded = 1
    except:
        file_uploaded = 0

    # begin gigantic if test
    if request.POST.has_key('sequ') and request.POST['sequ']:
        struct = structure.models.Structure()
        struct.owner = request.user
        struct.title = 'Custom amino acid sequence'

        dname = 'sequ'
        tmpdname = dname
        version = -1
        while os.path.exists(location + '/' + dname):
            version += 1
            dname = tmpdname + "-" + str(version)

        struct.name = dname
        struct.location = location + dname
        os.mkdir(struct.location)
        os.chmod(struct.location, 0775)
        fullpath = struct.location + '/sequ.pdb'
        struct.putSeqOnDisk(request.POST['sequ'], fullpath)
        try:
            pdb = pychm.io.pdb.PDBFile(fullpath)
            mol = pdb[0]
            mol.parse()
            getSegs(mol,struct,auto_append=True)
        except ValueError:
            messages.error(request, "Error parsing PDB file.")
            return HttpResponseRedirect("/charmming/fileupload/")
        # store the pickle file containing the structure
        #TODO: Add dname-pdpickle.dat as filename. Do once the rest of the infrastructure is complete...
        pfname = location + dname + '/' + 'pdbpickle.dat'
        pickleFile = open(pfname,'w')
        cPickle.dump(pdb,pickleFile)
        pickleFile.close()
        struct.pickle = pfname

        # unselect the existing structure
        try:
            oldfile = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
            oldfile.selected = ''
            oldfile.save()
        except:
            pass

        struct.selected = 'y'
        #YP lessons
        if lesson_obj:
            postdata = request.POST
            struct.lesson_type=lnum
            struct.lesson_id=lesson_obj.id
            struct.save()
            lessonaux.doLessonAct(struct,"onFileUpload",postdata,"")
        #end YP lessons
        struct.save()
        return HttpResponseRedirect('/charmming/buildstruct/')

    elif file_uploaded or request.POST.has_key('pdbid'):
        if file_uploaded:
            filename = request.FILES['pdbupload'].name
            specialchars_string = "#$/;\n\\_+=[]{}()&^%"
            specialchars = set(specialchars_string)
            if len(specialchars.intersection(filename)) > 0:
                messages.error(request, "Your structure name cannot contain any of the following characters: " + specialchars_string)
                return HttpResponseRedirect("/charmming/about/")
        elif request.POST.has_key('pdbid'):
            filename = request.POST['pdbid'].strip().lower()
        filename = filename.lower()

        # figure out where this file is going to go
        tmpdname = filename.split('.')[0]
        dname = tmpdname
        version = -1
        while os.path.exists(location + '/' + dname):
            version += 1
            dname = tmpdname + "-" + str(version)

        os.mkdir(location + '/' + dname, 0775)
        os.chmod(location + '/' + dname, 0775)
        fullpath = location + '/' + dname + '/' + filename


        # Put the initial PDB onto the disk
        if file_uploaded:
            temp = open(fullpath, 'w')
            for fchunk in request.FILES['pdbupload'].chunks():
                temp.write(fchunk)
            temp.close()
        elif request.POST.has_key('pdbid'):
            pdbid = request.POST['pdbid']
            fullpath += '.pdb'
            try:
                conn = HTTPConnection("www.pdb.org")
                conn.request("GET", "/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % pdbid)
                resp = conn.getresponse()
                outfp = open(fullpath, 'w+')
                if resp.status != 200:
                        prob_string = "The PDB server returned error code %d." % resp.status
                        messages.error(request, prob_string)
                        return HttpResponseRedirect("/charmming/fileupload/")
#                        return render_to_response('html/problem.html',{'prob_string': prob_string})
                pdb_file = resp.read()
                if "does not exist" in pdb_file:
                    outfp.close()
                    messages.error(request, "The PDB servers report that this structure does not exist.")
                    return HttpResponseRedirect("/charmming/fileupload/")
#                    return render_to_response('html/problem.html',{'prob_string': 'The PDB server reports that the structure does not exist.'})

                outfp.write(pdb_file)
                outfp.close()
                conn.close()
            except:
                # print an error page...
                messages.error(request, "There was an error processing your file. Please contact the server administrator for more details.")
                return HttpResponseRedirect("/charmming/fileupload/")
#                return render_to_response('html/problem.html',{'prob_string': 'There was an exception processing the file -- contact the server administrator for more details.'})


        if file_uploaded and ( filename.endswith('crd') or filename.endswith('cor') ):
            # set up the new structure object
            struct = structure.models.Structure()
            struct.name = dname
            #YP
            struct.owner = request.user
            #YP
            struct.location = location + dname
            #YP
            #pdb = pychm.io.crd.CRDFile(fullpath)
            try:
                pdb = pychm.io.pdb.PDBFile()
                thisMol = pychm.io.pdb.get_molFromCRD(fullpath)
                pdb._mols["model00"] = thisMol
                getSegs(mol,struct,auto_append=True)
            except ValueError:
                messages.error(request, "Error parsing PDB file.")
                return HttpResponseRedirect("/charmming/fileupload/")
        else:
            try:
                pdb = pychm.io.pdb.PDBFile(fullpath)
                mnum = 1
                struct = structure.models.Structure()
                struct.location = location + dname
                struct.name = dname
                struct.owner = request.user
                struct.getHeader(pdb.header)

                for thisMol in pdb.iter_models():
                    thisMol.parse()
                    if mnum == 1: getSegs(thisMol,struct,auto_append=True) #Leave these as True.
                    mnum += 1
            except ValueError:
                messages.error(request, "Error parsing PDB file.")
                return HttpResponseRedirect("/charmming/fileupload/")

        # store the pickle file containing the structure
        #TODO: Add dname-pdbpickle.dat as filename. DO once infrastructure is complete.
        pfname = location + dname + '/' + 'pdbpickle.dat'
        pickleFile = open(pfname,'w')
        cPickle.dump(pdb,pickleFile)
        pickleFile.close()
        struct.pickle = pfname


        # unselect the existing structure
        try:
            oldfile = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
            oldfile.selected = ''
            oldfile.save()
        except:
            pass

        struct.selected = 'y'
        #YP lessons
        if lesson_obj:
            postdata = request.POST
            struct.lesson_type=lnum
            struct.lesson_id=lesson_obj.id
            struct.save()
            lessonaux.doLessonAct(struct,"onFileUpload",postdata,"")
        #end YP lessons

        # set this structure as selected
        
        struct.save()
        return HttpResponseRedirect('/charmming/buildstruct/')

    # end of ye gigantic if test

    form = structure.models.PDBFileForm()

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/fileuploadform.html', {'form': form, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok, 'messages':all_messages} )


# This function populates the form for building a structure
def buildstruct(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    input.checkRequestData(request)

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please submit a structure first.")
        return HttpResponseRedirect("/charmming/fileupload/")
#        return HttpResponse("Please submit a structure first.")

    # pick out a proposed name for this working structure
    proposedname = struct.name
    existingWorkStructs = structure.models.WorkingStructure.objects.filter(structure=struct)

    if existingWorkStructs:
        usedNameList = [ews.identifier for ews in existingWorkStructs]

        while proposedname in usedNameList:
           m = re.search("([^-\s]+)-ws([0-9]+)$", proposedname)
           if m:
               basename = m.group(1)
               numbah   = int(m.group(2))
               proposedname = basename + "-ws" + str(numbah+1)
           else:
               proposedname += "-ws1"

    tdict = {}
    tdict['proposedname'] = proposedname
    all_messages = get_messages(request) #Get all error messages for display 
    tdict['messages'] = all_messages
    ws = []
    if len(existingWorkStructs) > 0:
        tdict['haveworkingstruct'] = True
        ws.extend(existingWorkStructs)
        ws.sort()
    else:
        tdict['haveworkingstruct'] = False
    tdict['built_list'] = ws

    # get list of all of the models that we can use
    st_models = []
    fp = open(struct.pickle, 'r')
    pdb = cPickle.load(fp)
    fp.close()
    tdict['model_list'] = pdb.keys()
    sl = []
    sl.extend(structure.models.Segment.objects.filter(structure=struct,is_working='n'))
    sl = sorted(sl,key=lambda x: x.name)
    customCount = 0 #Checks each seg for "is_custom", and if customCount > 1, make a fuss
    for seg in sl:
        if seg.is_custom:
            customCount += 1
    #By having the database field resName in each Segment object we can save a lot of grief
    #Somehow this doesn't work that well in normalmodes...
    tdict['customCount'] = customCount
    tdict['seg_list'] = sl
    tdict['disulfide_list'] = struct.getDisulfideList()
    tdict['proto_list'] = []
    tdict['super_user'] = request.user.is_superuser
    tdict['filepath'] = "/charmming/pdbuploads/" + struct.location.replace(charmming_config.user_home,'') + "/" + struct.name + ".pdb" #This assumes we have a PDB file! Please be careful with this.
    for seg in sl:
        if seg.type == "pro": #Good and badhets are not supported and will have trouble, especially if part of solvation since solvation changes the composition.
            tdict['proto_list'].extend(seg.getProtonizableResidues(pickleFile=pdb))
            #Why do we even pass good/bad hets into this if they don't have protonizable residues?
    #The following procedure is very similar to the one in selection.views. 
    #See selection.views.selectrstructure lines 180-196 for more data.
    #Since we hav emore than protein chains, we need to figure out the chain terminators in a more elaborate way
    if tdict['proto_list']: #We should only do these graphics if we have a proto list
        segmentlist = [] #Holds the segnames
        chain_terminators = []
        for segment in pdb['model00'].iter_seg(): #More slow code!
        #Connectivity difference in models is not a concern because this is BEFORE we build anything.
            atom = segment.iter_res().next().iter_atom().next().atomNum0
            segmentlist.append(segment)
            chain_terminators.append(atom)
        segmentlist = sorted(segmentlist,key=lambda x:x.iter_res().next().iter_atom().next().atomNum0)
        segmentlist = map(lambda x:x.segid,segmentlist)
        chain_terminators = sorted(chain_terminators)
        tdict['chain_terminators'] = json.dumps(chain_terminators)
        tdict['segmentlist'] = json.dumps(segmentlist)
#        obconv = openbabel.OBConversion()
#        mol = openbabel.OBMol()
#        obconv.SetInAndOutFormats("pdb","pdb") #This is silly but it allows us to add hydrogens
#        obconv.ReadFile(mol, (struct.location + "/" + struct.name + ".pdb").encode("utf-8"))
#        mol.AddHydrogens(False,True,7.4)
#        obconv.WriteFile(mol, (struct.location + "/" + struct.name + "_with_hydrogens.pdb").encode("utf-8"))
    tdict['structname'] = struct.name
    tdict['lesson_ok'], tdict['dd_ok'] = checkPermissions(request)
    return render_to_response('html/buildstruct.html', tdict)


def modstruct(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    input.checkRequestData(request)

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please submit a structure.")
        return HttpResponseRedirect('/charmming/fileupload/')


    if not request.POST.has_key('wsidentifier'):
        messages.error(request, "You must give this working structure an identifier.")
        return HttpResponseRedirect("/charmming/buildstruct/")
    else:
        if len(request.POST['wsidentifier']) > 20:
            messages.error(request, "Working structure identifier too long. Please make sure your identifier is under 20 characters long.")
            return HttpResponseRedirect("/charmming/buildstruct/")
    #You can't steal cookies in 20 characters. But just in case...
    if not request.POST['wsidentifier']: #Aren't these two equivalent?
        messages.error(request, "You must give this working structure an identifier.")
        return HttpResponseRedirect("/charmming/buildstruct/")

    if not request.POST.has_key('buildtype'):
        messages.error(request, "No build type specified.")
        return HttpResponseRedirect("/charmming/buildstruct/")

    if request.POST['buildtype'] == 'aa':
        segs = structure.models.Segment.objects.filter(structure=struct,is_working='n')
        seglist = []
        for s in segs:
            ulkey = 'select_' + s.name
            if request.POST.has_key(ulkey) and request.POST[ulkey] == 'y':
               seglist.append(s.name)
        if len(seglist) < 1:
            messages.error(request, "You must choose at least one segment!")
            return HttpResponseRedirect("/charmming/buildstruct/")

        tpdict = {}
        for seg in seglist:
            if request.POST.has_key('toppar_' + seg):
                tpdict[seg] = request.POST['toppar_' + seg]

                allowedList = ['standard','upload','autogen','redox']
                if request.user.is_superuser:
                    allowedList.extend(['dogmans','match','genrtf','antechamber'])

                if tpdict[seg] not in allowedList:
                    tpdict[seg] = 'standard'

                if tpdict[seg] == 'upload':
                    try:
                        uptop = request.FILES['topology_' + seg]
                        uppar = request.FILES['parameter_' + seg]
                    except:
                        messages.error(request, "Topology/parameter files not uploaded.")
                        return HttpResponseRedirect("/charmming/buildstruct/")

                    topfp = open(struct.location + '/' + request.POST['wsidentifier'] + '-' + seg + '.rtf', 'w+')
                    prmfp = open(struct.location + '/' + request.POST['wsidentifier'] + '-' + seg + '.prm', 'w+')
                    for chunk in uptop.chunks(): topfp.write(chunk)
                    for chunk in uppar.chunks(): prmfp.write(chunk)
            else:
                tpdict[seg] = 'standard'

        new_ws = structure.models.WorkingStructure()
        # to do, make this not contain spaces
        new_ws.identifier = request.POST['wsidentifier']
        new_ws.modelName = request.POST['basemodel']
        new_ws.associate(struct,seglist,tpdict)

        # Figure out terminal patching
        for segobj in new_ws.segments.all():
            fpvar = segobj.name + '_firstpatch'
            lpvar = segobj.name + '_lastpatch'

            if request.POST.has_key(fpvar):
                segobj.patch_first = request.POST[fpvar]
                segobj.save()
            if request.POST.has_key(lpvar):
                segobj.patch_last = request.POST[lpvar]
                segobj.save()

        # Figure out protonation and disulfide patching
        for pkey in request.POST.keys():
            if pkey.startswith('protostate_'):
                (junk,segid,resid) = pkey.split('_')

                if not segid in seglist: continue # not in the segments selected
                if request.POST[pkey] in ['hsd','lys','glu','asp']: continue # these are defaults, no need for a patch
                p = structure.models.Patch()
                p.structure = new_ws
                p.patch_segid = structure.models.Segment.objects.get(structure=struct,is_working='n',name=segid)
                p.patch_name = request.POST[pkey]
                p.patch_segres = "%s %s" % (segid,resid)
                p.save()

            if pkey.startswith('disul_'):
                (junk,segid1,resid1,segid2,resid2) = pkey.split('_')
                if not (segid1 in seglist and segid2 in seglist): continue # patch not valid for segments selected
                p = structure.models.Patch()
                p.structure = new_ws
                p.patch_name = 'disu'
                p.patch_segres = "%s %s %s %s" % (segid1,resid1,segid2,resid2)
                p.save()

    elif request.POST['buildtype'] == 'go':

        logfp = open('/tmp/build_go_model.txt', 'w')
        logfp.write('Build Go model.\n')
        logfp.close()

        seglist = []

        segs = structure.models.Segment.objects.filter(structure=struct,is_working='n')
        for s in segs:
            ulkey = 'go_select_' + s.name
            if request.POST.has_key(ulkey) and request.POST[ulkey] == 'y':
                seglist.append(s.name)

        new_ws = structure.models.CGWorkingStructure()
        new_ws.identifier = request.POST['wsidentifier']
        new_ws.modelName = request.POST['basemodel']
        new_ws.cg_type = 'go'
        new_ws.associate(struct,seglist,contactSet=request.POST['gm_contact_type'], nScale=request.POST['gm_nscale'], \
                         kBond=request.POST['gm_kbond'], kAngle=request.POST['gm_kangle'], contactrad=request.POST['gm_contactrad'])
        new_ws.save()

    elif request.POST['buildtype'] == 'bln':
        seglist = []

        segs = structure.models.Segment.objects.filter(structure=struct,is_working='n')
        for s in segs:
            ulkey = 'bln_select_' + s.name
            if request.POST.has_key(ulkey) and request.POST[ulkey] == 'y':
                seglist.append(s.name)

        new_ws = structure.models.CGWorkingStructure()
        new_ws.identifier = request.POST['wsidentifier']
        new_ws.modelName = request.POST['basemodel']
        new_ws.cg_type = 'bln'
        new_ws.associate(struct,seglist,kBondHelix=request.POST['bln_kbondhelix'],\
                         kAngleHelix=request.POST['bln_kanglehelix'],kBondSheet=request.POST['bln_kbondsheet'],kAngleSheet=request.POST['bln_kanglesheet'], \
                         kBondCoil=request.POST['bln_kbondcoil'],kAngleCoil=request.POST['bln_kanglecoil'])
        new_ws.save()
    else:
        messages.error(request, "Bad buildtype specified!")
        return HttpResponseRedirect("/charmming/buildstruct/")

    # ToDo: Figure out restraints

    # switch to using the new work structure
    try:
        old_ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        pass
    else:
        old_ws.selected = 'n'
        old_ws.save()

    new_ws.selected = 'y'
    new_ws.save()

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/built.html', {'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

def swap(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    input.checkRequestData(request)
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please select a structure.")
        return HttpResponseRedirect("/charmming/buildstruct/")
    if not request.POST.has_key('choosestruct'):
        messages.error(request, "Invalid structure choice.")
        return HttpResponseRedirect("/charmming/buildstruct/")

    try:
        old_ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        pass
    else:
        old_ws.selected = 'n'
        old_ws.save()

    try:
        new_ws = structure.models.WorkingStructure.objects.filter(structure=struct,identifier=request.POST['choosestruct'])[0]
    except:
        messages.error(request, "Invalid structure selected.")
        return HttpResponseRedirect("/charmming/buildstruct/")
    else:
        new_ws.selected = 'y'
        new_ws.save()

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/swapped.html', {'wsname': new_ws.identifier, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

#Why is this here?
def protonate(file):
    return render_to_response('html/protonate.html')
