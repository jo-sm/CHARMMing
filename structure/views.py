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
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from minimization.models import minimizeTask
from dynamics.models import mdTask, ldTask, sgldTask
from solvation.models import solvationTask
from account.models import *
from normalmodes.views import combineNmaPDBsForMovie
from normalmodes.aux import getNormalModeMovieNum
from normalmodes.models import nmodeTask
from apbs.models import redoxTask
from structure.qmmm import makeQChem, makeQChem_tpl, handleLinkAtoms, writeQMheader
from structure.models import Task, energyTask
from django.contrib.auth.models import User
from django.core.mail import mail_admins
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from account.views import isUserTrustworthy
from structure.aux import checkNterPatch
from lessons.models import LessonProblem
import output, lesson1, lesson2, lesson3, lesson4, lessonaux
import structure.models, input
import re
# import all created lessons by importing lesson_config
# also there is a dictionary called 'file_type' in lesson_config.py specififying the file type of files uploaded by the lessons  
from lesson_config import *
import os, sys, re, copy, datetime, time, stat
import mimetypes, string, random, glob, traceback, commands
import pychm.io, charmming_config, minimization
import cPickle

# problem during upload
def uploadError(request,problem):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    prob_string = "An unknown error occurred during file upload."
    if problem == 'overflow':
       prob_string = "The file you're trying to upload has more than 50,000 atoms."
    elif problem == 'parse_error':
       prob_string = "There was a problem trying to parse your structure file."
    return render_to_response('html/problem.html',{'prob_string': prob_string})

# Deletes user specified file from database and files relating to it
def deleteFile(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    #make sure request data isn't malicious
    input.checkRequestData(request)
    if request.POST.has_key('filename'):
        delete_filename = request.POST['filename']
    else:
        return HttpResponse('Bad')

    deleteAll = 0
    remove_extension = re.compile('\.\S*')
    #If user wants to remove all files, visit /deletefile/all_files/

    if(delete_filename == "all_files"):
        deleteAll = 1
        total_files = structure.models.Structure.objects.filter(owner=request.user)
    else:
        file = structure.models.Structure.objects.filter(owner=request.user,name=delete_filename)[0]
        total_files = [file]

    for s in total_files[::-1]:
        os.chdir(charmming_config.user_home + '/' + request.user.username)
        os.system("rm -rf " + s.name)

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

        s.delete()

    return HttpResponse('Done')


#Lets user view file in browser, puts it in iframe
def viewProcessContainer(request,filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    return render_to_response('html/viewprocessfilescontainer.html', {'filename':filename })

#Lets user view file in browser, puts it in iframe
def downloadFilesContainer(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    return render_to_response('html/downloadfilescontainer.html')

#Lets user view file in browser, puts it in iframe
def viewPDBsContainer(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    return render_to_response('html/viewpdbscontainer.html')

#Lets user view file in browser
def viewFiles(request,filename, mimetype = None):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    username = request.user.username

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return render_to_response('html/nopdbuploaded.html')
    try:
       os.stat("%s/%s" % (struct.location,filename))
    except:
       return HttpResponse("That file doesn't seem to exist. Maybe you deleted it?")

    mimetype = "Content-Type: text/richtext"
    response = HttpResponse(mimetype=mimetype)
    response.write(file("%s/%s" % (struct.location,filename), "rb").read())
    return response

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
        return HttpResponse("Please perform some operations on your structure")

    wsname = ws.identifier
    filelst = []
    headers = []

    # files from segment construction
    headers.append('Segment setup')
    for wseg in ws.segments.all():
        filelst.append(('Segment setup', wseg.builtPSF, 'PSF of segment %s' % wseg.name))
        filelst.append(('Segment setup', wseg.builtPSF.replace('.psf','.pdb'), 'PDB of segment %s' % wseg.name))
        filelst.append(('Segment setup', wseg.builtCRD, 'CRD of segment %s' % wseg.name))

        # there would be an input and output too
        filelst.append(('Segment setup', 'build-%s.inp' % wseg.name, 'Input file to build segment %s' % wseg.name))
        filelst.append(('Segment setup', 'build-%s.out' % wseg.name, 'Output file from building segment %s' % wseg.name))


    # now get all Tasks associated with the structure and get their workingFiles
    # hack alert: we only store CRDs now but we can extrapolate the PSFs, INPs, and OUTs
    # ToDo: store this info in the workingfile
    tasks = structure.models.Task.objects.filter(workstruct=ws)
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
    structure =  Structure.objects.filter(owner=request.user,selected='y')[0]
    tar_filename = structure.name + '.tar.gz'

    username = request.user.username
    os.chdir(charmming_config.user_home + username)

    try:
        os.unlink(tar_filename)
    except:
        pass

    # remove condor logs and error files
    os.system('rm -f ' + struct.name + '/*.{err,log}')

    os.system('cp ' + charmming_config.data_home + '/solvation/water.crd ' + struct.name)
    os.system('cp ' + charmming_config.data_home + '/calcewald.pl ' + struct.name)
    os.system('cp ' + charmming_config.data_home + '/savegv.py ' + struct.name)
    os.system('cp ' + charmming_config.data_home + '/savechrg.py ' + struct.name)
    os.system('cp ' + charmming_config.data_home + '/toppar/top_all27_prot_na.rtf ' + charmming_config.data_home + '/toppar/par_all27_prot_na.prm ' + struct.name)
    os.system('tar -czf ' + tar_filename + ' ' + struct.name)
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
        return HttpResponse("Hmmm ... that file doesn't seem to exist. Please file a bug report.")

    if mimetype is None:
        mimetype,encoding = mimetypes.guess_type(path + "/" + filename)
    response = HttpResponse(mimetype=mimetype)
    response['Content-Disposition'] = 'attachment; filename=%s' %filename
    response.write(file(path + "/" + filename, "rb").read())
    return response

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
        return HttpResponse('No structure uploaded')
    try:
        sval = os.stat("%s/%s" % (struct.location,filename))
    except:
        return HttpResponse('Oops ... that file no longer exists.')
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
        return HttpResponse("Please select a structure first.")
    try:
        ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        return HttpResponse('No structure uploaded')

    header_list = []
    file_list   = []
    tasks = structure.models.Task.objects.filter(workstruct=ws)
    for task in tasks:
        header_list.append(task.action)

        # get all input and output files associated with the task
        wfiles = structure.models.WorkingFile.objects.filter( Q(task=task), \
                   Q(type='inp') | Q(type='out'))

        for wf in wfiles:
            bn = wf.canonPath.split('/')[-1]
            ds = wf.description
            file_list.append((task.action,bn,ds))

    return render_to_response('html/viewprocessfiles.html', {'headers': header_list, 'files': file_list})

def visualize(request,filename):
    """
    Allows the user to visualize their structure via Jmol
    """

    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    username = request.user.username

    try:
        struct =  structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse('No structure uploaded')
    try:
        ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        return HttpResponse('No structure uploaded')

    filelst = []

    # files from segment construction
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

    return render_to_response('html/visualize.html', {'filelst': filelst})


#Let's user view PDB through ChemDoodle
def chemdoodle(request,filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse('No structure')
    orgfname = struct.location + '/' + filename
    newfname = struct.location + '/tmp-' + filename
    pdb = pychm.io.pdb.PDBFile(orgfname)
    pdb[0].write(newfname,outformat="pdborg",ter=True,end=True)

    mycontent = ''
    fp = open(newfname, 'r')
    for line in fp:
        line = line.strip() + '\\n'
        if line.startswith('REMARK') or line.startswith('Remark') or line.startswith('remark'): 
            continue
        mycontent += line
    fp.close()

    os.unlink(newfname)
    return render_to_response('html/chemdoodle.html', {'content': mycontent})


#Let's user view PDB through jmol
def jmol(request,filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse('No structure')
    filename = struct.location.replace(charmming_config.user_home,'') + '/' + filename

    return render_to_response('html/jmol.html', {'filepath': filename, 'segid':'NA', 'resid':'NA'})


def jmolHL(request,filename,segid,resid):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse('No structure')
    filename = struct.location.replace(charmming_config.user_home,'') + '/' + filename

    return render_to_response('html/jmol.html', {'filepath': filename, 'segid': segid, 'resid': resid })

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
    logfp = open('/tmp/viewstatus.txt', 'w')
    logfp.write('in viewstatus\n')
    if workingStruct:
        logfp.write('in ifworkingstruct\n')
        workingStruct.updateActionStatus()
        t = structure.models.Task.objects.filter(workstruct=workingStruct,active='y')
        for task in t:
            cankill = 0
            if task.status == 'R':
               cankill = 1
            statuses.append((task.action,task.get_status_display(),cankill,task.id))

    logfp.write('status = %s\n' % statuses)
    logfp.close()
    return render_to_response('html/statusreport.html', {'structure': file, 'haveWorkingStruct': True, 'statuses': statuses})

# kill a task
def killTask(request,taskid):
    logfp = open('/tmp/killTask.txt', 'w')
    logfp.write('In killTask\n')

    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)

    logfp.write('Got task ID = %s\n' % taskid)
    try:
        tid = int(taskid)
        t = Task.objects.get(id=tid)
    except:
        logfp.write('No such task.\n')
        logfp.close()
        return HttpResponse('No such task')

    logfp.write('Found task, status = %s!\n' % t.status)
    if t.workstruct.structure.owner != request.user:
        logfp.write('No permission\n')
        logfp.close()
        return HttpResponse('No permission')
    if t.status != 'R':
        logfp.write('Task not running!\n')
        logfp.close()
        return HttpResponse('Job not running')

    logfp.write('Calling task kill method.\n')
    logfp.close()
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
        else:
            return HttpResponse("You must type something into the form...\n")
    except:
        return HttpResponse("Bug failed to be sent properly. Please check all fields or E-mail btamiller@gmail.com.")



def fileuploadform(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    file = structure.models.PDBFileForm()
    lessonuser = 0
    groupList = request.user.groups.all()
    if request.user.groups.filter(name='lesson') or request.user.is_superuser:
        lessonuser = 1
    return render_to_response('html/fileuploadform.html', {'form':file,'lessonuser': lessonuser})

#Lets user see the list of PDBs they have uploaded
#THIS IS NOT FOR VIEWING JMOL/CHEMAXON STUFF
def viewpdbs(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    user_pdbs = structure.models.Structure.objects.filter(owner=request.user)
    return render_to_response('html/viewpdbs.html', {'user_pdbs': user_pdbs})

#This is for changing the currently selected PDB
#on the SELECT/EDIT PDBs page
def switchpdbs(request,switch_id):
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
    return render_to_response('html/switchpdb.html',{'oldfile':oldfile,'newfile':newfile})

#This calculates the energy
def energyform(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkRequestData(request)

    #chooses the file based on if it is selected or not
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    try:
         ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before attempting an energy calculation")

    os.chdir(struct.location)

    energy_lines = []
    try:
        energyfp = open('%s-energy.txt' % ws.identifier,'r')
        energy_lines = energyfp.readlines()
        energyfp.close()
    except:
        pass

    if request.POST.has_key('form_submit'):
        try:
            oldtsk = energyTask.objects.filter(workstruct=workstruct,active='y')[0]
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

        return calcEnergy_tpl(request,ws,pTaskID,et)

    else:
        # get all workingFiles associated with this struct
        tasks = Task.objects.filter(workstruct=ws,status='C',active='y')
        return render_to_response('html/energyform.html', {'ws_identifier': ws.identifier,'tasks': tasks, 'energy_lines': energy_lines})


def calcEnergy_tpl(request,workstruct,pTaskID,eobj): 
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    pTask = Task.objects.get(id=pTaskID)

    postdata = request.POST
    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'], template_dict['parameter_list'] = workstruct.getTopparList()
    template_dict['output_name'] = workstruct.identifier + '-ener'
    template_dict['input_file'] = workstruct.identifier + '-' + pTask.action

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
    dopbc = False
    if request.POST.has_key('usepbc'):
        if solvate_implicitly:
            return HttpResponse('Invalid options')

        # decide if the structure we're dealing with has
        # been solvated.
        solvated = False
        wfc = pstruct
        while True:
            if wfc.parentAction == 'solv':
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

        if template_dict['xtl_ucell'] == 'sphere':
            return HttpResponse('Cannot do PBC on a sphere')

    template_dict['dopbc'] = dopbc

    # lessons are still borked
    #if file.lesson_type:
    #    lessonaux.doLessonAct(file,"onEnergySubmit",request.POST)

    t = get_template('%s/mytemplates/input_scripts/calcEnergy_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))

    energy_input = workstruct.structure.location + "/" + workstruct.identifier + "-ener.inp"
    energy_output = workstruct.identifier + "-ener.out"
    inp_out = open(energy_input,'w')
    inp_out.write(charmm_inp)
    inp_out.close()
    eobj.scripts += ',%s' % energy_input    
    eobj.save()

    # go ahead and submit
    eobj.start()

    # loop until the energy calculation is finished
    while True:
        sstring = eobj.query()
        if "complete" in sstring:
            break
        if "failed" in sstring:
            return HttpResponse('Energy calculation failed. Please check output from appending the structure.')

    workstruct.save()

    # for now, I am not going to create a workstructure for the energy since it's effectively a
    # no-op; talk to Lee about whether this should be done. It probably should be, at one point

    return parseEnergy(workstruct,energy_output,eobj)

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

def getSegs(Molecule,Struct,auto_append=False):
    Struct.save()
    for seg in Molecule.iter_seg():
        reslist = seg.iter_res()
        firstres = reslist.next().resName

        newSeg = structure.models.Segment()
        newSeg.structure = Struct
        newSeg.is_working = 'n'
        newSeg.name = seg.segid
        newSeg.type = seg.segType

        if seg.segType in ['pro','rna','dna','good']:
            newSeg.rtf_list = charmming_config.data_home + '/toppar/top_all27_prot_na.rtf'
            newSeg.prm_list = charmming_config.data_home + '/toppar/par_all27_prot_na.prm'

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

    if request.POST.has_key('modeltype') and request.POST['modeltype'] == 'gomodel':
        makeGoModel = True
    if request.POST.has_key('modeltype') and request.POST['modeltype'] == 'blnmodel':
        makeBLNModel = True
    
    # check and see if this post is meant to be part of a lesson.
    #NOTE: Because Django is having issues with class inheritance, each file
    #Will have a lesson type (lesson1, lesson2, etc.) and a primary key id.
    #It's not a fun hack but until Django is changed it will have to be like this
    #YP
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

        pdb = pychm.io.pdb.PDBFile(fullpath)
        mol = pdb[0]
        mol.parse()
        getSegs(mol,struct)

        # store the pickle file containing the structure
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
        return HttpResponseRedirect('/charmming/buildstruct')

    elif file_uploaded or request.POST.has_key('pdbid'):
        if file_uploaded:
            filename = request.FILES['pdbupload'].name
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
                        return render_to_response('html/problem.html',{'prob_string': prob_string})
                pdb_file = resp.read()
                if "does not exist" in pdb_file:
                    outfp.close()
                    return render_to_response('html/problem.html',{'prob_string': 'The PDB server reports that the structure does not exist.'})

                outfp.write(pdb_file)
                outfp.close()
                conn.close()
            except:
                # print an error page...
                return render_to_response('html/problem.html',{'prob_string': 'There was an exception processing the file -- contact the server administrator for more details.'})


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
            pdb = pychm.io.pdb.PDBFile()
            thisMol = pychm.io.pdb.get_molFromCRD(fullpath)
            pdb._mols["model00"] = thisMol
            #YP
            getSegs(thisMol,struct,auto_append=True)
        
        else:
            pdb = pychm.io.pdb.PDBFile(fullpath)
            mnum = 1
            struct = structure.models.Structure()
            struct.location = location + dname
            struct.name = dname
            struct.owner = request.user
            struct.getHeader(pdb.header)

            for thisMol in pdb.iter_models():
                thisMol.parse()
                if mnum == 1: getSegs(thisMol,struct)
                mnum += 1

        # store the pickle file containing the structure
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
    return render_to_response('html/fileupload.html', {'form': form} )


# This function populates the form for building a structure
def buildstruct(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    input.checkRequestData(request)

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")

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
    sl.sort()
    tdict['seg_list'] = sl
    tdict['disulfide_list'] = struct.getDisulfideList()
    tdict['proto_list'] = []
    for seg in sl:
        tdict['proto_list'].extend(seg.getProtonizableResidues())

    return render_to_response('html/buildstruct.html', tdict)


def modstruct(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    input.checkRequestData(request)

    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure")


    if not request.POST.has_key('wsidentifier'):
        return HttpResponse("You must give this working structure an identifier")
    if not request.POST['wsidentifier']:
        return HttpResponse("You must give this working structure an identifier")

    if not request.POST.has_key('buildtype'):
        return HttpRespone("No build type specified")

    if request.POST['buildtype'] == 'aa':
        segs = structure.models.Segment.objects.filter(structure=struct,is_working='n')
        seglist = []
        for s in segs:
            ulkey = 'select_' + s.name
            if request.POST.has_key(ulkey) and request.POST[ulkey] == 'y':
               seglist.append(s.name)
        if len(seglist) < 1:
            return HttpResponse("You must choose at least one segment!")
    
        tpdict = {}
        for seg in seglist:
            if request.POST.has_key('toppar_' + seg):
                tpdict[seg] = request.POST['toppar_' + seg]
                if tpdict[seg] not in ['standard','upload','autogen','redox']:
                    tpdict[seg] = 'standard'

                if tpdict[seg] == 'upload':
                    try:
                        uptop = request.FILES['topology_' + seg]
                        uppar = request.FILES['parameter_' + seg]
                    except:
                        return HttpResponse("Topology/parameter files not uploaded")

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
                p.patch_name = 'disul'
                p.patch_segres = "%s %s %s %s" % (segid1,resid1,segid2,resid2)
                p.save()

    elif request.POST['buildtype'] == 'go':
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
                         domainScale=request.POST['gm_domainscale'], kBond=request.POST['gm_kbond'], kAngle=request.POST['gm_kangle'])
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
        new_ws.associate(struct,seglist,domainScale=request.POST['bln_domainscale'],kBondHelix=request.POST['bln_kbondhelix'],\
                         kAngleHelix=request.POST['bln_kanglehelix'],kBondSheet=request.POST['bln_kbondsheet'],kAngleSheet=request.POST['bln_kanglesheet'], \
                         kBondCoil=request.POST['bln_kbondcoil'],kAngleCoil=request.POST['bln_kanglecoil'])
        new_ws.save()
    else:
        return HttpResponse("Bad builttype specified!")

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

    return render_to_response('html/built.html')

def swap(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    input.checkRequestData(request)
    try:
        struct = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please select a structure")
    if not request.POST.has_key('choosestruct'):
        return HttpResponse("Invalid choice")

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
        return HttpResponse("Invalid structure selected")
    else:
        new_ws.selected = 'y'
        new_ws.save()

    return HttpResponse('Swapped')


def protonate(file):
    return render_to_response('html/protonate.html')    

def editstructurecontainer(request):
    return render_to_response('html/editstructurecontainer.html')
