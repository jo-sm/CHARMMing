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
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from minimization.models import minimizeParams
from minimization.views import append_tpl
from dynamics.models import mdParams, ldParams, sgldParams
from solvation.models import solvationParams
from account.models import *
from dynamics.views import combinePDBsForMovie
from normalmodes.views import combineNmaPDBsForMovie
from normalmodes.aux import getNormalModeMovieNum
from normalmodes.models import nmodeParams
from apbs.models import redoxParams
from structure.qmmm import makeQChem, makeQChem_tpl, handleLinkAtoms, writeQMheader
from django.contrib.auth.models import User
from django.core.mail import mail_admins
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from account.views import isUserTrustworthy
from structure.editscripts import generateHTMLScriptEdit
from structure.aux import checkNterPatch
import output, lesson1, lesson2, lesson3, lesson4, lessonaux
import structure.models, input
# import all created lessons by importing lesson_config
# also there is a dictionary called 'file_type' in lesson_config.py specififying the file type of files uploaded by the lessons  
from lesson_config import *
import os, sys, re, copy, datetime, time, stat
import mimetypes, string, random, glob, traceback, commands
import charmming.io, charmming_config
import cPickle

def getJobTime(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    file = structure.models.Structure.objects.filter(owner=request.user)[0]
    file.checkRequestData(request)
    os.chdir(file.location)
    if request.POST['time_filenames'] =="":
        return HttpResponse("")
    return HttpResponse(file.calcJobTime(request.POST['time_filenames'],request.POST['nstep']))

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
    structure.models.Structure.checkRequestData(request)
    delete_filename = request.POST['filename']
    deleteAll = 0
    remove_extension = re.compile('\.\S*')
    #If user wants to remove all files, visit /deletefile/all_files/

    if(delete_filename == "all_files"):
        deleteAll = 1
        total_files = Structure.objects.filter(owner=request.user)
    else:
       file =  Structure.objects.filter(owner=request.user,filename=delete_filename)[0]
       total_files = [file]
    for eachfile in total_files[::-1]:
        file = eachfile
        file_list = file.getFileList()
        input_list = file.getInputList()
        for item in file_list:
            try:
                os.remove(file.location + item)
                #Must also delete crd/psf files
                temp_item = remove_extension.sub('',item)
                psf = temp_item + ".psf"
                crd = temp_item + ".crd"
                os.remove(file.location + psf)
                os.remove(file.location + crd)
            except:
                pass
        for item in input_list:
            try:
                os.remove(file.location + item)
                #The submitJob routine creates a log and error file
                #so they must also be removed
                temp_item = remove_extension.sub('',item)
                err = temp_item + ".err"
                log = temp_item + ".log"
                os.remove(file.location + err)
                os.remove(file.location + log)
            except:
                pass
        os.chdir(file.location)
        os.system("rm -rf " + file.location  + "*" + file.stripDotPDB(file.filename) + "*")
        os.system("rm -f ANTECHAMBER* ATOMTYPE.INF")

        # clean up other models that are left behind.
        try:
            minimizations = minimizeParams.objects.filter(pdb=file,user=file.owner)
            for mnm in minimizations:
                mnm.delete()
        except:
            pass
        try:
            solvs = solvationParams.objects.filter(pdb=file)
            for mnm in solvs:
                mnm.delete()
        except:
            pass
        try:
            mds = mdParams.objects.filter(pdb=file)
            for mnm in mds:
                mnm.delete()
        except:
            pass
        try:
            lds = ldParams.objects.filter(pdb=file)
            for mnm in lds:
                mnm.delete()
        except:
            pass
        try:
            sglds = sgldParams.objects.filter(pdb=file)
            for mnm in sglds:
                mnm.delete()
        except:
            pass
        try:
            nmas = nmodeParams.objects.filter(pdb=file)
            for mnm in nmas:
                mnm.delete()
        except:
            pass
        try:
            eobs = energyParams.objects.filter(pdb=file)
            for eo in eobs:
                eo.delete()
        except:
            pass

        try:
            redoxparms = redoxParams.objects.filter(pdb=file)
            for ro in redoxparms:
                ro.delete()
        except:
            pass

        if file.selected == 'y' and deleteAll == 0:
            file.selected = ''
            allFileList = Structure.objects.filter(owner=request.user)
            if len(allFileList) >= 2:
                # We need to change the selection to another file, arbitrarilty we choose
                # the first one.
                for candidateSwitch in allFileList:
                     if candidateSwitch.filename != file.filename:
                         candidateSwitch.selected = 'y'
                         candidateSwitch.save()
                         break
        file.delete()
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
       os.stat("%s/%s/%s" % (charmming_config.user_home,username,filename))
    except:
       return HttpResponse("That file doesn't seem to exist. Maybe you deleted it?")

    mimetype = "Content-Type: text/richtext"
    response = HttpResponse(mimetype=mimetype)
    response.write(file("%s/%s/%s" % (charmming_config.user_home,username,filename), "rb").read())
    return response

#Wraps all files associated with the PDB into a tar and deletes it
def downloadFilesPage(request,mimetype=None):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        file =  Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return render_to_response('html/nopdbuploaded.html')
    os.chdir(file.location)
    filename_list = file.getDownloadSegmentList()
    appendfiles = []
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-final.pdb')
        appendfiles.append('new_' + file.stripDotPDB(file.filename) + '-final.pdb')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-final.psf')
        appendfiles.append('new_' + file.stripDotPDB(file.filename) + '-final.psf')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-final.crd')
        appendfiles.append('new_' + file.stripDotPDB(file.filename) + '-final.crd')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-append.inp')
        appendfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-append.inp')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-append.out')
        appendfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-append.out')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-rmsd.inp')
        appendfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-rmsd.inp')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-rmsd.out')
        appendfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-rmsd.out')
    except:
        pass

    minimizationfiles = []
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-min.pdb')
        minimizationfiles.append('new_' + file.stripDotPDB(file.filename) + '-min.pdb')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-min.psf')
        minimizationfiles.append('new_' + file.stripDotPDB(file.filename) + '-min.psf')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-min.crd')
        minimizationfiles.append('new_' + file.stripDotPDB(file.filename) + '-min.crd')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-min.inp')
        minimizationfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-min.inp')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-min.out')
        minimizationfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-min.out')
    except:
        pass
    solvationfiles = []
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-solv.pdb')
        solvationfiles.append('new_' + file.stripDotPDB(file.filename) + '-solv.pdb')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-solv.psf')
        solvationfiles.append('new_' + file.stripDotPDB(file.filename) + '-solv.psf')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-solv.crd')
        solvationfiles.append('new_' + file.stripDotPDB(file.filename) + '-solv.crd')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-solv.inp')
        solvationfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-solv.inp')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-solv.out')
        solvationfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-solv.out')
    except:
        pass
    try:
        os.stat('crystl_' + file.stripDotPDB(file.filename) + '.str')
        solvationfiles.append('crystl_' + file.stripDotPDB(file.filename) + '.str')
    except:
        pass

    neutralizedfiles = []
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-neutralized.pdb')
        neutralizedfiles.append('new_' + file.stripDotPDB(file.filename) + '-neutralized.pdb')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-neutralized.psf')
        neutralizedfiles.append('new_' + file.stripDotPDB(file.filename) + '-neutralized.psf')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-neutralized.crd')
        neutralizedfiles.append('new_' + file.stripDotPDB(file.filename) + '-neutralized.crd')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-neutralized.inp')
        neutralizedfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-neutralized.inp')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-neutralized.out')
        neutralizedfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-neutralized.out')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-addions.str')
        neutralizedfiles.append('new_' + file.stripDotPDB(file.filename) + '-addions.str')
    except:
        pass

    nmodesfiles = []
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-nmodes.txt')
        nmodesfiles.append('new_' + file.stripDotPDB(file.filename) + '-nmodes.txt')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-nmodes.inp')
        nmodesfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-nmodes.inp')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-nmodes.out')
        nmodesfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-nmodes.out')
    except:
        pass
    try:
        nmtrjnum = getNormalModeMovieNum(file)
        for trj in range(nmtrjnum):
            trj += 1
            nmodesfiles.append('new_' + file.stripDotPDB(file.filename) + '-mtraj_' + str(trj) + '.trj')
        for trj in range(nmtrjnum):
            trj += 1
            nmodesfiles.append('new_' + file.stripDotPDB(file.filename) + '-nma-mainmovie-' + str(trj) + '.pdb')
    except:
        pass
    mdfiles = []
    try:
        os.stat(file.stripDotPDB(file.filename) + '-mdproperties.dat')
        mdfiles.append(file.stripDotPDB(file.filename) + '-mdproperties.dat')
    except: 
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-md.pdb')
        mdfiles.append('new_' + file.stripDotPDB(file.filename) + '-md.pdb')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-mdavg.pdb')
        mdfiles.append('new_' + file.stripDotPDB(file.filename) + '-mdavg.pdb')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-md.psf')
        mdfiles.append('new_' + file.stripDotPDB(file.filename) + '-md.psf')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-md.crd')
        mdfiles.append('new_' + file.stripDotPDB(file.filename) + '-md.crd')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-md-mainmovie.pdb')
        mdfiles.append('new_' + file.stripDotPDB(file.filename) + '-md-mainmovie.pdb')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-md.inp')
        mdfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-md.inp')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-md.out')
        mdfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-md.out')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-md.dcd')
        mdfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-md.dcd')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-md.res')
        mdfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-md.res')
    except:
        pass
    ldfiles = []
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-ld.pdb')
        ldfiles.append('new_' + file.stripDotPDB(file.filename) + '-ld.pdb')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-ld.psf')
        ldfiles.append('new_' + file.stripDotPDB(file.filename) + '-ld.psf')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-ld.crd')
        ldfiles.append('new_' + file.stripDotPDB(file.filename) + '-ld.crd')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-ld-mainmovie.pdb')
        ldfiles.append('new_' + file.stripDotPDB(file.filename) + '-ld-mainmovie.pdb')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-ld.inp')
        ldfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-ld.inp')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-ld.out')
        ldfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-ld.out')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-ld.dcd')
        ldfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-ld.dcd')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-ld.res')
        ldfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-ld.res')
    except:
        pass
    sgldfiles = []
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-sgld.pdb')
        sgldfiles.append('new_' + file.stripDotPDB(file.filename) + '-sgld.pdb')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-sgld.psf')
        sgldfiles.append('new_' + file.stripDotPDB(file.filename) + '-sgld.psf')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-sgld.crd')
        sgldfiles.append('new_' + file.stripDotPDB(file.filename) + '-sgld.crd')
    except:
        pass
    try:
        os.stat('new_' + file.stripDotPDB(file.filename) + '-sgld-mainmovie.pdb')
        sgldfiles.append('new_' + file.stripDotPDB(file.filename) + '-sgld-mainmovie.pdb')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-sgld.inp')
        sgldfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-sgld.inp')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-sgld.out')
        sgldfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-sgld.out')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-sgld.dcd')
        sgldfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-sgld.dcd')
    except:
        pass
    try:
        os.stat('charmm-' + file.stripDotPDB(file.filename) + '-sgld.res')
        sgldfiles.append('charmm-' + file.stripDotPDB(file.filename) + '-sgld.res')
    except:
        pass
    systfiles = []
    done = re.compile('Done')
    fail = re.compile('Failed')
    try:
        solvation_params = solvationParams.objects.filter(pdb=file,selected='y')[0].statusHTML
    except:
        solvation_params = ''
    if fail.search(solvation_params) or done.search(solvation_params):
        systfiles = file.getNonProteinFiles()
    return render_to_response('html/downloadfiles.html',{'filename_list':filename_list,'append_files':appendfiles, 'minimization_files':minimizationfiles,'solvation_files':solvationfiles,'neutralized_files':neutralizedfiles,'md_files':mdfiles,'ld_files':ldfiles,'sgld_files':sgldfiles,'syst_files':systfiles,'nmodes_files': nmodesfiles})

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
    file.checkRequestData(request)
    file.checkForMaliciousCode(jobtype,request)
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
    file.checkRequestData(request)
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
        os.stat("%s/%s/%s" % (charmming_config.user_home,username,filename))
    except:
        return HttpResponse('Oops ... that file no longer exists.')
    response = HttpResponse(mimetype=mimetype)
    response['Content-Disposition'] = 'attachment; filename=%s' %filename
    response.write(file(charmming_config.user_home + "/" + username + "/" + filename, "rb").read())
    return response

#Loads template that contains all the inp/out files for download
def viewProcessFiles(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        file =  Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please select a structure first.")

    seg_list = []
    minsolv_list = []
    nmodes_list = []
    dyn_list = []
    seg_list = []
    redox_list = []

    for f in file.getInputList():
        if f.endswith(".inp"):
            descr = "Input for"
        else:
            descr = "Output of"

        ftype = ""
        try:
           if "-min." in f:
              descr += " minimization."
              minsolv_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
           elif "-solv." in f:
              descr += " solvation."
              minsolv_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
           elif "-neutralize." in f:
              descr += " neutralization."
              minsolv_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
           elif "-nmodes." in f:
              descr += " normal modes."
              nmodes_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
           elif "-md." in f:
              descr += " molecular dynamics (MD)." 
              dyn_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
           elif "-sgld." in f:
              descr += " self-guided Langevin dynamics (SGLD)."
              dyn_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
           elif "-ld." in f:
              descr += " Langevin dynamics (LD)."
              dyn_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
           elif "-rmsd." in f:
              descr += " RMSD calculation."
              minsolv_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
           elif "redox-" in f:
              descr += " oxidation/reduction calculations"
              redox_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
           else:
              descr += " segment setup and generation."
              seg_list.append( (f,time.ctime(os.stat(file.location + '/' + f)[stat.ST_ATIME])[4:],descr) )
        except:
           pass

    return render_to_response('html/viewprocessfiles.html', {'minsolv_list': minsolv_list,'dyn_list': dyn_list,'seg_list': seg_list,'nmodes_list': nmodes_list,'username': request.user.username, 'redox_list': redox_list})

def getProtoPDB(file):
    pre = re.compile("view-([a-z])-(\d+)-([a-z]+).pdb")
    plist = []

    pdbdir = file.location + "proto_" + file.stripDotPDB(file.filename)
    pdblist = glob.glob(pdbdir + "/view-*.pdb")
    for pdb in pdblist:
        m = pre.search(pdb)
        if not m:
            continue
        spdb = pdb.replace(file.location, "")
        if m.group(3) in ["hsd","hse","hsp"]:
            restp = "Histadine"
        elif m.group(3) == "lys":
            restp = "Lysine"
        elif m.group(3) == "asp":
            restp = "Aspartic acid"
        elif m.group(3) == "glu":
            restp = "Glutamic acid"
        plist.append( (spdb,m.group(1),m.group(2),restp) )

    return plist

#Lets user view PDB through jmol/Chemaxon
def visualize(request,filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        username = request.user.username
        file =  Structure.objects.filter(owner=request.user,selected='y')[0]
        file_list = file.getFileList()
        het_list = file.getNonGoodHetPDBList()
        tip_list = file.getGoodHetPDBList()
        dashf_list = file.getDashFPDBList()
        protein_list = file.getProteinSegPDBList()
    except:
        return HttpResponse("No PDB Uploaded")

    # sometimes, the same file name can occur in both the protein list and
    # the dashf list ... this needs to get cleaned up. We *should* use a set
    # union for this, but to maintain Python 2.3 compatibility we use this
    # hackish method instead.
    for fname in dashf_list:
        if fname in protein_list:
            protein_list.remove(fname)
        elif fname in tip_list:
            tip_list.remove(fname)
        elif fname in het_list:
            het_list.remove(fname)

    append_pdb = "new_" + file.stripDotPDB(file.filename) + "-final.pdb" 
    solv_pdb = "new_" + file.stripDotPDB(file.filename) + "-solv.pdb" 
    neu_pdb = "new_" + file.stripDotPDB(file.filename) + "-neutralized.pdb" 
    min_pdb = "new_" + file.stripDotPDB(file.filename) + "-min.pdb" 
    md_pdb = "new_" + file.stripDotPDB(file.filename) + "-md.pdb" 
    md_movie_pdb = "new_" + file.stripDotPDB(file.filename) + "-md-mainmovie.pdb" 
    ld_pdb = "new_" + file.stripDotPDB(file.filename) + "-ld.pdb" 
    ld_movie_pdb = "new_" + file.stripDotPDB(file.filename) + "-ld-mainmovie.pdb" 
    sgld_pdb = "new_" + file.stripDotPDB(file.filename) + "-sgld.pdb" 
    sgld_movie_pdb = "new_" + file.stripDotPDB(file.filename) + "-sgld-mainmovie.pdb" 
    nmtrjnum = getNormalModeMovieNum(file)
    nmodes_list = []
    for trj in range(nmtrjnum):
        trj += 1
        nmodes_list.append('new_' + file.stripDotPDB(file.filename) + '-nma-mainmovie-' + str(trj) + '.pdb')

    proto_list = getProtoPDB(file)
    return render_to_response('html/visualize.html', {'filename': filename,'username':username,'file_list':file_list,\
                             'het_list':het_list,'tip_list':tip_list,'protein_list':protein_list, 'min_pdb':min_pdb,\
                             'solv_pdb':solv_pdb,'neu_pdb':neu_pdb,'md_pdb':md_pdb,'ld_pdb':ld_pdb,'sgld_pdb':sgld_pdb,\
                             'append_pdb':append_pdb, 'dashf_list':dashf_list,'md_movie_pdb':md_movie_pdb,\
                             'ld_movie_pdb':ld_movie_pdb,'sgld_movie_pdb':sgld_movie_pdb,'proto_list': proto_list,'nmodes_list':nmodes_list})

#Let's user view PDB through jmol
def jmol(request,filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        username = request.user.username
        return render_to_response('html/jmol.html', {'filename': filename,'username':username,'segid':'NA','resid':'NA' })
    except:
        return HttpResponse("No PDB Uploaded")

def jmolHL(request,filename,segid,resid):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    fp = open("/tmp/foo", "a+")
    fp.write("Got %s %s %s\n" % (filename,segid,resid))
    fp.close()
    try:
        username = request.user.username
        return render_to_response('html/jmol.html', {'filename': filename,'username':username, 'segid': segid, 'resid': resid })
    except:
        return HttpResponse("No PDB Uploaded")

#Let's user view PDB through chemaxon
def chemaxon(request,filename):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        username = request.user.username
        return render_to_response('html/chemaxon.html', {'filename': filename,'username':username })
    except:
        return HttpResponse("No PDB Uploaded")

#Allows the user to see what processes their PDBs are undergoing
def viewstatus(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    try:
        file = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return render_to_response('html/statusreport.html', {'structure': None })
    done = re.compile('Done')

    try:
        md_param = mdParams.objects.filter(pdb=file,selected='y')[0]
    except:
        md_param = None
    try:
        ld_param = ldParams.objects.filter(pdb=file,selected='y')[0]
    except:
        ld_param = None
    try:
        sgld_param = sgldParams.objects.filter(pdb=file,selected='y')[0]
    except:
        sgld_param = None
    try:
        nma_param = nmodeParams.objects.filter(pdb=file,selected='y')[0]
    except:
        nma_param = None
    try:
        redox_param = redoxParams.objects.filter(pdb=file,selected='y')[0]
    except:
        redox_param = None

    if(nma_param and nma_param.make_nma_movie):
        combineNmaPDBsForMovie(file)
        if (nma_param.nma_movie_status):
            nma_param.make_nma_movie = False
            nma_param.make_nma_req = False
            nma_param.save()
    file.updateActionStatus()
    if(md_param and md_param.make_md_movie): 
        combinePDBsForMovie(file,'md')
        #If the movie status has been updated then reset the req
        #and movie stuff
        if (md_param.md_movie_status):
            md_param.make_md_movie = False
            md_param.make_md_req = False
            md_param.save()
        #If there is a ld_param and if the main job is finished and the movie has not been created yet
    if(ld_param and ld_param.make_ld_movie):
        combinePDBsForMovie(file,'ld')
        if (ld_param.ld_movie_status):
            ld_param.make_ld_movie = False
            ld_param.make_ld_req = False
            ld_param.save()
    if(sgld_param and sgld_param.make_sgld_movie):
        combinePDBsForMovie(file,'sgld')
        if (sgld_param.sgld_movie_status):
            sgld_param.make_sgld_movie = False
            sgld_param.make_sgld_req = False
            sgld_param.save()
    #file.updateActionStatus()
    try:
        mini_param = minimizeParams.objects.filter(pdb=file,selected='y')[0]
    except:
        mini_param = None
    try:
        solv_param = solvationParams.objects.filter(pdb=file,selected='y')[0]
    except:
        solv_param = None
    return render_to_response('html/statusreport.html', {'structure': file, 'mini_param':mini_param, 'solv_param':solv_param, 'md_param':md_param, \
                                                         'nma_param':nma_param,'ld_param':ld_param,'sgld_param':sgld_param, 'redox_param': redox_param})

#Used to report an error. Used in calling problem.html
def reportError(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    structure.models.Structure.checkRequestData(request)
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



# edit the protonation state of the currently selected protein
# Tim Miller, 04-02-2008
def editProtonation(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    file = Structure.objects.filter(owner=request.user,selected='y')[0]
    file.checkRequestData(request)
    plines = {}
    hist_list = {}

    # step 1, make list of residues to change...  histadine changes must be done right in the
    # PDB file -- all others have a patch added to a stream file.

    patchre = re.compile("patch_([A-Z]+)_([0-9]+)")
    for kname in request.POST.keys():
        m = patchre.match(kname)
        if m:
            newres = request.POST[kname]
            segid  = m.group(1)
            resid  = int(m.group(2))
            if newres in ['HSD','HSE','HSP']:
                try:
                    hist_list[segid].append((resid,newres))   
                except:
                    hist_list[segid] = []
                    hist_list[segid].append((resid,newres))
            elif newres in ['ASPP','GLUP','LSN']:
                try:
                    plines[segid].append("PATCH %s %s %d\n" % (newres,segid,resid))
                except:
                    plines[segid] = []
                    plines[segid].append("PATCH %s %s %d\n" % (newres,segid,resid))
    for seg in plines.keys():
        patchfp = open(file.location + "pproto_" + file.stripDotPDB(file.filename) + "-" + seg.lower() + ".str", "w+")
        patchfp.write("* Patch user-specified residues to protonate in seg %s.\n" % seg)
        patchfp.write("* Note: histadine changes are done by editing the sequence directly.\n")
        patchfp.write("*\n\n")
        for pl in plines[seg]:
            patchfp.write(pl)
        patchfp.write("return\n")
        patchfp.close()

    # step 2: any Histadine patches must be done in the correct PDB file.
    for seg in hist_list.keys():
        ifname = file.location + "new_" + file.stripDotPDB(file.filename) + "-" + seg.lower() + ".pdb"
        bfname = file.location + "bak-new_" + file.stripDotPDB(file.filename) + "-" + seg.lower() + ".pdb"
        ofname = file.location + "prot_" + file.stripDotPDB(file.filename) + "-" + seg.lower() + ".pdb"

        ifp = open(ifname, 'r')
        ofp = open(ofname, 'w+')

        resids = [z[0] for z in hist_list[seg]]
        newval = [z[1] for z in hist_list[seg]]
        for line in ifp:
            line = line.strip()

            if not line.startswith('ATOM'):
                ofp.write("%s\n" % line)
                continue

            larr = string.splitfields(line)
            try:
                cresid = int(larr[5])
                curres = larr[3]
                ridx = resids.index(cresid)
            except:
                ofp.write("%s\n" % line)
                continue
            # If we've made it this far, we know that we have a histadine, check if the
            # user wants to change it.
            if curres != newval[ridx]:
                # they most certainly do
                line = line.replace(curres,newval[ridx])
            ofp.write("%s\n" % line)         


        ifp.close()
        ofp.close()
        os.system("mv %s %s" % (ifname,bfname))
        os.system("mv %s %s" % (ofname,ifname))
    return HttpResponse("Sequence edit and patches prepared according to the selections given.")

# respond for request...
def viewprotores(request,segid,resid):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    return HttpResponse("To do...")

# This subroutine kicks of jobs to visualize the protonizable residues
def subProtoVis(plist,epdb,request):
    if epdb:
        file = Structure.objects.filter(owner=request.user,filename=epdb)[0]
    else:
        file = Structure.objects.filter(owner=request.user,selected='y')[0]
    try:
        os.stat(file.location + "proto_" + file.stripDotPDB(file.filename))
    except:
        try:
            os.mkdir(file.location + "proto_" + file.stripDotPDB(file.filename))
            os.system("chmod g+w " + file.location + "proto_" + file.stripDotPDB(file.filename))
        except:
            # we won't be able to create views, but that's not worth wrecking the train for
            # so just return
            return
    else:
        # OK directory already exists ... assume everything is there and just return
        return

    # OK, for each stinking protonizable residue, we need to create an input script that
    # cuts a 5 A radius around it and writes out the resultant PDB

    charmm_inp = file.makeCHARMMInputHeader("Make PDB for visualizing titratable group",{})
    # delete all BUT the surrounding 8 A
    for item in plist:
        segid = item[1]
        resid = item[2]
        charmm_inp += "\n\nopen read card unit 20 name " + file.location + "new_" + file.stripDotPDB(file.filename) + "-" + segid.lower() + ".pdb\n"
        charmm_inp += "read sequ pdb unit 20\n"
        charmm_inp += "rewind unit 20\n"
        charmm_inp += "generate pseg setu\n"
        charmm_inp += "read coor pdb unit 20\n"
        charmm_inp += "close unit 20\n\n"
        charmm_inp += "! get rid of hydrogens so they don't clutter up the view...\n"
        charmm_inp += "dele atom sort sele hydrogen end\n\n"
        charmm_inp += "ic para\n"
        charmm_inp += "ic fill preserve\n"
        charmm_inp += "ic build\n\n"
        charmm_inp += "! cut a sphere of radius 3 around the central point in the residue to be protonated\n"
        charmm_inp += "coor stat sele resid " + str(resid) + " end\n"
        charmm_inp += "! write out PDB\n"
        charmm_inp += "open write card unit 20 name " + file.location + "proto_" + file.stripDotPDB(file.filename) + "/view-" + segid.lower() + "-" + str(resid) + "-" + item[0] + ".pdb\n"
        charmm_inp += "write coor pdb unit 20 sele .byres. (point ?xave ?yave ?zave cut 8.0) end\n"
        charmm_inp += "close unit 20\n"
        charmm_inp += "dele atom sort sele all end\n\n"

    charmm_inp += "stop\n"

    fname = file.location + "proto_" + file.stripDotPDB(file.filename) + "/charmm-makepdb.inp"
    outfp = open(fname, "w")
    outfp.write(charmm_inp)
    outfp.close()
    si = schedInterface()
    si.submitJob(file.owner.id,file.location + "proto_" + file.stripDotPDB(file.filename),[fname])

#Used in editing PDB information, receives a GET request and
#changes the PDB information based on what was changed
def editPDB(request,edit_pdb):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    input.checkForMaliciousCode(edit_pdb,request)
    if(edit_pdb):
        file = structure.models.Structure.objects.filter(owner=request.user,name=edit_pdb)[0]
        warnings = file.ifWarningsExist()
        try:
            oldfile = Structure.objects.filter(owner=request.user,selected='y')[0]
            oldfile.selected = ''
            oldfile.save()
            file.selected='y'
            file.save()
        except:
            pass
    else:
        file = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
        warnings = file.ifWarningsExist()
    try:
       if request.GET:
           temp_handle = request.GET['fieldname']
           if(temp_handle == 'title'):
               file.title = request.GET['content']
               file.save()
               return HttpResponse(file.title)
           elif(temp_handle == 'author'):
               file.author = request.GET['content']
               file.save()
               return HttpResponse(file.author)
           elif(temp_handle == 'journal'):
               file.journal = request.GET['content']
               file.save()
               return HttpResponse(file.journal)
    except:
       pass
    if file.append_status != "Done":
       do_editprot = 1
       proto_list = get_proto_res(file)
       if len(proto_list) == 0:
           do_editprot = 0
       else:
           subProtoVis(proto_list,edit_pdb,request)
    else:
       do_editprot = 0
       proto_list = []
    #if file.lessons:
    #    file.lessons.onEditPDB()
    return render_to_response('html/editpdbinfo.html', {'pdb': file, 'proto': do_editprot, 'proto_list': proto_list,'warnings':warnings})

# code for editing multi-model PDBs, should be very similar to editPDB above, except we have a list of
# files that can be edit
def editMultiModel(request,modlist,makeGo,makeBLN):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    # deselect the current PDB and arbitrarily select the first PDB in the list...
    oldfile = Structure.objects.filter(owner=request.user,selected='y')[0]
    oldfile.selected = ''
    oldfile.save()
    file = modlist[0]
    file.selected = 'y'
    file.save()
    # the dictionary 'file_type' in lesson_config.py takes lesson number as its key and the corresponding value is the file type of the files the lesson uploads
    # here we create lesson_obj only if the lesson uploads CRD/PSF files
    try: 
        if file_type[file.lesson_type[6:]].find("CRD") != -1 or file_type[file.lesson_type[6:]].find("PSF") != -1:
            try:
                lesson_obj = eval(file.lesson_type+'.models.'+file.lesson_type.capitalize()+'.objects.filter(user = file.owner,id=file.lesson_id)[0].onFileUpload(request.POST)')
            except:
                pass
    except:
        pass   

    # take care of generating protonation states. NB this could wind up being
    # a bit intense if there are many models so we may want to rethink this
    # strategy...
    for file in modlist:
       #each model must get a filename for its topology/param files
       file.getRtfPrm(request,makeGo,makeBLN)
       proto_list = get_proto_res(file)
       if len(proto_list) > 0:
           # We don't actually let the user modify protonation on the multimodel
           # page (just yet), but we do need toput in the request to generate the
           # visualization PDBs in case they decide to edit later.
           subProtoVis(proto_list,file.filename,request)
    npdbs = len(modlist)
    return render_to_response('html/editmultimodel.html', {'pdbs': modlist, 'npdbs': npdbs})



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
    user_pdbs = Structure.objects.filter(owner=request.user)
    return render_to_response('html/viewpdbs.html', {'user_pdbs': user_pdbs})

#This is for changing the currently selected PDB
#on the SELECT/EDIT PDBs page
def switchpdbs(request,switch_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        oldfile = Structure.objects.filter(owner=request.user,selected='y')[0]
        oldfile.selected = ''
        oldfile.save()
    except:
        pass
    newfile = Structure.objects.filter(owner=request.user,filename=switch_id)[0] 
    newfile.selected = 'y'
    newfile.save()
    return render_to_response('html/switchpdb.html',{'oldfile':oldfile,'newfile':newfile})

#This calculates the energy
def energyform(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    structure.models.Structure.checkRequestData(request)
    #chooses the file based on if it is selected or not
    try:
        file =  Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    os.chdir(file.location)
    temp = file.stripDotPDB(file.filename)
    #creates a list of filenames associated with the PDB
    filename_list = file.getLimitedFileList("blank")
    #filename_list = file.getFileList()
    solv_pdb = "new_" + file.stripDotPDB(file.filename) + "-solv.pdb"
    neu_pdb = "new_" + file.stripDotPDB(file.filename) + "-neutralized.pdb"
    min_pdb = "new_" + file.stripDotPDB(file.filename) + "-min.pdb"
    md_pdb = "new_" + file.stripDotPDB(file.filename) + "-md.pdb"
    ld_pdb = "new_" + file.stripDotPDB(file.filename) + "-ld.pdb"
    sgld_pdb = "new_" + file.stripDotPDB(file.filename) + "-sgld.pdb"
    seg_list = file.segids.split(' ')
    het_list = file.getNonGoodHetPDBList()
    tip_list = file.getGoodHetPDBList()
    seg2_list = file.getProteinSegPDBList()
    protein_list = [x for x in seg2_list if x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb")]
    go_list = [x for x in seg2_list if x.endswith("-go.pdb") or x.endswith("-go-final.pdb")]
    bln_list = [x for x in seg2_list if x.endswith("-bln.pdb") or x.endswith("-bln-final.pdb")]
    nucleic_list = [x for x in seg2_list if x.endswith("-dna.pdb") or x.endswith("-dna-final.pdb") or x.endswith("-rna.pdb") or x.endswith("-rna-final.pdb")]
    userup_list = [x for x in seg2_list if not ( x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb") or x.endswith("-dna.pdb") or x.endswith("-dna-final.pdb") \
                     or x.endswith("-rna.pdb") or x.endswith("-rna-final.pdb") or x.endswith("het.pdb") or x.endswith("het-final.pdb") or x.endswith("-go.pdb") \
                     or x.endswith("-go-final.pdb") or x.endswith("-bln.pdb")  or x.endswith("-bln-final.pdb") ) ]
    disulfide_list = file.getPDBDisulfides()

    # Tim Miller: 03-09-2009 we need to pass a segpatch_list (containing
    # both the list of segments and their default NTER patches) to the
    # template.
    propatch_list = []
    nucpatch_list = []
    for seg in seg_list:
       	if seg.endswith("-pro"):
            defpatch = checkNterPatch(file,seg)
            propatch_list.append((seg,defpatch))
        elif seg.endswith("-dna") or seg.endswith("-rna"):
            nucpatch_list.append((seg,"5TER","3TER"))

    energy_lines = ''
    try:
        energyfp = open('energy-' + file.stripDotPDB(file.filename) + '.txt','r')
        for line in energyfp:
            energy_lines += line
        energyfp.close()
    except:
        pass
    done = re.compile('Done')
    for i in range(len(filename_list)):
        #First check and see if the selected choices are segids
        #Otherwise see if it is a solvated/minimized PDB
        try:
            tempid = request.POST[filename_list[i]]
            filename = filename_list[i]
        except:
            try:
                tempid = request.POST['solv_or_min']
                filename = request.POST['solv_or_min']
            except:
                tempid = "null"
        if(tempid!="null"):
            seg_list = file.segids.split(' ')
            try:
                if(request.POST['usepatch']):
                    file.handlePatching(request.POST)
            except:
                #If there is no patch, make sure patch_name is zero 
                file.patch_name = ""
                file.save()
            scriptlist = []
            if(filename != min_pdb and filename != solv_pdb and filename != md_pdb and filename != ld_pdb and filename != sgld_pdb and filename != neu_pdb):
                seg_list = append_tpl(request.POST,filename_list,file,scriptlist)
                energy_this_file = 'new_' + file.stripDotPDB(file.filename) + "-final.pdb"
                return calcEnergy_tpl(request,file,seg_list,energy_this_file,scriptlist)
            else:
                return calcEnergy_tpl(request,file,seg_list,filename,scriptlist)

    doCustomShake = 1
    if file.ifExistsRtfPrm() < 0:
        doCustomShake = 0
    trusted = isUserTrustworthy(request.user)
    return render_to_response('html/energyform.html', {'filename_list': filename_list,'propatch_list':propatch_list,'solv_pdb':solv_pdb,'neu_pdb':neu_pdb, \
      'min_pdb':min_pdb,'md_pdb':md_pdb,'ld_pdb':ld_pdb,'sgld_pdb':sgld_pdb,'het_list':het_list,'tip_list':tip_list,'protein_list':protein_list, 'file': file, \
      'docustshake': doCustomShake,'energy_lines':energy_lines,'trusted':trusted,'disulfide_list': disulfide_list, 'seg_list': seg_list, \
      'nucpatch_list': nucpatch_list, 'nucleic_list': nucleic_list, 'userup_list': userup_list, 'go_list': go_list, 'bln_list': bln_list})


def calcEnergy_tpl(request,file,seg_list,min_pdb,scriptlist): 
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    postdata = request.POST
    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'] = file.getTopologyList()
    template_dict['parameter_list'] = file.getParameterList()
    template_dict['filebase'] = file.stripDotPDB(file.filename)
    template_dict['input_file'] = file.stripDotPDB(min_pdb)
    template_dict['file_location'] = file.location

    try:
       solvate_implicitly = postdata['solvate_implicitly']
    except:
       solvate_implicitly = 0


    try:
        oldeo = energyParams.objects.filter(pdb = file, selected = 'y')[0]
        oldeo.selected = 'n'
        oldeo.save()
    except:
        pass
    eobj = energyParams()
    eobj.pdb = file
    eobj.selected = 'y'
    eobj.finale = '0.0' # needed to shut up an error from Django

    #If the user wants to solvate implicitly the scpism line is needed
    #84 will be the scpism number in this program
    solvate_implicitly = 0
    try:
        if(postdata['solvate_implicitly']):
            solvate_implicitly = 1
    except:
        pass
    template_dict['solvate_implicitly'] = solvate_implicitly
    template_dict['fixnonh'] = 0

    try:
        postdata['fixnonh']
        template_dict['fixnonh'] = 1

    except:
        #If there is a topology or paramter file then don't constrain anything
        if file.ifExistsRtfPrm() < 0:
            template_dict['fixnonh'] = 2
    #handles shake
    template_dict['shake'] = request.POST.has_key('apply_shake')
    if request.POST.has_key('apply_shake'):
        template_dict['which_shake'] = postdata['which_shake']
        template_dict['qmmmsel'] = postdata['qmsele']

        if postdata['which_shake'] == 'define_shake':
            template_dict['shake_line'] = postdata['shake_line']

            if postdata['shake_line'] != '':
                file.checkForMaliciousCode(postdata['shake_line'],postdata)

    template_dict['restraints'] = ''
    try:
        postdata['apply_restraints']
        template_dict['restraints'] = file.handleRestraints(request)
    except:
        pass

    # check to see if PBC needs to be used -- if so we have to set up Ewald
    template_dict['usepbc'] = ''
    template_dict['dopbc'] = 0
    try:
        if postdata['usepbc']:
            template_dict['usepbc'] =  postdata['usepbc']
            eobj.usepbc = 'y'
            if file.solvation_structure != 'sphere':
                dopbc = 1
                template_dict['dopbc'] = 1
            else:
                dopbc = 0
        else:
            eobj.usepbc = 'n'
            dopbc = 0
    except:
        eobj.usepbc = 'n'
        dopbc = 0

    template_dict['solvation_structure'] = file.solvation_structure
    template_dict['relative_boundary'] = 0
    if dopbc:
        relative_boundary = 0
        if file.solvation_structure != '' and file.crystal_x < 0:
            relative_boundary = 1

        template_dict['relative_boundary'] = relative_boundary
        template_dict['dim_x'] = str(file.crystal_x)
        template_dict['dim_z'] = str(file.crystal_z)
        #template_dict['greaterval'] = str(greaterval)

        # set up images
        if file.solvation_structure == '' or solvate_implicitly:
             pass
        else:
            # we should have a solvation file to read from
            try:
                os.stat(file.location + "new_" + file.stripDotPDB(file.filename) + ".xtl")
            except:
                # need to throw some sort of error ... for now just toss cookies
                return HttpResponse("Oops ... transfer file not found.")

    template_dict['useqmmm'] = postdata.has_key("useqmmm")
    if postdata.has_key("useqmmm"):
        eobj.useqmmm = 'y'

        # validate input

        if postdata['qmmm_exchange'] in ['HF','B','B3']:
            exch = postdata['qmmm_exchange']
        else:
            exch = 'HF'
        if postdata['qmmm_correlation'] in ['None','LYP']:
            corr = postdata['qmmm_correlation']
        else:
            corr = 'None'
        if postdata['qmmm_basisset'] in ['STO-3G','3-21G*','6-31G*']:
            bs = postdata['qmmm_basisset']
        else:
            bs = 'sto3g'
        file.checkForMaliciousCode(postdata['qmsele'],postdata)
        qmsel = postdata['qmsele']
        if qmsel == '':
            qmsel = 'resid 1'
        if postdata['qmmm_charge'] in ['-5','-4','-3','-2','-1','0','1','2','3','4','5']:
            charge = postdata['qmmm_charge']
        else:
            charge = '0'
        if postdata['qmmm_multiplicity'] in ['0','1','2','3','4','5','6','7','8','9','10']:
            multi = postdata['qmmm_multiplicity']
        else:
            multi = '0'
        try:
            if int(postdata['num_linkatoms']) > 0:
                linkatoms = handleLinkAtoms(file,postdata)
            else:
                linkatoms = None
        except:
            linkatoms = None
        eobj.qmmmsel = qmsel

#        charmm_inp = makeQChem(charmm_inp, file, exch, corr, bs, qmsel, "SP", charge, multi, file.stripDotPDB(min_pdb) + ".crd", linkatoms)
        template_dict = makeQChem_tpl(template_dict, file, exch, corr, bs, qmsel, "SP", charge, multi, file.stripDotPDB(min_pdb) + ".crd", linkatoms)

#        if "achtung, bad QM" in charmm_inp:
#            return HttpResponse("<b>Bad QM parameter specified!</b>\n")
#        elif "you changed the QM region" in charmm_inp:
#            return HttpResponse("<b>You changed the QM region from the previous structure. That's not allowed!</b>\n")
    else:
        eobj.useqmmm = 'n'

    eobj.save()
    if file.lesson_type:
        lessonaux.doLessonAct(file,"onEnergySubmit",request.POST)

    template_dict['headqmatom'] = 'blankme'
    if eobj.useqmmm == 'y':
        headstr = writeQMheader("", "SELE " + qmsel + " END")
        template_dict['headqmatom'] = headstr.strip()

    template_dict['output_name'] = "new_" + file.stripDotPDB(file.filename) + "-energy"
    t = get_template('%s/mytemplates/input_scripts/calcEnergy_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    if postdata.has_key("useqmmm"):
        if "achtung, bad QM" in charmm_inp:
            return HttpResponse("<b>Bad QM parameter specified!</b>\n")
        elif "you changed the QM region" in charmm_inp:
            return HttpResponse("<b>You changed the QM region from the previous structure. That's not allowed!</b>\n")

    user_id = file.owner.id
    os.chdir(file.location)
    energy_input = "charmm-" + file.stripDotPDB(file.filename) + "-energy.inp"
    scriptlist.append(energy_input)
    energy_output = "charmm-" + file.stripDotPDB(file.filename) + "-energy.out"
    inp_out = open(file.location + energy_input,'w')
    inp_out.write(charmm_inp)
    inp_out.close()
    if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
        return HttpResponse(generateHTMLScriptEdit(charmm_inp,scriptlist,'energy'))
    else:
        si = schedInterface()
        enerjobID = si.submitJob(user_id,file.location,scriptlist)
        #Append/energy needs to be doen before energy can be parsed
        done = re.compile('Done')
        failed = re.compile('Failed')
        sstring = si.checkStatus(enerjobID)
        while(not done.search(statsDisplay(sstring,enerjobID)) and not failed.search(statsDisplay(sstring,enerjobID))):
            sstring = si.checkStatus(enerjobID)
        file.save()
        return parseEnergy(file,energy_output,eobj)

def parseEnergy(file,output_filename,enerobj=None):
    time.sleep(2)

    outfp = open(file.location + output_filename,'r')
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
            try:
                enerobj.save()
            except:
                # usually means the FP was rounded; should be harmless
                pass
        if(dash_occurance == 2 or (initiate > 0 and line.strip().startswith('CHARMM>'))):
            if(line.startswith('CHARMM>')):
                break
            energy_lines += line
            break
        if(ener.search(line) or initiate > 0):
            if(initiate == 0):
                initiate = 1
            energy_lines += line
        if(dash.search(line) and initiate > 0):
            dash_occurance += 1
    outfp.close()
    writefp = open(file.location + 'energy-' + file.stripDotPDB(file.filename) + '.txt','w')
    writefp.write(energy_lines)
    writefp.close()

    if file.lesson_type:
        lessonaux.doLessonAct(file=file,function="onEnergyDone",finale=enerobj.finale)

    return render_to_response('html/displayenergy.html',{'linelist':energy_lines})

def getSegs(Molecule,Struct,auto_append=False):
    Struct.save()
    for seg in Molecule.iter_seg():
        newSeg = structure.models.Segment()
        newSeg.structure = Struct
        if auto_append:
            newSeg.is_appended = 'y'
        else:
            newSeg.is_appended = 'n'
        newSeg.name = seg.segid
        newSeg.type = seg.segType
        newSeg.save()

def newupload(request, template="html/fileupload.html"):
    doc =\
    """
    Handle a new file being uploaded.
    """

    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        tempfile.checkRequestData(request)
    except:
        pass
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
    # end of checking for lesson   

    try:
        request.FILES['pdbupload'].name
        pdb_uploaded_by_user = 1
    except:
        pdb_uploaded_by_user = 0

    if pdb_uploaded_by_user:
        filename = request.FILES['pdbupload'].name
        filename = filename.lower()
        if lesson_obj:
            struct.lesson_type = lessontype
            struct.lesson_id = lessonid

        # figure out where this file is going to go and put it there.
        tmpdname = filename.split('.')[0]
        dname = tmpdname
        version = -1
        while os.path.exists(location + '/' + dname):
            version += 1
            dname = tmpdname + "-" + str(version)

        os.mkdir(location + '/' + dname)
        fullpath = location + '/' + dname + '/' + filename
        temp = open(fullpath, 'w')
        for fchunk in request.FILES['pdbupload'].chunks():
            temp.write(fchunk)
        temp.close()

        if filename.endswith('crd') or filename.endswith('cor'):
            # set up the new structure object
            struct = structure.models.Structure()
            struct.name = dname

            crd = charmming.io.crd.CRDFile(fullpath)
            thisMol = crd.iter_models.next()
            getSegs(thisMol,struct,auto_append=True)
        else:
            pdb = charmming.io.pdb.PDBFile(fullpath)
            if len(pdb.keys()) > 1:
                mnum = 1
                for model in pdb.iter_models():
                    mdname = dname + "-mod-%d" % mnum
                    os.mkdir(location + '/' + mdname)
 
                    struct = structure.models.Structure()
                    struct.name = mdname
                    struct.append_status = 'n'
                    struct.getHeader(pdb.header)
                    struct.owner = request.user

                    model.parse()
                    getSegs(model,struct)
                    
                    pfname = location + dname + '/pdbpickle.dat'
                    pickleFile = open(pfname,'w')
                    cPickle.dump(pdb['model%d' % (mnum-1)],pickleFile)
                    pickleFile.close()
                    struct.pickle = pfname

                    struct.save()
                    mnum += 1
            else:
                struct = structure.models.Structure()
                struct.name = dname
                struct.append_status = 'n'
                struct.getHeader(pdb.header)
                struct.owner = request.user

                thisMol = pdb.iter_models().next()
                thisMol.parse()
                getSegs(thisMol,struct)

                pfname = location + dname + 'pdbpickle.dat'
                pickleFile = open(pfname,'w')
                cPickle.dump(pdb['model0'],pickleFile)
                pickleFile.close()
                struct.pickle = pfname


        # unselect the existing structure
        try:
            oldfile = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
            oldfile.selected = ''
            oldfile.save()
        except:
            pass

        # set this structure as selected
        struct.selected = 'y'
        struct.save()
        return HttpResponseRedirect('/charmming/editpdbinfo/'+struct.name)

    #The below code handles custom sequences
    elif request.POST.has_key('sequ'):
        file = structure.models.Structure()
        file.owner = u
        file.location = location

        if lesson_obj:
            file.lesson_type = lessontype
            file.lesson_id = lessonid
            file.save()
        sequ = request.POST['sequ']
        sequ = sequ.upper()
        file.filename = 'sequ-' + str(filepid) + "-1.pdb"
        #check of parameter/topology files
        file.getRtfPrm(request,makeGoModel,makeBLNModel)
        file.segids = 'sequ-pro'
        sequ_residues = sequ.split()
        plines = []
        resid = 0

        for res in sequ_residues:
            resid += 1
            if res in ['ASPP','GLUP','LSN']:
                try:
                    plines.append("PATCH %s %s %d\n" % (res,'sequ',resid))
                except:
                    plines = []
                    plines.append("PATCH %s %s %d\n" % (res,'sequ',resid))
                # Standard CHARMM topology files do not understand ASPP,GLUP,
                # or LSN so they must get changed to ASP, GLU, and LYS
                # respectively. Afterwards the patches can be applied.
                if res in ['ASPP']:
                    sequ_residues[(resid-1)] = 'ASP'
                elif res in ['GLUP']:
                    sequ_residues[(resid-1)] = 'GLU'
                elif res in ['LSN']:
                    sequ_residues[(resid-1)] = 'LYS'
        patchfp = open(file.location + "pproto_" + file.stripDotPDB(file.filename) + "-" + 'sequ' + ".str", "w+")
        patchfp.write("* Patch user-specified residues to protonate in seg %s.\n" % 'sequ')
        patchfp.write("* Note: histadine changes are done by editing the sequence directly.\n")
        patchfp.write("*\n\n")
        for pl in plines:
            patchfp.write(pl)
        patchfp.write("return\n")
        patchfp.close()

        #The written file for a sequence PDB should not match file.filename
        #Doing this creates the best compatibility with other parts
        #of the interface. The parser usually renames the PDB in the format
        #below but the parser does not get called on custom sequencing
        temp_sequ_handle = open(file.location + "new_" + 'sequ-' + str(file.pid) + "-1-sequ-pro.pdb",'w')
        sequ = " ".join(sequ_residues)
        temp_sequ_handle.write(sequ)
        temp_sequ_handle.close()
        file.title = "Custom Sequence"
        file.author = file.owner.username
        file.journal = "None"
        # The newly uploaded PDB becomes selected by default
        # but first the previously selected PDB must become unselected
        try:
            test_select = structure.models.Structure.objects.filter(owner=file.owner,selected='y')[0]
            test_select.selected=''
            test_select.save()
        except:
            pass
        file.selected='y'
        file.save()
        if lesson_obj:
            file.lesson_type = lessontype
            file.lesson_id = lessonid
            file.save()
            try:
                lesson_obj.onFileUpload(request.POST)
            except:
                # TODO -- handle failures
                pass
	return HttpResponseRedirect('/charmming/editpdbinfo/'+file.filename)

    # handle downloading the PDB from PDB.org
    elif request.POST.has_key('pdbid'):
        pdbid = request.POST['pdbid']
        pdbid = pdbid.strip()
        pdbid = pdbid.lower()
        if(error.search(pdbid)):
            return HttpResponse('PDB filename invalid')
        cif = re.compile('\Zcif')
        filename = pdbid

        #Because .cif files are standard and PDBs are not, the following block
        #downloads a cif file then converts it to PDB
        dne = re.compile("does not exist")
        try:
                conn = HTTPConnection("www.pdb.org")
                conn.request("GET", "/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s" % pdbid)
                resp = conn.getresponse()
                outfp = open(location + "/" + filename + ".pdb", 'w+')
                if resp.status != 200:
                        prob_string = "The PDB server returned error code %d." % resp.status
                        return render_to_response('html/problem.html',{'prob_string': prob_string})
                cif_file = resp.read()
                if dne.search(cif_file):
                        outfp.close()
                        return render_to_response('html/problem.html',{'prob_string': 'The PDB server reports that the structure does not exist.'})

                outfp.write(cif_file)
                outfp.close()
                conn.close()
        except:
                # print an error page...
                return render_to_response('html/problem.html',{'prob_string': 'There was an exception processing the file -- contact the server administrator for more details.'})

        os.chdir(location)
        filename = filename + ".pdb"

        modelinfo = getModels(filename,u,location,lessontype,lessonid)
        models = modelinfo['models']
        errors = modelinfo['errors']
        if len(models) < 1:
            # parse error
            return render_to_response('html/problem.html', {'prob_string': 'The PDB parser did not find any valid models.'})
        elif len(models) == 1:
            file = models[0]
            #check of parameter/topology files
            file.getRtfPrm(request,makeGoModel,makeBLNModel)
            if lesson_obj:
                file.lesson_type = lessontype
                file.lesson_id = lessonid
                file.save()
                try:
                    lesson_obj.onFileUpload(request.POST)
                except:
                    # TODO, scream
                    pass

            return HttpResponseRedirect('/charmming/editpdbinfo/'+file.filename)
        else:
            # editMultiModel does not return, but renders the HTML to the user
            return editMultiModel(request,models,makeGoModel,makeBLNModel)

    form = structure.models.PDBFileForm()
    return render_to_response('html/fileupload.html', {'form': form} )

#Returns a list of segids in a PSF as a list
def scanPsfForSegs(file):
    psf = file.location + 'psf-' + file.stripDotCRD(file.filename) + '.psf'
    psf_file = open(psf,'r')
    natom = re.compile('NATOM')
    nbond = re.compile('NBOND')
    reached_natom = 0
    base_segid = "not_a_real_segid"
    test_segid = "this_is_also_not_a_real_segid"
    segid_list = []
    for line in psf_file:
        if reached_natom == 0:
            if(natom.search(line)):
                reached_natom = 1
            continue
        psf_fields = line.split()
        if psf_fields == []:
            continue
        if nbond.search(line):
            break

        test_segid = psf_fields[1].lower()
        if test_segid != base_segid:
            base_segid = test_segid
            segid_list.append(base_segid)
    file.segids = " ".join(segid_list)
    file.save()
    return segid_list


# This algorithm isn't efficient b/c it checks all atoms unnecessarily
def get_proto_res(file):
    protonizable = []

    # all residues that we know how to protonize...
    possible_p = ['HIS','HSD','HSE','HSP','ASP','GLU','LYS']

    seg_list = file.segids.split(' ')
    for seg in seg_list:
        fpath = file.location + "new_" + file.stripDotPDB(file.filename) + "-" + seg + ".pdb"
        try:
            fp = open(fpath, "r")
            for line in fp:
                line = line.strip()
                if not line.startswith('ATOM'):
                    continue
                larr = string.splitfields(line)

                # Is this bad boy something we can protonate?
                if not larr[3] in possible_p:
                    continue

                # nested exceptions suck, but we're going to do them anyway.
                try:
                    restp = larr[3]
                    segid = larr[4]
                    resid = int(larr[5])
                except:
                    continue
                if not (restp,segid,resid) in protonizable:
                    protonizable.append((restp,segid,resid))
            fp.close()
        except:
            pass
    return protonizable

def protonate(file):
    return render_to_response('html/protonate.html')    
