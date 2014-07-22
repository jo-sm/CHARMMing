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
import openbabel
from httplib import HTTPConnection
from django import forms
from django.db import transaction
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from django.contrib import messages
from django.contrib.messages import get_messages
#from pdbinfo.models import PDBFile, PDBFileForm, ParseException, energyParams
#from minimization.models import minimizeParams
#from minimization.views import append_tpl
#from dynamics.models import mdParams, ldParams, sgldParams
#from solvation.models import solvationParams
from account.models import *
from account.views import checkPermissions
#from dynamics.views import combinePDBsForMovie
#from normalmodes.views import combineNmaPDBsForMovie
#from normalmodes.aux import getNormalModeMovieNum
#from normalmodes.models import nmodeParams
#from apbs.models import redoxParams
#from pdbinfo.qmmm import makeQChem, makeQChem_tpl, handleLinkAtoms, writeQMheader
from django.contrib.auth.models import User
#from django.core import validators 
#from django.core.mail import mail_admins
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from account.views import isUserTrustworthy
#from pdbinfo.editscripts import generateHTMLScriptEdit
#from pdbinfo.aux import checkNterPatch
#dd stuff
from dd_substrate.models import ligands, ligand_sets, ligand_sets_ligands, poses
from dd_infrastructure.models import files,files_objects,file_types, sources
from dd_target.models import proteins
###
import output
import lesson1
import lesson2
import lesson3
import lesson4
import lessonaux
from structure.models import Structure, WorkingStructure, WorkingFile, Task, Segment
# import all created lessons by importing lesson_config
# also there is a dictionary called 'file_type' in lesson_config.py specififying the file type of files uploaded by the lessons  
from lesson_config import *
import os
import re
import copy
import datetime
import time
import mimetypes
import time
import stat
import string
import random
import glob
import sys, traceback
import commands
import charmming_config
import settings
from dd_substrate import common
#import common
import MySQLdb
import MySQLdb.cursors

#@transaction.commit_manually
def setLigands(request,ligandset_id,addedids,removedids):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    pclog=open("/tmp/sllog.log",'w')
    pclog.write("addedids in:" + addedids +'\n')
    public_user=User.objects.get(username='public')
    if (ligandset_id!='0') and (ligandset_id!=''):
        if (addedids!="") and (addedids!='0'):
            addedids_list=[int(n) for n in addedids.split(',')]
            pclog.write("addedids:"+str(addedids_list)+'\n')
            added_ligands=ligands.objects.filter(id__in=addedids_list)
            pclog.write("addedc:" + str(added_ligands)+'\n')

            for ligand in added_ligands:
                set_ligand=ligand_sets_ligands\
                                    (owner_id=request.user.id,ligands_set_id=ligandset_id,
                                    ligands_id=ligand.id)
                set_ligand.save()
                pclog.write("pc saved:" + str(set_ligand) +'\n')
            #transaction.commit()
    
        if (removedids!="") and (removedids!='0'):                      
            pclog.write("removedids: %s\n" % (removedids))
            removedids_list=[int(n) for n in removedids.split(',')]
            pclog.write("removedids_list: %s\n" % (removedids_list))
            ligand_sets_ligands.objects.filter(id__in=removedids_list).delete()
            #transaction.commit()
            
         #try:
        lig_set=ligand_sets.objects.get(id=ligandset_id)
        setligands=ligand_sets_ligands.objects.filter(owner=request.user,ligands_set=lig_set).select_related()
 
        ligand_info_list=[]

        for setligand in setligands:
            ligand_info=common.LigandInfo()
            ligand_info.id=setligand.ligands.id
            ligand_info.ligandsetid=setligand.id
            ligand_info.name=setligand.ligands.ligand_name
            ligand_info.description=setligand.ligands.description
            ligand_info_list.append(ligand_info)

        for ligand_info in ligand_info_list:
            pclog.write("id in (select file_id from dd_infrastructure_files_objects where object_table_name='dd_substrate_ligands' " \
            " and object_id=%s)\n" % (ligand_info.id))
            try:
                ligandfile=files.objects.filter().extra \
                (where=["id in (select file_id from dd_infrastructure_files_objects where object_table_name='dd_substrate_ligands' " \
                " and object_id=%s)" % (ligand_info.id)])[0]
                ligand_info.file_id=ligandfile.id
                ligand_info.file_location = ligandfile.file_location.replace("/home/schedd/","/charmming/pdbuploads/") + ligandfile.file_name
            except:
                ligand_info.file_id=0
            #for ligandfile in ligandfiles:                                                                                                                                                                              
            #ligand_info.file=ligandfile

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddsetligands.html', \
                              {'setligands':ligand_info_list,'ligandcount':len(setligands), 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

#@transaction.commit_manually
def updateAvailableSetLigands(request,selected_set_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')

    
    ligand_info_list=[]
    public_user=User.objects.get(username='public')
    log=open("/tmp/allog.log",'w')
    log.write("set id:" + selected_set_id + "\n")
    if (str(selected_set_id)=='0'): #show all ligands
    
        user_ligands = ligands.objects.filter(owner=request.user).select_related()
        
        log.write("public user id:" + str(public_user.id) + "\n")
        
        public_ligands = ligands.objects.filter(owner=public_user)
        
        
        for user_ligand in user_ligands:
            ligand_info=common.LigandInfo()
            ligand_info.id=user_ligand.id
            log.write("user ligand id:" + str(user_ligand.id) + "\n")
            ligand_info.name=user_ligand.ligand_name
            ligand_info.description=user_ligand.description
            #log.write("id in (select file_id from dd_infrastructure_files_objects where object_table_name='dd_substrate_ligands' "\
            #"and object_id=%s)" % (user_ligand.id))
            #ligand_info.file_objects_list = files_objects.objects.filter\
            #                            (owner=request.user,object_table_name="dd_substrate_ligands", \
            #                            object_id=ligand.id).select_related()
            #ligandfileobject=files_objects.objects.get(object_table_name="dd_substrate_ligands", object_id=user_ligand.id)
            #ligandfile=files.objects.filter().extra\
            #(where=["id in (select file_id from dd_infrastructure_files_objects where object_table_name='dd_substrate_ligands' "\
            #"and object_id=%s)" % (ligand.id)])[0]
            ligand_info_list.append(ligand_info)
            
        for ligand in public_ligands:
            log.write("public ligand id:" + str(ligand.id) + "\n")
            ligand_info=common.LigandInfo()
            ligand_info.id=ligand.id
            ligand_info.name=ligand.ligand_name
            ligand_info.description=ligand.description
            #ligand_info.file_objects_list = files_objects.objects.filter\
            #                            (owner=public_user,object_table_name="dd_substrate_ligands", \
            #                            object_id=ligand.id).select_related()
            #ligandfile=files.objects.filter().extra\
            #(where=["id in (select file_id from dd_infrastructure_files_objects where object_table_name='dd_substrate_ligands' "\
            #"and object_id=%s)" % (ligand.id)])[0]
            ligand_info_list.append(ligand_info)
            
    else: #show selected set ligands
        selected_set_ligands=ligands.objects.filter(owner__in=[request.user,public_user]).extra \
        (where=["id IN (select ligands_id from dd_substrate_ligand_sets_ligands where owner_id in (%s,%s) " \
        " and ligands_set_id=%s)" % (request.user.id,public_user.id,selected_set_id)])
        for ligand in selected_set_ligands:
            ligand_info=common.LigandInfo()
            ligand_info.id=ligand.id
            ligand_info.name=ligand.ligand_name
            ligand_info.description=ligand.description                                
            #ligandfile=files.objects.filter().extra\
            #(where=["id in (select file_id from dd_infrastructure_files_objects where object_table_name='dd_substrate_ligands' "\
            #"and object_id=%s)" % (ligand.id)])[0]
            ligand_info_list.append(ligand_info)
        

    for ligand_info in ligand_info_list:
        #log.write("ligand_info.id = %s\n" % (ligand_info.id))
        #log.write("id in (select file_id from dd_infrastructure_files_objects where object_table_name='dd_substrate_ligands' "\
        #" and object_id=%s) and file_name like '%%.mol2%%'\n" % (ligand_info.id))
        try:
            ligandfile=files.objects.filter().extra\
            (where=["id in (select file_id from dd_infrastructure_files_objects where object_table_name='dd_substrate_ligands' "\
            "and object_id=%s)" % (ligand_info.id)])[0]
            ligand_info.file_id=ligandfile.id
            ligand_info.file_location = ligandfile.file_location.replace("/home/schedd/","/charmming/pdbuploads/") + ligandfile.file_name
        except:
            ligand_info.file_id=0
        #for ligandfile in ligandfiles:
        
        #log.write("ligand file id: %s" % (ligandfile.id))
        #ligand_info.file_id=ligandfile.id

    log.write("list:" + str(ligand_info_list) + "\n")

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddavailablesetligands.html', {'available_ligands':ligand_info_list, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

@transaction.commit_manually    
def deleteLigandSet(request,set_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')

   
    log=open("/tmp/deleteset.log", 'w')
    log.write("deleting set: %s\n" % (set_id))
    #if 1==1:
    try:
        setligands=ligand_sets_ligands.objects.filter(ligands_set=ligand_sets.objects.get(id=set_id))
        log.write("deleting set ligands: %s\n" % (str(setligands)))
        setligands.delete()

        deleteset=ligand_sets.objects.get(id=set_id)
        log.write("deleteset: %s\n" % (str(deleteset)))
        deleteset.delete()

        transaction.commit()
    #else:
    except:
        transaction.rollback()

    return HttpResponse('Done')

def ligandSetDetails(request,ligandset_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    
    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    try:
        dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)

    except:
        pass
       
    pdlog=open("/tmp/pdlog.log",'w')
    pdlog.write("set:" + ligandset_id + "\n")
    public_user=User.objects.get(username='public')
    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    setssql="select '0' as setid, 'All Available' as setname, 'All Ligands available to the user' as setdescription, " \
            "(select count(distinct id) from dd_substrate_ligands where owner_id in (%s,%s)) as ligandcount " \
	        " union " \
	        "select ls.id as setid, ls.ligand_set_name as setname, ls.description as setdescription, " \
            "(select count(distinct ligands_id) from dd_substrate_ligand_sets_ligands where " \
            "ligands_set_id=ls.id and owner_id in (%s,%s)) as ligandcount " \
	        "from dd_substrate_ligand_sets ls where ls.owner_id in (%s,%s) and ls.ligand_set_name<>'All Available' and ls.id <> %s" \
	        % (request.user.id,public_user.id,request.user.id,public_user.id,request.user.id,public_user.id,ligandset_id)
    #setssql="select '0' as setid, gs.name as setname from vcs_gridsets gs"
    pdlog.write("setssql:" + setssql + "\n")
    cursor.execute(setssql)
    rows = cursor.fetchall()
    sets=[]

    for row in rows:
        set=common.LigandSetObj()
        set.id=row['setid']
        set.name=row['setname']
        set.ligand_count=row['ligandcount']
        sets.append(set)

    set_id=ligandset_id.replace("/","")
    
    lesson_ok, dd_ok = checkPermissions(request)    
    return render_to_response('html/ddligandsetdetails.html', \
                              {'existingsets':sets, \
                              'set_id':set_id, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
"""                              
def viewUpdateLigandSetInfo(request,ligandsetid,action):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    pilog=open("/tmp/lslog.log",'w')
    pilog.write("ls\n")
    username = request.user.username
    owner = User.objects.get(username=username)

    lesson_ok, dd_ok = checkPermissions(request)    

    #update or add is taking place
    if ((action=='update') or (action=='addnew')):

        if (len(request.POST['name'].strip())==0): 
            message="Operation Failed!<br>Ligand Set Name Cannot be Blank"
            if (request.POST['action']=='update'):
                ligandset = ligand_sets.objects.get(owner=request.user,id=ligandset)
                return render_to_response('html/ddligandsetinfo.html', {'ligandset':ligandset, 'message':message, 'messagecolor':'Red', 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
            else:
                return render_to_response('html/ddnewligandsetinfo.html', {'description':request.POST['description'], 'name':request.POST['name'], 'message':message, 'messagecolor':'Red', 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
        else:
            if (action=='update'):
                ligandset = ligand_sets.objects.get(owner=request.user,id=ligandsetid)
                ligandset.ligand_set_name=request.POST['name'].strip()
                ligandset.description=request.POST['description'].strip()
                ligandset.save()
                ligandset_id=ligandsetid
                message="Ligand Set Info was successfully updated"
                
            elif (action=='addnew'):

                if common.LigandSetExists(request.POST['name'],request.user)==True:
                    return render_to_response('html/ddnewligandsetinfo.html', {'description':request.POST['description'], 'name':request.POST['name'], 'message':'Operation Failed!<br>Duplicate Ligand Set Found.', 'messagecolor':'Red', 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

                ligandset=ligand_sets()
                ligandset.owner=owner
                #ligandset.date_created=datetime.datetime.now().isoformat(' ')
                ligandset.ligand_set_name=request.POST['name'].strip()
                ligandset.description=request.POST['description'].strip()
                ligandset.save()
                #ligandset_id = ligand_set.id
                message="Ligand Set was successfully created"
                
        return render_to_response('html/ddligandsetinfo.html', {'ligandset':ligandset, 'message':message, 'messagecolor':'Green', 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

    #no action, just display blank or populated form
    else: 
        if (ligandsetid!='0'):
            pilog.write("ligandsetid:" + ligandsetid + "\n")
            ligandset = ligand_sets.objects.get(owner=request.user,id=ligandsetid)
            pilog.write("ligandset name:" + str(ligandset.ligand_set_name) + "\n")
            return render_to_response('html/ddligandsetinfo.html', {'ligandset':ligandset, 'message':'', 'messagecolor':'Red', 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
        else:
            return render_to_response('html/ddnewligandsetinfo.html', {'description':'', 'message':'', 'messagecolor':'Red', 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
"""
def ligandSetsAddNew(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')

    
    lesson_ok, dd_ok = checkPermissions(request)    
    return render_to_response('html/ddligandsetsaddnew.html', {'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
    
#def viewLigandSetsContainer(request):
#    if not request.user.is_authenticated():
#        return render_to_response('html/loggedout.html')
#    return render_to_response('html/ddviewligandsetscontainer.html')


def viewLigandSets(request):
    log=open("/tmp/ligandsets.log",'w')
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    try:
        dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)

    except:
       pass
    public_user=User.objects.get(username='public')
    log.write("user:" + str(request.user.id))
    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    setssql="select ls.id as setid, ls.ligand_set_name as setname, ls.description as setdescription, " \
            "(select count(distinct ligands_id) from dd_substrate_ligand_sets_ligands where " \
            "ligands_set_id=ls.id and owner_id in (%s)) as ligandcount " \
	        "from dd_substrate_ligand_sets ls where ls.owner_id in (%s) and ls.ligand_set_name<>'All Available'" \
	        % (request.user.id,request.user.id)
    log.write(setssql)
    cursor.execute(setssql)
    rows = cursor.fetchall()
    sets=[]
    for row in rows:
        set=common.LigandSetObj()
        set.id=row['setid']
        set.name=row['setname']
        set.description=row['setdescription']
        set.grid_count=row['ligandcount']
        sets.append(set)     


    cursor.close()
    dba.close()
    log.close()

    lesson_ok, dd_ok = checkPermissions(request)    
    return render_to_response('html/ddviewligandsets.html', {'ligandsets':sets, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
    
#def WNPRefreshLigands(request):
#    
#    log=open("/tmp/ligand_batch_upload.log",'w')
#    location=charmming_config.user_dd_ligands_home + '/public/'
#    wnpscript=open(location + "/wnpscript.inp",'w')
#    wnpscript.write("cd %s\n" % (location))
#    #os.system("touch %swnp_fix" % (location))
#    wnpscript.write("touch %swnp_fix\n" % (location))
#    wnpscript.write('for f in %s*.mol2\n' % (location))
#    wnpscript.write('do\n')
#    wnpscript.write('    echo "read mol2 ${f%.*}" >> wnp_fix\n')
#    #wnpscript.write('    echo "modi hydrogen *" >> wnp_fix_partial_charges\n')
#    wnpscript.write('    echo "modi " >> wnp_fix\n')
#    wnpscript.write('    echo "atom lig current *" >> wnp_fix\n')
#    wnpscript.write('    echo "atom charge auto *" >> wnp_fix\n')
#    wnpscript.write('    echo "atom type auto charmm *" >> wnp_fix\n')
#    wnpscript.write('    echo "done" >> wnp_fix\n')
#    wnpscript.write('    echo "atom name frag -done -auto seq" >> wnp_fix\n')
#    wnpscript.write('    echo "atom q * -done mpeoe" >> wnp_fix\n')
#    wnpscript.write('    echo "done" >> wnp_fix\n')
#    """
#    submitscript.write('    echo "modify atom q **** -done MPEOE" >> wnp_fix\n')
#    submitscript.write('    echo "done" >> wnp_fix\n')
#    """ 
#    wnpscript.write('    bf=${f##*/}\n')
#    wnpscript.write('    molecule=`echo | awk "NR==2 {print;exit}" $f`\n')
#    wnpscript.write('    echo "write mol2 ${f%.*}  $molecule" >> wnp_fix\n')
#    wnpscript.write('    echo "done" >> wnp_fix\n')
#    wnpscript.write('    echo "delete molecules /$molecule " >> wnp_fix\n')
#    #wnpscript.write('    %s host none < wnp_fix\n' % (charmming_config.witnotp_exe))
#    #wnpscript.write('    echo `mv ${f%.*}.mol ${f%.*}.mol2`\n')
#    #wnpscript.write("    echo `sed -i 's/N.?/NP /' ${f%.*}.mol2`\n")
#    wnpscript.write('done\n')
#    wnpscript.write('    %s host none < wnp_fix\n' % (charmming_config.witnotp_exe))
#    os.system("cd %s" % (location))
#    os.system("sh wnpscript.inp")
#    """
#    u = User.objects.get(username="public")
#    for file in glob.glob(os.path.join(location, '*.mol2')):
#        file_obj=open(file,'r')
#        mol_title=getMoleculeTitle(file_obj)
#        newligand = ligands()
#        newfile = files()
#        newfileobject= files_objects()
#        try: 
#            newligand=ligands.objects.filter(owner=u,name=mol_title)[0]
#        except:
#            newligand.owner=u
#            newligand.ligand_name=mol_title
#            newligand.ligand_owner_index=str(common.getNewLigandOwnerIndex(u))
#            newligand.description="Preloaded ZINC lead-like CHARMm compound"
#            source=sources.objects.filter(source_name='Zinc Database')[0]
#            newligand.source=source
#            log.write("saving ligand with owner id %s and name %s\n" % (newligand.ligand_owner_index,newligand.ligand_name))
#            newligand.save()
#       
#            new_location=charmming_config.user_dd_ligands_home + '/' + 'public' + '/' + 'ligand_' + str(newligand.ligand_owner_index) + '/'
#            filename='ligand_' + str(newligand.ligand_owner_index) + '_1.mol2'
#            ligand_file_type=file_types.objects.filter(file_type_name="Ligand Structure File")[0]
#            newfile.owner=u
#            newfile.file_name=filename
#            newfile.file_location=new_location
#            newfile.file_type=ligand_file_type
#            newfile.description="Structure File for Ligand %s, %s" % (newligand.ligand_name,newligand.description)
#            log.write("saving file with location %s and name %s\n" % (newfile.file_location,newfile.file_name))
#            newfile.save()
#            
#            newfileobject.owner=u
#            newfileobject.file=newfile
#            newfileobject.object_table_name='dd_substrate_ligands'
#            newfileobject.object_id=newligand.id
#            log.write("saving new fileobject\n")
#            newfileobject.save()
#            
#            os.system("mkdir " + new_location + "\n")
#            os.system("cp %s %s" % (file,new_location + newfile.file_name))
#            log.write("mkdir " + new_location + "\n")
#            log.write("cp %s %s\n" % (file,new_location + newfile.file_name))
#    """    
#    log.close()
#    return HttpResponse('Done')

def getMoleculeTitle(file):
    while 1:
        line = file.readline()
        if not line:
            break
        if line.strip() == "@<TRIPOS>MOLECULE":
            titleline=file.readline()
            title=titleline.strip()
            break
    return title
    
#@transaction.commit_manually
def newLigandUpload(request, template="html/ddligandfileupload.html"):

    command_list=[]
    uploadligandlog=open("/tmp/ligandupload.log",'w')
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')

    try:
        struct = Structure.objects.get(owner=request.user,selected='y')
    except:
        return HttpResponse("Please submit a structure first.")

    try:
        ws = WorkingStructure.objects.get(structure=struct,selected='y')
    except:
        return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    username=request.user.username
    u = User.objects.get(username=username)
    #transaction.commit()
    newligand = ligands()
    #newfile = files()
    #newfileobject= files_objects()
    newligandpose=poses()
    newfileposeobject=files_objects()
    ligand_file_type=file_types.objects.filter(file_type_name="Ligand Structure File")[0] 
    psf_file_type=file_types.objects.get(file_type_name="CHARMM Structure File")
    str_file_type=file_types.objects.get(file_type_name="CHARMM Param Stream File")
    #transaction.commit()
    user_source=sources.objects.filter(source_name='User')[0]
    #transaction.commit()
    form=ligandFileForm()
    try:
        request.FILES['ligand_file'].name
        ligand_uploaded_by_user = 1
        uploadligandlog.write("user uploaded file\n")
    except:
        ligand_uploaded_by_user = 0

    if ligand_uploaded_by_user==1:
        if request.POST['pdbcode']<>"":
            sourceproteins=proteins.objects.filter(owner=u,pdb_code=request.POST['pdbcode'])
        #transaction.commit() 
        #this commented block of code is needed when checking for existing ligand is coded in
        #if 1==1:
        #try: 
        #    newligand=ligands.objects.filter(owner=u,id=request.POST[''])[0]
        #    transaction.commit()
        ####### Atomtyping check
        tempname='dd_user_ligand'
        os.system("rm /tmp/%s.mol2" % (tempname))
        destination = open("/tmp/%s.mol2" % (tempname),'w')
        for fchunk in request.FILES['ligand_file'].chunks():
            destination.write(fchunk)
        destination.close()
        
        if request.POST['tpmeth'] == 'autogen':
            logfp = open('/tmp/dd_autogen.txt', 'w')
            success = False

            for tpmeth in charmming_config.dd_toppar_generators.split(','):
                logfp.write('trying method: %s\n' % tpmeth)
                if tpmeth == 'genrtf':
                    rval = -1#self.makeGenRTF(bhResList)
                elif tpmeth == 'antechamber':
                    rval = -1#common.makeAntechamber(tempname)
                elif tpmeth == 'cgenff':
                    rval = common.makeCGenFF(ws,tempname)
                    if rval ==0:
                        os.system("cp %s/build_%s.inp /tmp" % (charmming_config.dd_scripts_home,tempname))
                        os.chdir("/tmp")
                        os.system("%s < build_%s.inp > build_%s.out" % (charmming_config.charmm_exe,tempname,tempname))
                        if "NORMAL TERMINATION BY NORMAL STOP" in open("/tmp/build_%s.out" % (tempname)).read():
                            rval=0
                        else:
                            rval=-1
                elif tpmeth == 'match':
                    rval = common.makeMatch(tempname)
                    if rval ==0:
                        os.system("cp %s/build_%s.inp /tmp" % (charmming_config.dd_scripts_home,tempname))
                        os.chdir("/tmp")
                        os.system("%s < build_%s.inp > build_%s.out" % (charmming_config.charmm_exe,tempname,tempname))
                        if "NORMAL TERMINATION BY NORMAL STOP" in open("/tmp/build_%s.out" % (tempname)).read():
                            rval=0
                        else:
                            rval=-1

                logfp.write('got rval = %d\n' % rval)
                if rval == 0:
                    success = True
                    break

            logfp.close()
           
            if not success:
                messages.error(request,"Could not build topology/parameters for your uploaded file using MATCH and CGenFF. Please verify that you uploaded a properly-formatted .MOL2 file.")
                return render_to_response('html/ddligandfileupload.html', {'form':form, 'super_user':request.user.is_superuser,'messages':get_messages(request)})
#                raise AssertionError('Unable to build topology/parameters/structure file')
          
            #if common.checkDAIM(tempname)!=0:
            #    raise AssertionError('With the current decomposition settings the uploaded ligand cannot be prepped for docking')
            #    success = False
                
        ###### Atomtyping check end
        if success:
        
        #except:
        
        #try:

            newligand.owner=u
            newligand.ligand_name=request.POST['ligandname']
            newligand.ligand_owner_index=str(common.getNewLigandOwnerIndex(u))
            newligand.description=request.POST['liganddescription']
            newligand.source=user_source

            #newligand.id=str(996)
            newligand.save()
            uploadligandlog.write("saving new ligand: %s\n" % (str(newligand.id) + " " + newligand.ligand_name))

                #transaction.commit()
            #newligand = ligands.objects.filter(owner=u).order_by("-id")[0]
        #if 1==1:
        #try:
            os.system("mkdir " + charmming_config.user_dd_ligands_home + '/' + username + '/')
            location = charmming_config.user_dd_ligands_home + '/' + username + '/' + 'ligand_' + str(newligand.ligand_owner_index) + '/'
            filename='ligand_' + str(newligand.ligand_owner_index) + '_1.mol2'
            uploadligandlog.write("mkdir " + location)
            os.system("mkdir " + location) 
            destination = open(location + filename,'w')
            for fchunk in request.FILES['ligand_file'].chunks():
               destination.write(fchunk)
            destination.close()
            newfile=files()
            newfile.owner=u
            newfile.file_name=filename
            newfile.file_location=location
            newfile.file_type=ligand_file_type
            newfile.description="Structure File for Ligand %s, %s" % (newligand.ligand_name,newligand.description)
            #newfile.id=str(999)
            newfile.save()
            uploadligandlog.write("saving new ligand: %s\n" % (str(newligand.id) + " " + newligand.ligand_name))
            #newfile = files.objects.filter(owner=u).order_by("-id")[0]
            
            
            newfileobject=files_objects()
            newfileobject.owner=u
            newfileobject.file=newfile
            newfileobject.object_table_name='dd_substrate_ligands'
            newfileobject.object_id=newligand.id
            newfileobject.save()
            uploadligandlog.write("saving new file object: %s\n" % (newfileobject.object_table_name + " " + str(newfileobject.object_id)))

            psffilename='ligand_' + str(newligand.ligand_owner_index) + '_1.psf'
            os.system("cp /tmp/dd_user_ligand.psf %s%s" % (location,psffilename))
            newpsffile=files()
            newpsffile.owner=u
            newpsffile.file_name=psffilename
            newpsffile.file_location=location
            newpsffile.file_type=psf_file_type
            newpsffile.description="CHARMM Structure File for Ligand %s, %s" % (newligand.ligand_name,newligand.description)
            #newfile.id=str(999)                                                                                                                              
            newpsffile.save()
            uploadligandlog.write("saving new psf file for ligand: %s\n" % (str(newligand.id) + " " + newligand.ligand_name))
            #newfile = files.objects.filter(owner=u).order_by("-id")[0]                                                                                       
            
            newpsffileobject=files_objects()
            newpsffileobject.owner=u
            newpsffileobject.file=newpsffile
            newpsffileobject.object_table_name='dd_substrate_ligands'
            newpsffileobject.object_id=newligand.id
            newpsffileobject.save()
            uploadligandlog.write("saving new file object: %s\n" % (newpsffileobject.object_table_name + " " + str(newpsffileobject.object_id)))

            strfilename='ligand_' + str(newligand.ligand_owner_index) + '_1.str'
            os.system("cp /tmp/dd_user_ligand.str %s%s"% (location,strfilename))
            os.system("cp /tmp/dd_user_ligand.rtf %s%s"% (location,strfilename.replace("str","rtf")))
            
            newstrfile=files()
            newstrfile.owner=u
            newstrfile.file_name=strfilename
            newstrfile.file_location=location
            newstrfile.file_type=ligand_file_type
            newstrfile.description="CHARMM Param Stream File for Ligand %s, %s" % (newligand.ligand_name,newligand.description)
            #newfile.id=str(999)                                                                                                                              
            newstrfile.save()
            uploadligandlog.write("saving new stream file for ligand: %s\n" % (str(newligand.id) + " " + newligand.ligand_name))
            #newfile = files.objects.filter(owner=u).order_by("-id")[0]                                                                                       
            
            newstrfileobject=files_objects()
            newstrfileobject.owner=u
            newstrfileobject.file=newstrfile
            newstrfileobject.object_table_name='dd_substrate_ligands'
            newstrfileobject.object_id=newligand.id
            newstrfileobject.save()
            uploadligandlog.write("saving new stream file object: %s\n" % (newstrfileobject.object_table_name + " " + str(newstrfileobject.object_id)))


            if request.POST['pdbcode']<>"":
                pose_source=sources()

                if len(sourceproteins)!=0:
                    pose_source.source_object_table_name="dd_target_proteins"
                    pose_source.source_object_id=sourceproteins[0].id

                pose_source.description="Protein Crystal Structure " + request.POST['pdbcode']
                pose_source.source_name="PDB " + request.POST['pdbcode']
                uploadligandlog.write("saving new source: %s\n" % (pose_source.source_name + ", " + str(pose_source.source_object_id)))
                #pose_source.id=str(997)
                pose_source.save()
            else:
                pose_source=user_source

            newligandpose.owner=u
            newligandpose.pose_object_table_name='dd_substrate_ligands'
            newligandpose.pose_object_id=newligand.id
            newligandpose.name='Default Pose'
            newligandpose.description='Initially uploaded pose'
            newligandpose.source_id=pose_source.id
            #newligandpose.id=str(998)
            newligandpose.save()
            uploadligandlog.write("saving new ligand pose: %s\n" % (str(newligandpose.pose_object_id) + " " + str(newligandpose.source_id)))

            newfileposeobject.owner=u
            newfileposeobject.file=newfile
            newfileposeobject.object_table_name='dd_substrate_poses'
            newfileposeobject.object_id=newligandpose.id
            newfileposeobject.save()
            uploadligandlog.write("saving new file pose object: %s\n" % (newfileposeobject.file.file_name + " " + str(newfileposeobject.object_id)))
            
            #transaction.commit()
        #except:
        #    transaction.rollback()
        
        #transaction.commit()
        
        uploadligandlog.close()
	return HttpResponseRedirect("/charmming/dd_infrastructure/dsfform/")

    lesson_ok, dd_ok = checkPermissions(request)        
    return render_to_response('html/ddligandfileupload.html', {'form':form, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok, 'super_user':request.user.is_superuser})

#def ligandFileUploadForm(request):
#    if not request.user.is_authenticated():
#        return render_to_response('html/loggedout.html')
#    
#    form = ligandFileForm()
#    return render_to_response('html/ddligandfileuploadform.html', {'form':form})

def ligandDropdown(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')

    
    ligandlist=ligands.objects.filter(owner=request.user).order_by("ligand_name")

    lesson_ok, dd_ok = checkPermissions(request)    
    return render_to_response('html/ddliganddropdown.html', {'ligandlist':ligandlist, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

class ligandFileForm(forms.Form):
    ligand_file = forms.FileField()
    
    
# main routine
def RefreshLigands(request):

    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')


    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    log=open("/tmp/ligand_batch_upload.log",'w')
    location=charmming_config.user_dd_ligands_home + '/public/zinctest3/'
    log.write("starting ligand upload\n")
    os.system("cd %s\n" % (location))
    os.environ["MATCH"]="/usr/local/charmming/MATCH_RELEASE/MATCH"
    os.environ["PerlChemistry"]="/usr/local/charmming/MATCH_RELEASE/PerlChemistry"
    u = User.objects.get(username='public')
    log.write("environ match %s\n" % (os.environ["MATCH"]))
    log.write("newindex is %s: \n" % (str(common.getNewLigandOwnerIndex(u))))
    count=0
    newligands=[]
    for file in glob.glob(os.path.join(location, '*.mol2')):
        os.system("/usr/local/charmming/MATCH_RELEASE/MATCH/scripts/MATCH.pl %s > %smatchresult\n" % (file,location))
        log.write("/usr/local/charmming/MATCH_RELEASE/MATCH/scripts/MATCH.pl %s > %smatchresult\n" % (file,location))
        #time.sleep(10)
        while os.path.exists(location+"matchresult") == False: time.sleep(3)
        if "Success!" in open(location+'matchresult').read():
            log.write("match success for ligand %s\n" % (file))
            file_obj=open(file,'r')
            mol_title=getMoleculeTitle(file_obj)
            newligand = ligands()
            newfile = files()
            newfileobject= files_objects()
            try: 
                newligands=ligands.objects.filter(owner=u,ligand_name=mol_title)
                log.write("name %s exists\n" % (mol_title))
            except:
                pass
            if len(newligands)==0:
                log.write("name %s doesn't exist\n" % (mol_title))
                newligand.owner=u
                newligand.ligand_name=mol_title
                newligand.ligand_owner_index=str(common.getNewLigandOwnerIndex(u))
                newligand.description="System preloaded drug-like compound"
                source=sources.objects.get(source_name='Zinc Database')
                newligand.source=source
                log.write("saving ligand with owner id %s and name %s\n" % (newligand.ligand_owner_index,newligand.ligand_name))
                newligand.save()

                new_location=charmming_config.user_dd_ligands_home + '/' + 'public' + '/' + 'ligand_' + str(newligand.ligand_owner_index) + '/'
                filename='ligand_' + str(newligand.ligand_owner_index) + '_1.mol2'
                ligand_file_type=file_types.objects.get(file_type_name="Ligand Structure File")
                newfile.owner=u
                newfile.file_name=filename
                newfile.file_location=new_location
                newfile.file_type=ligand_file_type
                newfile.description="Structure File for Ligand %s, %s" % (newligand.ligand_name,newligand.description)
                log.write("saving file with location %s and name %s\n" % (newfile.file_location,newfile.file_name))
                newfile.save()

                newfileobject.owner=u
                newfileobject.file=newfile
                newfileobject.object_table_name='dd_substrate_ligands'
                newfileobject.object_id=newligand.id
                log.write("saving new fileobject\n")
                newfileobject.save()

                os.system("mkdir " + new_location)
                os.system("cp %s %s%s" % (file,new_location,filename))
                log.write("mkdir %s\n" % (new_location))
                log.write("cp %s %s%s\n" % (file,new_location,filename))
            else:
                log.write("ligand %s exists\n" % (mol_title))    
        else:
            log.write("match didn't work on ligand %s\n" % (file))
    
        count=count+1
        if count>1000:
            break

    return HttpResponse('Done')        
    

def RefreshLigands_back(request):

    log=open("/tmp/ligand_batch_upload.log",'w')
    location=charmming_config.user_dd_ligands_home + '/public/zinc_initial_test_compounds/'
    u = User.objects.get(username='public')
    log.write("newindex is %s: \n" % (str(common.getNewLigandOwnerIndex(u))))
    for file in glob.glob(os.path.join(location, '*.mol2')):
        file_obj=open(file,'r')
        mol_title=getMoleculeTitle(file_obj)
        newligand = ligands()
        newfile = files()
        newfileobject= files_objects()
        try: 
            newligand=ligands.objects.get(owner=u,name=mol_title)
        except:
            newligand.owner=u
            newligand.ligand_name=mol_title
            newligand.ligand_owner_index=str(common.getNewLigandOwnerIndex(u))
            newligand.description="System preloaded drug-like compound"
            source=sources.objects.get(source_name='Zinc Database')
            newligand.source=source
            log.write("saving ligand with owner id %s and name %s\n" % (newligand.ligand_owner_index,newligand.ligand_name))
            newligand.save()
       
            new_location=charmming_config.user_dd_ligands_home + '/' + 'public' + '/' + 'ligand_' + str(newligand.ligand_owner_index) + '/'
            filename='ligand_' + str(newligand.ligand_owner_index) + '_' + str(newligand.ligand_owner_index)+'.mol2'
            ligand_file_type=file_types.objects.get(file_type_name="Ligand Structure File")
            newfile.owner=u
            newfile.file_name=filename
            newfile.file_location=new_location
            newfile.file_type=ligand_file_type
            newfile.description="Structure File for Ligand %s, %s" % (newligand.ligand_name,newligand.description)
            log.write("saving file with location %s and name %s\n" % (newfile.file_location,newfile.file_name))
            newfile.save()
            
            newfileobject.owner=u
            newfileobject.file=newfile
            newfileobject.object_table_name='dd_substrate_ligands'
            newfileobject.object_id=newligand.id
            log.write("saving new fileobject\n")
            newfileobject.save()
            
            os.system("mkdir " + new_location)
            os.system("cp %s %s%s" % (file,new_location,filename))
            log.write("mkdir %s\n" % (new_location))
            log.write("cp %s %s%s\n" % (file,new_location,filename))
    
    return HttpResponse('Done')        
    
def RefreshMaybhitLigands(request):

    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')


    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    log=open("/tmp/maybhit_batch_upload.log",'w')
    location=charmming_config.user_dd_ligands_home + '/public/for_charmming/'
    log.write("starting ligand upload\n")
    os.system("cd %s\n" % (location))
    u = User.objects.get(username='public')
    log.write("newindex is %s: \n" % (str(common.getNewLigandOwnerIndex(u))))
    count=0
    newligands=[]
    try:
        maybhit_set=ligand_sets.objects.get(ligand_set_name="Maybridge Diversity Set", owner=u)
        log.write("set found\n")
    except:
        maybhit_set=ligand_sets()
        maybhit_set.ligand_set_name="Maybridge Diversity Set"
        maybhit_set.description="A subset of the Maybridge HitFinder set limited to molecules that had met docking by decomposition protocol requirements and limitations"
        maybhit_set.owner=u
        log.write("saving new set with name %s\n" % (maybhit_set.name))
        maybhit_set.save()
    log.write("############\n")

    for file in glob.glob(os.path.join(location, 'ZINC????????.mol2')):
        file_obj=open(file,'r')
        basename=os.path.basename(file)
        mol_title=os.path.splitext(basename)[0]
        chemical_name=getChemicalName(file_obj)
        str_file=glob.glob(os.path.join(location, mol_title + ".str"))[0]
        psf_file=glob.glob(os.path.join(location, mol_title + ".psf"))[0]
        newligand = ligands()
        
        
        try: 
            newligands=ligands.objects.filter(owner=u,ligand_name=mol_title)
            
        except:
            pass

        if len(newligands)==0:
            log.write("name %s doesn't exis...t proceeding with new ligand\n" % (mol_title))
            newligand.owner=u
            newligand.ligand_name=mol_title
            newligand.ligand_owner_index=str(common.getNewLigandOwnerIndex(u))
            
            if chemical_name=="":
                newligand.description="System preloaded drug-like compound"
            else:
                newligand.description=chemical_name
            
            source=sources.objects.get(source_name='Zinc Database')
            newligand.source=source
            log.write("saving ligand with owner id %s and name %s and description %s\n" % (newligand.ligand_owner_index,newligand.ligand_name,newligand.description))
            newligand.save()

            new_location=charmming_config.user_dd_ligands_home + '/' + 'public' + '/' + 'ligand_' + str(newligand.ligand_owner_index) + '/'
            filename='ligand_' + str(newligand.ligand_owner_index) + '_1.mol2'
            str_filename = 'ligand_' + str(newligand.ligand_owner_index) + '_1.str'
            psf_filename = 'ligand_' + str(newligand.ligand_owner_index) + '_1.psf'
            #min_filename = 'ligand_' + str(newligand.ligand_owner_index) + '_2.mol2'
            
            os.system("mkdir " + new_location)
            os.system("chmod g+w " + new_location)            
            log.write("mkdir %s\n" % (new_location))            

            ligand_file_type=file_types.objects.get(file_type_name="Ligand Structure File")
            newfile=files()
            newfile.owner=u
            newfile.file_name=filename
            newfile.file_location=new_location
            newfile.file_type=ligand_file_type
            newfile.description="Structure File for Ligand %s" % (newligand.ligand_name)
            log.write("saving file with location %s and name %s\n" % (newfile.file_location,newfile.file_name))
            newfile.save()
            os.system("cp %s %s%s" % (file,new_location,filename))
            log.write("cp %s %s%s\n" % (file,new_location,filename))

            str_file_type=file_types.objects.get(file_type_name="CHARMM Param Stream File")
            str_newfile=files()
            str_newfile.owner=u
            str_newfile.file_name=str_filename
            str_newfile.file_location=new_location
            str_newfile.file_type=str_file_type
            str_newfile.description="Stream File for Ligand %s" % (newligand.ligand_name)
            str_newfile.save()            
            log.write("saving str file with location %s and name %s and description %s\n" % (str_newfile.file_location,str_newfile.file_name,str_newfile.description))
            os.system("cp %s %s%s" % (str_file,new_location,str_filename))     
            log.write("cp %s %s%s\n" % (str_file,new_location,str_filename))

            psf_file_type=file_types.objects.get(file_type_name="CHARMM Structure File")
            psf_newfile=files()
            psf_newfile.owner=u
            psf_newfile.file_name=psf_filename
            psf_newfile.file_location=new_location
            psf_newfile.file_type=psf_file_type
            psf_newfile.description="CHARMM Structure File for Ligand %s" % (newligand.ligand_name)
            log.write("saving psf file with location %s and name %s and description %s\n" % (psf_newfile.file_location,psf_newfile.file_name,psf_newfile.description))
            psf_newfile.save()
            os.system("cp %s %s%s" % (psf_file,new_location,psf_filename))
            log.write("cp %s %s%s\n" % (psf_file,new_location,psf_filename))
            
            newfileobject= files_objects()
            newfileobject.owner=u
            newfileobject.file=newfile
            newfileobject.object_table_name='dd_substrate_ligands'
            newfileobject.object_id=newligand.id            
            newfileobject.save()
            log.write("saving new fileobject with fileid: %s and objectid: %s of table %s\n" % (newfileobject.file_id, newfileobject.object_id, newfileobject.object_table_name))
 
            newfileobject= files_objects()
            newfileobject.owner=u
            newfileobject.file=str_newfile
            newfileobject.object_table_name='dd_substrate_ligands'
            newfileobject.object_id=newligand.id
            newfileobject.save()
            log.write("saving new str fileobject with fileid: %s and objectid: %s of table %s\n" % (newfileobject.file_id, newfileobject.object_id, newfileobject.object_table_name))
            
            newfileobject= files_objects()
            newfileobject.owner=u
            newfileobject.file=psf_newfile
            newfileobject.object_table_name='dd_substrate_ligands'
            newfileobject.object_id=newligand.id
            newfileobject.save() 
            log.write("saving new psf fileobject with fileid: %s and objectid: %s of table %s\n" % (newfileobject.file_id, newfileobject.object_id, newfileobject.object_table_name))
 
            ######include ligand in the set
            set_ligand=ligand_sets_ligands()
            set_ligand.ligands_set=maybhit_set
            set_ligand.ligands=newligand
            set_ligand.owner=u
            set_ligand.save()
            log.write("ligand %s is included in the set %s\n" % (set_ligand.ligands.ligand_name,set_ligand.ligands_set.ligand_set_name))


            ###end include in the set 
        else:
            log.write("ligand %s exists... moving on to the next molecule\n" % (mol_title))    
    
        log.write("########################\n")
        count=count+1
        #if count>19:
        #    break

    return HttpResponse('Done')

def getMoleculeTitle(file):
    while 1:
        line = file.readline()
        if not line:
            break
        if line.strip() == "@<TRIPOS>MOLECULE":
            titleline=file.readline()
            title=titleline.strip()
            break
    return title

def getChemicalName(file):
    chem_name=""
    while 1:
        line = file.readline()
        if not line:
            break
        if line.strip() == "USER_CHARGES":
            chem_line=file.readline()
            chem_name=chem_line.strip()
            break
    return chem_name

def viewligands(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    try:
        dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)

    except:
       pass

    public_user=User.objects.get(username='public')
    #log.write("user:" + str(request.user.id))
    cursor = dba.cursor(MySQLdb.cursors.DictCursor)
    ligandssql="select l.id,l.ligand_owner_index, l.ligand_name,l.description, s.source_name " \
                "from dd_substrate_ligands l left join dd_infrastructure_sources s on l.source_id=s.id " \
                "where l.owner_id in (%s) order by l.ligand_name" \
                % (request.user.id)
    #log.write(setssql)
    cursor.execute(ligandssql)
    ligands = cursor.fetchall()

    lesson_ok, dd_ok = checkPermissions(request)       
    return render_to_response('html/ddviewligands.html', {'ligands': ligands, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

@transaction.commit_manually
def deleteligand(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    public_user=User.objects.get(username='public')
    transaction.commit()
    deleteliglog=open("/tmp/deletelig.log",'w')
    deleteliglog.write("starting deletion\n")
    #make sure request data isn't malicious
    #input.checkRequestData(request)
    if request.POST.has_key('id'):
        delete_id = request.POST['id']
    else:
        #return HttpResponse('Bad')
        delete_id="1"

    deleteAll = 0
    #If user wants to remove all files, visit /deletefile/all_files/
     
    if(delete_id == "all_ligands"):
        deleteAll = 1
        #total_files = structure.models.Structure.objects.filter(owner=request.user)
    else:
        try:
            #file = structure.models.Structure.objects.filter(owner=request.user,name=delete_filename)[0]
            #total_files = [file]
            deleteliglog.write("deleting ligand with id %s\n" % (delete_id))
            ligand=ligands.objects.get(owner__in=[request.user],id=delete_id)
            deleteliglog.write("found ligand: %s, %s\n" % (ligand.id,ligand.ligand_name))
            try:
                ligand_files_objects=files_objects.objects.filter(owner__in=[request.user],object_table_name='dd_substrate_ligands',object_id=delete_id)
                for ligand_file_object in ligand_files_objects:
                    deleteliglog.write("    found ligand file object: %s\n" % (ligand_file_object.id))
                    try:
                        ligand_files=files.objects.filter(owner__in=[request.user],id=ligand_file_object.file_id)
                        for ligand_file in ligand_files:
                            deleteliglog.write("        found ligand file: %s, %s\n" % (ligand_file.id, ligand_file.file_name))
                            try:
                                ligand_associations=files_objects.objects.filter(owner__in=[request.user], file=ligand_file).exclude(id=ligand_file_object.id)
                                for ligand_association in ligand_associations:
                                    deleteliglog.write("            found ligand association: %s, %s, %s\n" % (ligand_association.id, ligand_association.object_table_name, ligand_association.object_id))
                                ligand_associations.delete()
                            except:
                                transaction.rollback()
                                return HttpResponse('Bad')

                            ligand_files.delete()
                    except:
                        transaction.rollback()
                        return HttpResponse('Bad')

                try: 
                    ligand_poses=poses.objects.filter(pose_object_table_name="dd_substrate_ligands",pose_object_id=ligand.id)
                    for ligand_pose in ligand_poses:
                        deleteliglog.write("    found ligand pose: %s, %s\n" % (ligand_pose.id,ligand_pose.pose_name))
                        if common.DeletePose(request.user,ligand_pose.id)==False:
                            transaction.rollback()
                            return HttpResponse('Bad')
                except:
                    transaction.rollback()
                    return HttpResponse('Bad')     
                ligand_files_objects.delete()
            except:
                transaction.rollback()
                return HttpResponse('Bad')

            ligand.delete()
        except:
            transaction.rollback()
            return HttpResponse('Bad')

    return HttpResponse('Done')



def viewUpdateLigandSetInfo(request,ligandsetid,action):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')



    pilog=open("/tmp/silog.log",'w')
    pilog.write("si %s %s\n" % (ligandsetid,action))
    username = request.user.username
    owner = User.objects.get(username=username)

    #update or add is taking place
    if ((action=='update') or (action=='addnew')):
        pilog.write("name: %s, description: %s\n" % (request.POST["name"],request.POST["description"]))
        if (len(request.POST['name'].strip())==0): 
            message="Operation Failed!<br>Set Name Cannot be Blank"
            if (action=='update'):
                ligand_set = ligand_sets.objects.get(owner=request.user,id=ligandsetid)
                return render_to_response('html/ddligandsetinfo.html', {'ligand_set':ligand_set, 'message':message, 'messagecolor':'Red'})
            else:
                return render_to_response('html/ddnewligandsetinfo.html', {'description':request.POST['description'], 'name':request.POST['name'], 'message':message, 'messagecolor':'Red'})
        else:
            if (action=='update'):
                ligand_set = ligand_sets.objects.get(owner=request.user,id=ligandsetid)
                ligand_set.ligand_set_name=request.POST['name'].strip()
                ligand_set.description=request.POST['description'].strip()
                pilog.write("saving update\n")
                ligand_set.save()
                set_id=ligandsetid
                message="Ligand Set Info was successfully updated"
                
            elif (action=='addnew'):

                if common.LigandSetExists(request.POST['name'],request.user)==True:
                    return render_to_response('html/ddnewligandsetinfo.html', {'description':request.POST['description'], 'name':request.POST['name'], 'message':'Operation Failed!<br>Duplicate Project Found.', 'messagecolor':'Red'})

                ligand_set=ligand_sets()
                ligand_set.owner=owner
                ligand_set.ligand_set_name=request.POST['name'].strip()
                ligand_set.description=request.POST['description'].strip()
                ligand_set.save()
            
                added_ligand_set = ligand_sets.objects.filter(owner=request.user).order_by('-id')[0] #filter[0] is intended
                set_id = added_ligand_set.id
                message="Ligand Set was successfully created"

        ligand_set = ligand_sets.objects.get(owner=request.user,id=set_id)
        pilog.write("redirecting with message: %s\n" % (message))
        return render_to_response('html/ddligandsetinfo.html', {'ligand_set':ligand_set, 'message':message, 'messagecolor':'Green'})

    #no action, just display blank or populated form
    else: 
        if (ligandsetid!='0'):
            pilog.write("ligandsetid:" + ligandsetid + "\n")
            ligand_set = ligand_sets.objects.get(owner=request.user,id=ligandsetid)
            pilog.write("noaction ligandset:" + str(ligand_set.ligand_set_name) + "\n")
            return render_to_response('html/ddligandsetinfo.html', {'ligand_set':ligand_set, 'message':'', 'messagecolor':'Red'})
        else:
            return render_to_response('html/ddnewligandsetinfo.html', {'description':'', 'message':'', 'messagecolor':'Red'})

def viewLigandGLmol(request, ligand_file_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')
    
    username = request.user.username
    u = User.objects.get(username=username)
    public_user=User.objects.get(username='public')

    ligandfile = files.objects.get(owner__in=[u,public_user],id=ligand_file_id.replace("/",""))
    filename=ligandfile.file_location + "/" + ligandfile.file_name

#    ligandfile=filename.replace("/home/schedd/","/charmming/pdbuploads/")
    #Add logic here for getting the PDB file somehow. Let's do an OBConv.

    obconv = openbabel.OBConversion()
    obconv.SetInAndOutFormats("mol2", "sdf")
    mol = openbabel.OBMol()
    obconv.ReadFile(mol, filename.encode("utf-8"))
    ligand_data = obconv.WriteString(mol)

    try:
        username = request.user.username #???
        return render_to_response('html/ddviewligandglmol.html', {'ligand_data': ligand_data })
    except:
        return HttpResponse("No ligand file found")

#This has been made obsolete by JSmol modal window.
#def viewLigandJmol(request,ligand_file_id):
#
#    if not request.user.is_authenticated():
#        return render_to_response('html/loggedout.html')
#
#    lesson_ok, dd_ok = checkPermissions(request)
#    if not dd_ok:
#        return render_to_response('html/unauthorized.html')
#
#
#
#    #if not request.user.is_authenticated():
#    #    return render_to_response('html/loggedout.html')
#    
#    username = request.user.username
#    u = User.objects.get(username=username)
#    public_user=User.objects.get(username='public')
#    #ligand_filename="vcsdrawn30"
#    #job=jobs.objects.get(owner=request.user,job_scheduler_id=job_id.replace("/",""))
#    #job_folder = charmming_config.user_dd_jobs_home + '/' + username + '/' + 'dd_job_' + str(job.job_owner_index)
#    ligandfile = files.objects.get(owner__in=[u,public_user],id=ligand_file_id.replace("/",""))
#    filename=ligandfile.file_location + "/" + ligandfile.file_name
#
#
#    viewligand=open("/tmp/viewligand.log",'w')
#    #viewpose.write("jobid:%s\n" % job_id)
#    #Getting the pdb
#    
#    #conformation_file=files.objects.filter(owner=request.user).extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_infrastructure_jobs' and object_id=%s)" % (request.user.id,job.id)]).extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_target_protein_conformations')" % (request.user.id)])[0]
#    
#    ligandfile=filename.replace("/home/schedd/","/charmming/pdbuploads/")
#    #proteinfile=conformation_file.file_location+conformation_file.file_name
#    #viewpose.write("babel -imol2 %s -opdb %s" % (proteinfile,proteinfile.replace(".mol2",".pdb")))
#    #os.system("babel -imol2 %s -opdb %s" % (proteinfile,proteinfile.replace(".mol2",".pdb")))
#    #proteinfile=conformation_file.file_location.replace("/var/tmp/","/charmming/")+conformation_file.file_name
#    #proteinfile=proteinfile.replace(".mol2",".pdb").replace("/var/tmp/","/charmming/")
#    viewligand.write("ligandfile: %s\n" % (ligandfile))
#    
#    try:
#        username = request.user.username
#        return render_to_response('html/ddviewligandjmol.html', {'ligandfile': ligandfile })
#    except:
#        return HttpResponse("No ligand file found")
