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
from django.db import transaction
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
#from pdbinfo.models import PDBFile, PDBFileForm, ParseException, energyParams
#from minimization.models import minimizeParams
#from minimization.views import append_tpl
#from dynamics.models import mdParams, ldParams, sgldParams
#from solvation.models import solvationParams
#from account.models import *
#from dynamics.views import combinePDBsForMovie
#from normalmodes.views import combineNmaPDBsForMovie
#from normalmodes.aux import getNormalModeMovieNum
#from normalmodes.models import nmodeParams
#from apbs.models import redoxParams
#from pdbinfo.qmmm import makeQChem_tpl, handleLinkAtoms, writeQMheader

from structure.models import Structure, WorkingStructure, WorkingFile, Task, Segment
from django.contrib.auth.models import User
from django.core import validators 
from django.core.mail import mail_admins
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from account.views import isUserTrustworthy
from account.views import checkPermissions
#from pdbinfo.editscripts import generateHTMLScriptEdit
#from pdbinfo.aux import checkNterPatch
#dd stuff
from dd_infrastructure.models import projects,files,files_objects, jobs, job_types, sources, dsfdockingtask
from dd_target.models import proteins, protein_conformations, projects_protein_conformations
from dd_substrate.models import ligands
from dd_substrate.models import poses
from dd_target import common as target_common
from dd_substrate import common as substrate_common
#import dd_analysis
from dd_analysis.models import attributes, object_attributes
import common
###
import output
import lesson1
import lesson2
import lesson3
import lesson4
import lessonaux
# import all created lessons by importing lesson_config
# also there is a dictionary called 'file_type' in lesson_config.py specififying the file type of files uploaded by the lessons  
from lesson_config import *
import MySQLdb
import MySQLdb.cursors
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
import sys, traceback, shutil
import commands
import charmming_config
import settings
import chemspipy


#Lets user view file in browser, puts it in iframe
#def viewProjectsContainer(request):
#    if not request.user.is_authenticated():
#        return render_to_response('html/loggedout.html')
#
#    return render_to_response('html/ddviewprojectscontainer.html')

#Lets user see the list of projects they have created
def viewProjects(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')

    user_projects = projects.objects.filter(owner=request.user)

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddviewprojects.html', {'user_projects': user_projects,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

def viewUpdateProjectInfo(request,projectid,action):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    pilog=open("/tmp/pilog.log",'w')
    pilog.write("pi\n")
    username = request.user.username
    owner = User.objects.get(username=username)
    lesson_ok, dd_ok = checkPermissions(request)
    #update or add is taking place
    if ((action=='update') or (action=='addnew')):

        if (len(request.POST['name'].strip())==0): 
            message="Operation Failed!<br>Project Name Cannot be Blank"
            if (action=='update'):
                project = projects.objects.get(owner=request.user,id=projectid)
                return render_to_response('html/ddprojectinfo.html', {'project':project, 'message':message, 'messagecolor':'Red','lesson_ok': lesson_ok, 'dd_ok': dd_ok})
            else:
                return render_to_response('html/ddnewprojectinfo.html', {'description':request.POST['description'], 'name':request.POST['name'], 'message':message, 'messagecolor':'Red','lesson_ok': lesson_ok, 'dd_ok': dd_ok})
        else:
            if (action=='update'):
                project = projects.objects.get(owner=request.user,id=projectid)
                project.project_name=request.POST['name'].strip()
                project.description=request.POST['description'].strip()
                project.save()
                project_id=projectid
                message="Project Info was successfully updated"
                
            elif (action=='addnew'):

                if common.ProjectExists(request.POST['name'],request.user)==True:
                    return render_to_response('html/ddnewprojectinfo.html', {'description':request.POST['description'], 'name':request.POST['name'], 'message':'Operation Failed!<br>Duplicate Project Found.', 'messagecolor':'Red','lesson_ok': lesson_ok, 'dd_ok': dd_ok})

                project=projects()
                project.owner=owner
                project.date_created=datetime.datetime.now().isoformat(' ')
                project.project_name=request.POST['name'].strip()
                project.description=request.POST['description'].strip()
                project.save()
            
                added_project = projects.objects.filter(owner=request.user).order_by('-id')[0] #filter[0] is intended
                project_id = added_project.id
                message="Project was successfully created"

        project = projects.objects.get(owner=request.user,id=project_id)
        return render_to_response('html/ddprojectinfo.html', {'project':project, 'message':message, 'messagecolor':'Green'})
    #no action, just display blank or populated form
    else: 
        if (projectid!='0'):
            pilog.write("projectid:" + projectid + "\n")
            project = projects.objects.get(owner=request.user,id=projectid)
            pilog.write("project:" + str(project.project_name) + "\n")
            return render_to_response('html/ddprojectinfo.html', {'project':project, 'message':'', 'messagecolor':'Red','lesson_ok': lesson_ok, 'dd_ok': dd_ok})
        else:
            return render_to_response('html/ddnewprojectinfo.html', {'description':'', 'message':'', 'messagecolor':'Red','lesson_ok': lesson_ok, 'dd_ok': dd_ok})

def projectsaddnew(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')
    user_projects = projects.objects.filter(owner=request.user)
    
    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddprojectsaddnew.html', {'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

#This is for changing the currently selected project
#on the manage projects page
def switchProjects(request,switch_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    try:
        oldproject =  projects.objects.filter(owner=request.user,selected='y')[0]
        oldproject.selected = ''
        oldproject.save()
    except:
        pass
    newproject = projects.objects.filter(owner=request.user,id=switch_id)[0] 
    newproject.selected = 'y'
    newproject.save()

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddswitchproject.html',{'oldproject':oldproject,\
                              'newproject':newproject,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

def projectDetails(request,project_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')
    user_projects = projects.objects.filter(owner=request.user)
    
    try:
       dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)
    except:
       pass
    pdlog=open("/tmp/pdlog.log",'w')
    pdlog.write("prid:" + project_id + "\n")
    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    
    user_proteins=proteins.objects.filter(owner=request.user).order_by('pdb_code')
    """
    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    
    project_conformations_sql = "select fo.file_id as file_id, " \
      "pc.id as conformation_id, pc.conformation_name as conformation_name, " \
      "p.protein_name as protein_name, p.pdb_code as protein_pdb_code " \
      "from dd_infrastructure_files_objects fo " \
      "left join dd_target_protein_conformations pc on fo.object_id=pc.id " \
      "left join dd_target_proteins p on p.id=pc.protein_id " \
      "where fo.object_table_name = 'dd_target_protein_conformations' " \
      "and pc.id in " \
        "(select protein_conformation_id from " \
        "dd_target_projects_protein_conformations where owner_id=%s and " \
        "project_id=%s) " \
      "and pc.owner_id=%s and fo.owner_id=%s" \
      % (request.user.id,project_id,request.user.id,request.user.id)
    pdlog.write("project id:" + project_id + "\n")
    pdlog.write(project_conformations_sql + "\n")
    cursor.execute(project_conformations_sql)
    rows = cursor.fetchall()
    project_conformations=[]
    for row in rows:
        conformation=common.objProteinConformationFile()
        conformation.file_id=row['file_id']     
        conformation.conformation_name=row['conformation_name']
        conformation.protein_pdb_code=row['protein_pdb_code']
        project_conformations.append(conformation)  
   
    cursor.close()
    dba.close()
    """

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddprojectdetails.html', \
                              {'targets':user_proteins, \
                              'project_id':project_id,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

def projectAvailableConformations(request,protein_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')
    paclog=open("/tmp/paclog.log",'w')
    if (protein_id!='0') and (protein_id!=''):
        #try:
        dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)
        cursor = dba.cursor(MySQLdb.cursors.DictCursor)
        """
        conformations_sql = "select fo.file_id as file_id, " \
          "pc.id as conformation_id, pc.conformation_name as conformation_name, " \
          "p.protein_name as protein_name, p.pdb_code as protein_pdb_code " \
          "from dd_infrastructure_files_objects fo " \
          "left join dd_target_protein_conformations pc on fo.object_id=pc.id " \
          "left join dd_target_proteins p on p.id=pc.protein_id " \
          "where fo.object_table_name = 'dd_target_protein_conformations' " \
          "and pc.protein_id = %s " \
          "and pc.owner_id=%s and fo.owner_id=%s" \
          % (protein_id,request.user.id,request.user.id)
        """
        conformations_sql = "select pc.id as conformation_id, " \
                            "pc.conformation_name as conformation_name, " \
                            "p.pdb_code as protein_pdb_code " \
                            "from dd_target_protein_conformations pc " \
                            "left join dd_target_proteins p " \
                            "on pc.protein_id=p.id "\
                            "where p.id = %s and p.owner_id=%s and pc.owner_id=%s" \
                            % (protein_id,request.user.id,request.user.id)
        paclog.write("" + conformations_sql + "\n")
        cursor.execute(conformations_sql)
        rows = cursor.fetchall()
        conformations=[]
        for row in rows:
            conformation=target_common.objProteinConformation()
            conformation.conformation_id=row['conformation_id']
            #conformation.file_id=row['file_id']     
            conformation.conformation_name=row['conformation_name']
            conformation.protein_pdb_code=row['protein_pdb_code']
            conformations.append(conformation)   

        cursor.close()
        dba.close()
        #except:
            #pass
    
    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddavailableconformations.html', \
                              {'available_conformations':conformations,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

@transaction.commit_manually
def projectConformations(request,project_id,addedids,removedids):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    pclog=open("/tmp/pclog.log",'w')
    pclog.write("addedids in:" + addedids +'\n')
    if (project_id!='0') and (project_id!=''):
        if (addedids!="") and (addedids!='0'):
            addedids_list=[int(n) for n in addedids.split(',')]
            pclog.write("addedids:"+str(addedids_list)+'\n')
            added_conformations=protein_conformations.objects.filter(id__in=addedids_list)
            pclog.write("addedc:" + str(added_conformations)+'\n')

            for conformation in added_conformations:
                project_conformation=projects_protein_conformations\
                                    (owner_id=request.user.id,project_id=project_id,
                                    protein_conformation_id=conformation.id)
                project_conformation.save()
                pclog.write("pc saved:" + str(project_conformation) +'\n')
            transaction.commit()
    
        if (removedids!="") and (removedids!='0'):                      
            removedids_list=[int(n) for n in removedids.split(',')]
            projects_protein_conformations.objects.filter(id__in=removedids_list).delete()
            transaction.commit()
            
        
        #try:
        dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)

        cursor = dba.cursor(MySQLdb.cursors.DictCursor)

        
        conformations_sql="select ppc.id as project_conformation_id, " \
                          "pc.id as conformation_id, " \
                          "pc.conformation_name as conformation_name, " \
                          "p.pdb_code as protein_pdb_code " \
                          "from dd_target_projects_protein_conformations ppc " \
                          "left join dd_target_protein_conformations pc on " \
                          "ppc.protein_conformation_id=pc.id " \
                          "left join dd_target_proteins p on " \
                          "pc.protein_id=p.id " \
                          "where ppc.project_id = %s and ppc.owner_id = %s" \
                          % (project_id, request.user.id) 
                        
        pclog.write(conformations_sql + "\n")
        cursor.execute(conformations_sql)
        rows = cursor.fetchall()
        conformations=[]
        for row in rows:
            conformation=common.objProjectProteinConformation()
            conformation.conformation_id=row['conformation_id']
            conformation.project_conformation_id=row['project_conformation_id']
            conformation.conformation_name=row['conformation_name']
            conformation.protein_pdb_code=row['protein_pdb_code']
            conformations.append(conformation)   

        transaction.rollback()
        cursor.close()
        dba.close()
        #except:
            #pass

    
    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddprojectconformations.html', \
                              {'conformations':conformations,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

#def projectD(request,projectid):
#    if not request.user.is_authenticated():
#        return render_to_response('html/loggedout.html')
    

#    pdlog=open("/tmp/plog.log",'w')
#    pdlog.write("prid:" + projectid + "\n")

#    lesson_ok, dd_ok = checkPermissions(request)
#    return render_to_response('html/ddprojectdetails.html'{'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
    
    
#def projectDet(request):
#    return render_to_response('html/ddprojectdet.html')
    

#def DSFContainer(request):
#    if not request.user.is_authenticated():
#        return render_to_response('html/loggedout.html')
#    return render_to_response('html/dddsfcontainer.html')

#displays the daim/seed/ffld job initiation form
def DSFFormDisplay_old(request):
    dsflog=open("/tmp/dsf.log",'w')
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    dsflog.write("authenticated\n") 
    message=""
    username=request.user.username
    u = User.objects.get(username=username)   

    lesson_ok, dd_ok = checkPermissions(request)
    try:
        conformation_info_list=[]
        project = projects.objects.filter(owner=request.user,selected='y')[0]
        project_conformations = projects_protein_conformations.objects.filter(owner=request.user,project=project).select_related()
        for project_conformation in project_conformations:
            conformation_info=target_common.TargetConformationInfo()
            conformation_info.id=project_conformation.protein_conformation.id
            conformation_info.name=project_conformation.protein_conformation.conformation_name
            conformation_info.description=project_conformation.protein_conformation.description
            conformation_info.target_pdb_code=project_conformation.protein_conformation.protein.pdb_code
            
            conformation_info.file_objects_list = files_objects.objects.filter\
                                        (owner=request.user,object_table_name="dd_target_protein_conformations", \
                                        object_id=project_conformation.protein_conformation.id).select_related()
                                        
            conformation_info_list.append(conformation_info)
    except:
        pass
        
    #dsflog.write("got conformations\n")
    try:
        dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)

    except:
        pass
       
    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    public_user=User.objects.get(username='public')
    setssql="select '00' as setid, 'All User Uploaded' as setname, 'All Ligands available to the user' as setdescription, " \
            "(select count(distinct id) from dd_substrate_ligands where owner_id in (%s)) as ligandcount " \
            " union " \
            "select '0' as setid, 'All Available' as setname, 'All Ligands available to the user' as setdescription, " \
            "(select count(distinct id) from dd_substrate_ligands where owner_id in (%s,%s)) as ligandcount " \
	        " union " \
	        "select ls.id as setid, ls.ligand_set_name as setname, ls.description as setdescription, " \
            "(select count(distinct ligands_id) from dd_substrate_ligand_sets_ligands where " \
            "ligands_set_id=ls.id and owner_id in (%s,%s)) as ligandcount " \
	        "from dd_substrate_ligand_sets ls where ls.owner_id in (%s,%s) and ls.ligand_set_name<>'All Available'" \
	        % (request.user.id,request.user.id,public_user.id,request.user.id,public_user.id,request.user.id,public_user.id)
    #setssql="select '0' as setid, gs.name as setname from vcs_gridsets gs"
    
    cursor.execute(setssql)
    rows = cursor.fetchall()
    sets=[]

    for row in rows:
        set=substrate_common.LigandSetObj()
        set.id=row['setid']
        set.name=row['setname']
        set.ligand_count=row['ligandcount']
        sets.append(set)
        
    if not (request.POST):

        try:
            struct = Structure.objects.get(owner=request.user,selected='y')
        except:
            return HttpResponse("Please submit a structure first.")
        
        try:
            ws = WorkingStructure.objects.get(structure=struct,selected='y')
        except:
            return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

        #try:
        abadfiles=[]
        os.chdir(struct.location)
        start=""
        end=""
        nativefilename=""
        os.system("rm *-native-*.pdb")
        for file in glob.glob(struct.location + '/*-bad*.pdb'):
            dsflog.write("file is %s\n" % file)
            if not "badres" in file and not "segment" in file:
                
                os.system("awk '{if($4!=$res && $1==\"ATOM\") {res=$4; print >> \"%s-native-\"res\".pdb\";}}' %s" % (os.path.basename(file)[0:1],os.path.basename(file)))
        for native_lig_file in glob.glob(struct.location + '/*-native-*.pdb'):
            nativefilename=os.path.basename(native_lig_file)
            abadfile=common.objFile()
            abadfile.fullpath=native_lig_file
            abadfile.name=nativefilename
            abadfile.tag=nativefilename[0+len("a-native-"):nativefilename.index(".pdb",0+len("a-native-"))]
            #abadfile.tag=nativefilename[nativefilename.index("a-native-")+len("a-native-"):nativefilename.index(".pdb",nativefilename.index("a-native-")+len("a-native-"))]
            dsflog.write("adding native ligand: %s\n" % abadfile.tag)
            abadfiles.append(abadfile)
            
        dsflog.write("native ligands: %s\n" % str(abadfiles[0].name))
        #except:
        #    return HttpResponse("No native ligands found in the structure. Cannot identify binding site")

        tasks = Task.objects.filter(workstruct=ws,status='C',active='y').exclude(action='energy')
        return render_to_response('html/dddsfform.html', {'abadfiles':abadfiles,'ws_identifier':ws.identifier,'tasks':tasks, 'conformations':conformation_info_list,'existingsets':sets,'message':message,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
    
    
    #####this block contains the actions taking place when the form gets submitted
    ##### i.e. the job gets launched
    tar_files=[]
    job_owner_id=common.getNewJobOwnerIndex(u)
    dsflog.write("job_owner_id: %s\n" % (job_owner_id) )
    job_basename='dd_job_' + str(job_owner_id)
    job_folder = charmming_config.user_dd_jobs_home + '/' + username + '/' + job_basename
    os.system("mkdir %s" % (job_folder))
    os.system("chmod -R g+w %s" % (job_folder))
    os.chdir(job_folder)
    #os.system("pwd > pwd") 
    os.system("touch %s/ligands" % (job_folder))
    os.system("mkdir %s/Inputs" % (job_folder))
    os.system("chmod -R g+w %s/Inputs" % (job_folder))
    os.system("chown schedd %s/Inputs" % (job_folder))
    os.system("mkdir %s/Results" % (job_folder))
    os.system("chmod -R g+w %s/Results" % (job_folder))
    os.system("chown schedd %s/Results" % (job_folder))
    os.system("mkdir %s/ligands_to_dock" % (job_folder))
    os.system("chmod -R g+w %s/ligands_to_dock" % (job_folder))
    os.system("chown schedd %s/ligands_to_dock" % (job_folder))
    os.system("cp %s/charmm_prot_convert.sh %s\n" % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/cgenff_convert.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/targetprep.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/typeassign.py %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/ffld_paramconvert.pl %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/cgenff.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/cgenff_lig.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/add_substructure.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/add_lig_substructure.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/Run_FFLD_FLEA_yp.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/num_substructures.pl %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/pdb2mol.in %s\n' % (charmming_config.dd_scripts_home,job_folder))

    #####charmm stuff
    os.system("mkdir %s/charmm_mini" % (job_folder))
    os.system("mkdir %s/charmm_mini/lig" % (job_folder))
    os.system("mkdir %s/charmm_mini/combined" % (job_folder))
    os.system("mkdir %s/charmm_mini/target" % (job_folder))
    os.system("mkdir %s/charmm_mini/poses" % (job_folder))
    os.system("mkdir %s/charmm_mini/minimized" % (job_folder))
    os.system("cp %s* %s/charmm_mini/\n" % (charmming_config.charmm_files,job_folder))
    os.system("cp %s/*.inp %s/charmm_mini/\n" % (charmming_config.dd_scripts_home,job_folder))
    os.system("cp %s/mol2crd %s/charmm_mini/\n" % (charmming_config.dd_scripts_home,job_folder))
    os.system("cp %s/run_mini.sh %s/charmm_mini/\n" % (charmming_config.dd_scripts_home,job_folder))
    os.system("cp %s/pdb2mol.sh  %s/charmm_mini/\n" % (charmming_config.dd_scripts_home,job_folder))
    ####end charmm stuff

    ####ffld scoring stuff
    os.system("mkdir %s/ffld_eval" % (job_folder))
    os.system("mkdir %s/ffld_eval/poses" % (job_folder))
    os.system("mkdir %s/ffld_eval/ligands" % (job_folder))
    os.system("cp %s/run_ffld_eval.sh  %s/ffld_eval/\n" % (charmming_config.dd_scripts_home,job_folder))
    ####

    ####seed scoring stuff
    os.system("mkdir %s/seed_eval" % (job_folder))
    os.system("mkdir %s/seed_eval/poses" % (job_folder))
    os.system("cp %s/run_seed_eval.sh  %s/seed_eval/\n" % (charmming_config.dd_scripts_home,job_folder))
    os.system("cp %s/make_inp.sh  %s/seed_eval/\n" % (charmming_config.dd_scripts_home,job_folder))
    ####

    os.system("chmod -R g+w %s" % (job_folder))
    os.system("chmod -R g+x %s" % (job_folder))
    #copy native ligand into the job folder
    #native_ligand_destination = open(job_folder + "/native_ligand.pdb", 'w')
    try:
        struct = Structure.objects.get(owner=request.user,selected='y')
    except:
        return HttpResponse('No structure')
    
    filename = struct.location + '/' + request.POST['native_ligand']
    shutil.copy2(filename,job_folder + "/native_ligand.pdb")
    
    #for fchunk in request.FILES['native_ligand_file'].chunks():
    #for fchunk in request.POST['native_ligand'].chunks():
    #native_ligand_destination.write(fchunk)
    
    daim_ligands_list_file=open(job_folder + "/ligands",'w')

    user_ligands=ligands.objects.filter(Q(owner=request.user) | Q(owner=public_user))
    user_ligand_ids=[]
    wnp_ligprep_string=""
    for user_ligand in user_ligands:
        user_ligand_ids.append(user_ligand.id)
        

    dsflog.write("userligands: %s\n" % user_ligand_ids)
    ligandcount=0
    ligandfile_list=[]
    try:
        user_ligand_file_objects=files_objects.objects.filter(Q(owner=request.user) | Q(owner=public_user), \
                                                     object_table_name='dd_substrate_ligands', \
                                                     object_id__in=user_ligand_ids)
        dsflog.write("lignad file objects: %s\n" % user_ligand_file_objects)                                             
        job_ligand_file_ids=[]                                      
        for user_ligand_file_object in user_ligand_file_objects:
            try:
                dsflog.write("dbligandfile: %s\n" % str(user_ligand_file_object.file_id))
                dsflog.write("request ligand %s: %s\n" % (str(user_ligand_file_object.file_id),request.POST["ligand_" + str(user_ligand_file_object.file_id)]))
                if (request.POST["ligand_" + str(user_ligand_file_object.file_id)]=='on'):
                    ligand_file=files.objects.get(id=user_ligand_file_object.file_id)
                    job_ligand_file_ids.append(ligand_file.id)
                    dsflog.write("cp %s%s %s/%s\n" % (ligand_file.file_location, ligand_file.file_name,\
                                            job_folder, ligand_file.file_name))
                    os.system("cp %s%s %s/%s_%s" % (ligand_file.file_location, ligand_file.file_name,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name))
                    strfile=("%s%s" % (ligand_file.file_location, ligand_file.file_name.replace(".mol2",".str")))
                    psffile=("%s%s" % (ligand_file.file_location, ligand_file.file_name.replace(".mol2",".psf")))
                    rtffile=("%s%s" % (ligand_file.file_location, ligand_file.file_name.replace(".mol2",".rtf")))
                    dsflog.write("checking for str: %s\n" % (strfile))
                    if (os.path.isfile(strfile)):
                        os.system("cp %s %s/%s_%s" % (strfile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".str")))
                        dsflog.write("cp %s %s/%s_%s" % (strfile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".str")))
                        os.system("cp %s %s/%s_%s" % (rtffile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".rtf")))
                        dsflog.write("cp %s %s/%s_%s" % (rtffile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".rtf")))
                        os.system("cp %s %s/charmm_mini/lig/%s_%s" % (strfile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".str")))
                        dsflog.write("cp %s %s/charmm_mini/lig/%s_%s" % (strfile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".str")))
                        os.system("cp %s %s/charmm_mini/lig/%s_%s" % (psffile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".psf")))
                        dsflog.write("cp %s %s/charmm_mini/lig/%s_%s" % (psffile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".psf")))
                       

                    ligandcount=ligandcount+1
                    ligandfile_list.append(str(ligand_file.owner_id)+"_"+ligand_file.file_name)
                    daim_ligands_list_file.write(str(ligand_file.owner_id)+"_"+ligand_file.file_name + '\n')
                    # """
                    # wnp_ligprep_string=wnp_ligprep_string + "read mol2 /var/tmp/dd/jobs/root/dd_job_1/%s\n" % (ligand_file.file_name.replace(".mol2",""))
                    # wnp_ligprep_string=wnp_ligprep_string + "modi hydrogen *
                    # wnp_ligprep_string=wnp_ligprep_string + "atom lig current *
                    # wnp_ligprep_string=wnp_ligprep_string + "atom charge auto *
                    # wnp_ligprep_string=wnp_ligprep_string + "atom type auto charmm *
                    # wnp_ligprep_string=wnp_ligprep_string + "done
                    # wnp_ligprep_string=wnp_ligprep_string + "atom name frag -done -auto seq
                    # wnp_ligprep_string=wnp_ligprep_string + "atom q * -done mpeoe
                    # wnp_ligprep_string=wnp_ligprep_string + "done
                    # wnp_ligprep_string=wnp_ligprep_string + "write mol2 /var/tmp/dd/jobs/root/dd_job_1/3_ligand_1_1 ZINC02796227
                    # """
            except:
                pass
    except:
        pass
          
    daim_ligands_list_file.close()                   
    submitscript=open(job_folder + "/submitscript.inp",'w')
                                       
    submitscript.write("#! /bin/bash\n")
    submitscript.write("cd %s\n" % (job_folder))
    submitscript.write("sleep 5\n")
    submitscript.write("while [ ! -s %s/native_ligand.pdb ]\n" % (job_folder))
    submitscript.write("do\n")
    submitscript.write("sleep 5\n")
    submitscript.write("echo 'waiting for native ligand'\n") 
    submitscript.write("done\n")
    submitscript.write("cp %s %s/.daim/\n" % (charmming_config.daim_param,charmming_config.user_home))
    submitscript.write("cp %s %s/.daim/\n" % (charmming_config.daim_prop,charmming_config.user_home))
    submitscript.write("babel -ipdb %s/native_ligand.pdb -omol2 %s/native_ligand.mol2\n" % (job_folder,job_folder))
    current_project=projects.objects.get(owner=request.user,selected='y')
    
    
    
    #project_conformations=projects_protein_conformations.objects.filter(owner=request.user, \
    #                                                          project=current_project)
    conformation_ids=[]
    for conformation in project_conformations:
        conformation_ids.append(conformation.protein_conformation_id)
                                                                  
    dsflog.write("confids: " + str(conformation_ids) + "\n")
    #DAIM
    job_conformation_file_ids=[]                                                 
    #try:
    '''
    conformations_file_objects=files_objects.objects.filter(owner=request.user, \
                                                 object_table_name='dd_target_protein_conformations', \
                                                 object_id__in=conformation_ids)
    for conformation_file_object in conformations_file_objects:
        try:
            dsflog.write("dbconffile: %s\n" % str(conformation_file_object.file_id))
            if (request.POST["conformation_" + str(conformation_file_object.file_id)]=='on'):
                conformation_file=files.objects.get(id=conformation_file_object.file_id)
                job_conformation_file_ids.append(conformation_file.id)
                dsflog.write("cp %s%s %s/%s" % (conformation_file.file_location, conformation_file.file_name,\
                                            job_folder, conformation_file.file_name))
                dsflog.write("cp %s%s %s/%s" % (conformation_file.file_location, conformation_file.file_name.replace(".pdb",".psf"),\
                                            job_folder, conformation_file.file_name.replace(".pdb",".psf")))
                os.system("cp %s%s %s/%s" % (conformation_file.file_location, conformation_file.file_name,\
                                            job_folder, conformation_file.file_name))
                os.system("cp %s%s %s/%s" % (conformation_file.file_location, conformation_file.file_name.replace(".pdb",".psf"),\
                                            job_folder, conformation_file.file_name.replace(".pdb",".psf")))
                #SEED wants the receptor file to be in ./Inputs
                #os.system("cp %s%s %s/Inputs/%s" % (conformation_file.file_location, conformation_file.file_name,\
                #                            job_folder, conformation_file.file_name))
                #####os.system("cp %s%s %s/Inputs/%s\n" % (conformation_file.file_location, conformation_file.file_name,\
                #####                   job_folder, conformation_file.file_name))                            
                ###submitscript.write(charmming_config.daim_exe + " --seed --receptor " + \
                ###                  conformation_file.file_name + " --native " + str(ligand_file.owner_id)+"_"+ligand_file.file_name + \
                ###                  " -c " + str(ligandcount) + " < ligands\n")
                submitscript.write("%s/targetprep.sh %s/%s\n" % (job_folder,job_folder,os.path.splitext(conformation_file.file_name)[0]))
                submitscript.write("cd %s\n" % (job_folder))
                submitscript.write("while [ ! -s %s/%s.mol2 ]\n" % (job_folder,os.path.splitext(conformation_file.file_name)[0]))
                submitscript.write("do\n")
                submitscript.write("sleep 5\n")
                submitscript.write("echo 'waiting for target'\n") 
                submitscript.write("done\n")
                submitscript.write("echo 'found %s.mol2'\n"  % (os.path.splitext(conformation_file.file_name)[0]))
                submitscript.write(charmming_config.daim_exe + " --seed --receptor " + \
                                  os.path.splitext(conformation_file.file_name)[0] + ".mol2 --native " + "native_ligand.mol2" + \
                                  " -c " + str(ligandcount) + " < " + job_folder + "/ligands\n")
                submitscript.write ("sed -i 's/\\sp\\s/ b /g' targ_seed.inp\n")
                submitscript.write("cp %s/%s.mol2 %s/Inputs/\n" % (job_folder,os.path.splitext(conformation_file.file_name)[0],job_folder))
                submitscript.write("%s/add_lig_substructure.sh\n" % (job_folder))  
                               
        except:
            pass
    
    '''
    ##################
    #CHARMMING integrated target


    try:
        struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    try:
        ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
       return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    
    try:
        oldtsk = dsfdockingtask.objects.filter(workstruct=ws,active='y')[0]
        oldtsk.active = 'n'
        oldtsk.save()
    except:
         pass

    dsftask = dsfdockingtask()
    dsftask.setup(ws)
    dsftask.active = 'y'
    dsftask.action = 'dsfdocking'
    dsftask.save()

    if ws.isBuilt != 't':
        isBuilt = False
        pTask = ws.build(dsftask)
        pTaskID = pTask.id
    else:
        isBuilt = True
        pTaskID = int(request.POST['ptask'])


    pTask=Task.objects.get(id=pTaskID)
    dsflog.write("location : %s\n" % (dsftask.workstruct.structure.location))
    dsflog.write("identifier: %s\n" % (dsftask.workstruct.identifier))
    dsflog.write("identifier action: %s\n" % (dsftask.workstruct.identifier + '-' + pTask.action))
    dsflog.write("job folder: %s\n" % (job_folder))
    #struct=Structure.objects.get(id=t.id)
    
    pro_segment=Segment.objects.filter(structure=dsftask.workstruct.structure,type="pro",is_working="y")[0] #filter[0] is on purpose
    pro_name=pro_segment.name + '-' + str(pro_segment.id)
    #pro_name=dsftask.workstruct.identifier + '-' + pTask.action
    dsflog.write("cp %s/%s.pdb %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    dsflog.write("cp %s/%s.psf %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    dsflog.write("cp %s/%s.crd %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    dsflog.write("cp %s/%s.psf %s/charmm_mini/target/" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    dsflog.write("cp %s/%s.crd %s/charmm_mini/target/" % (dsftask.workstruct.structure.location,pro_name,job_folder))

    os.system("cp %s/%s.pdb %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    os.system("cp %s/%s.psf %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    os.system("cp %s/%s.crd %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    os.system("cp %s/%s.psf %s/charmm_mini/target/" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    os.system("cp %s/%s.crd %s/charmm_mini/target" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    #os.system("cp %s/%s_vmd.mol2 %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    #os.system("cp %s/%s.pdb %s" % (dsftask.workstruct.structure.location,pro_segment.name + '-' + str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.psf %s" % (dsftask.workstruct.structure.location,pro_segment.name + '-'+ str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.crd %s" % (dsftask.workstruct.structure.location,pro_segment.name + '-'+ str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.psf %s/charmm_mini/target/" % (dsftask.workstruct.structure.location,pro_segment.name + '-'+ str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.crd %s/charmm_mini/target" % (dsftask.workstruct.structure.location,pro_segment.name + '-'+ str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.pdb %s" % (dsftask.workstruct.structure.location, dsftask.workstruct.identifier + '-' + pTask.action,job_folder))
    #os.system("cp %s/%s.psf %s" % (dsftask.workstruct.structure.location, dsftask.workstruct.identifier + '-' + pTask.action,job_folder))
    
    submitscript.write("%s -dispdev text -e pdb2mol.in -args %s >> pdb2mol.out\n" % (charmming_config.vmd_exe,pro_name))
    submitscript.write("%s/targetprep.sh %s/%s\n" % (job_folder,job_folder,pro_name))
    submitscript.write("cd %s\n" % (job_folder))
    submitscript.write("while [ ! -s %s/%s.mol2 ]\n" % (job_folder,pro_name))
    submitscript.write("do\n")
    submitscript.write("sleep 5\n")
    submitscript.write("echo 'waiting for target'\n")
    submitscript.write("done\n")
    submitscript.write("echo 'found %s.mol2'\n"  % (pro_name))
    submitscript.write(charmming_config.daim_exe + " --seed --receptor " + \
                      pro_name + ".mol2 --native " + "native_ligand.mol2" + \
                      " -c " + str(ligandcount) + " < " + job_folder + "/ligands\n")
    submitscript.write ("sed -i 's/\\sp\\s/ b /g' %s_seed.inp\n" % (pro_name[:4]))
    submitscript.write("cp %s/%s.mol2 %s/Inputs/\n" % (job_folder,pro_name,job_folder))
    submitscript.write("%s/add_lig_substructure.sh\n" % (job_folder))


    #end Charmming integrated target
    #################


    #except:
    #    pass
        
    #os.system("cp %s %s" % (charmming_config.daim_param,job_folder))
    #os.system("cp %s %s" % (charmming_config.daim_prop,job_folder))
    
    os.system("cp %s %s\n" % (charmming_config.daim_param,job_folder))
    os.system("cp %s %s\n" % (charmming_config.daim_prop,job_folder))
    
    
    #os.system("cp %s %s/Inputs/seed.par" % (charmming_config.seed_param, job_folder))
    #os.system("cp %s %s/Inputs/FFLD_param" % (charmming_config.ffld_param, job_folder))
    os.system("cp %s %s/Inputs/seed.par\n" % (charmming_config.seed_param, job_folder))
    #os.system("cp %s %s/FFLD_param_cgenff\n" % (charmming_config.ffld_param, job_folder))
    os.system("cp %s %s\n" % (charmming_config.flea_param, job_folder))
    os.system("cp %s %s/\n" % (charmming_config.charmm_param, job_folder))    
    #os.system("chmod -R g+w %s/Inputs" % (job_folder))
    
    
                                                         
    #submitscript.write('ls Inputs\n')
    #submitscript.write('for f in %s/Inputs/*_frag_*.mol2\n' % (job_folder))
    #submitscript.write('do\n')
    #submitscript.write('    echo "read mol2 ${f%.*}" > wnp_fix_partial_charges\n')
    #submitscript.write('    echo "modi hydrogen *" >> wnp_fix_partial_charges\n')
    #submitscript.write('    echo "atom lig current *" >> wnp_fix_partial_charges\n')
    #submitscript.write('    echo "atom charge auto *" >> wnp_fix_partial_charges\n')
    #submitscript.write('    echo "atom type auto charmm *" >> wnp_fix_partial_charges\n')
    #submitscript.write('    echo "done" >> wnp_fix_partial_charges\n')
    #submitscript.write('    echo "atom name frag -done -auto seq" >> wnp_fix_partial_charges\n')
    #submitscript.write('    echo "atom q * -done mpeoe" >> wnp_fix_partial_charges\n')
    #submitscript.write('    echo "done" >> wnp_fix_partial_charges\n')
    # """
    # submitscript.write('    echo "modify atom q **** -done MPEOE" >> wnp_fix_partial_charges\n')
    # submitscript.write('    echo "done" >> wnp_fix_partial_charges\n')
    # """ 
    #submitscript.write('    bf=${f##*/}\n')
    #submitscript.write('    echo "write mol2 ${f%.*} ${bf%.*} " >> wnp_fix_partial_charges\n')
    #submitscript.write('    %s host none < wnp_fix_partial_charges\n' % (charmming_config.witnotp_exe))
    #submitscript.write('    echo `mv ${f%.*}.mol ${f%.*}.mol2`\n')
    #submitscript.write("    echo `sed -i 's/N.?/NP /' ${f%.*}.mol2`\n")
    #submitscript.write('done\n')
    
   
    submitscript.write('%s/cgenff.sh %s\n' % (job_folder,job_folder))
    submitscript.write('%s/cgenff_lig.sh %s\n' % (job_folder,job_folder))
    submitscript.write("perl ffld_paramconvert.pl > FFLD_param_cgenff\n")    
    submitscript.write('%s %s_seed.inp > %s/seed.out\n' % (charmming_config.seed_exe,pro_name[:4],job_folder))
    submitscript.write('cp %s/outputs/* %s/Inputs/\n' % (job_folder,job_folder))
    submitscript.write('cp `more %s/ligands` %s/ligands_to_dock/\n' % (job_folder,job_folder))
    submitscript.write('%s/Run_FFLD_FLEA_yp.sh\n' % (job_folder))
    submitscript.write('cd %s/charmm_mini/\n' % (job_folder))
    submitscript.write('sh run_mini.sh %s %s\n' % (pro_name,charmming_config.charmm_exe))
    submitscript.write('cd %s/seed_eval/\n' % (job_folder))
    submitscript.write('sh run_seed_eval.sh %s %s\n' % (charmming_config.seed_exe,pro_name))
    submitscript.write('cd %s/ffld_eval/\n' % (job_folder))
    submitscript.write('sh run_ffld_eval.sh %s %s\n' % (charmming_config.ffld_exe,pro_name))
   

    #submitscript.write('if [ "$(ls -A \'%s/clustering/clusters\')" ]; then\n' % (job_folder))
    #submitscript.write('    echo "NORMAL TERMINATION"\n')
    #submitscript.write('fi\n')


    #for ligandfile in ligandfile_list:
        #submitscript.write('cp %s/INFO/%s %s/%s\n' % (job_folder,ligandfile.replace(".mol2","_mis.info"), \
        #                                           job_folder,ligandfile.replace(".mol2",".info")))
        #submitscript.write('%s --outputname ffld_out_%s.pdb %s\n' % (charmming_config.ffld_exe,ligandfile.replace(".mol2",""),ligandfile))
        
        #submitscript.write("count=`tail -1 %s/ffld_out_%s.pdb | awk '{print $5}'`\n" % (job_folder,ligandfile.replace(".mol2","")))
        #submitscript.write("mkdir %s/Results/%s\n" % (job_folder,ligandfile.replace(".mol2","")))
        #submitscript.write('for i in `seq 0 1 $count`; do grep "     *$i      " %s/ffld_out_%s.pdb > %s/Results/%s/ffld_out_$i.pdb; done\n' \
        #                % (job_folder,ligandfile.replace(".mol2",""),job_folder,ligandfile.replace(".mol2","")))
    
    
    submitscript.close()
    
    os.system("chmod g+w -R %s" % (job_folder))
    os.system("chmod g+x %s" % (job_folder + "/submitscript.inp"))

    
    scriptlist=[]
    scriptlist.append(job_folder + "/submitscript.inp")
    exedict={job_folder + "/submitscript.inp": 'sh'}
    
    ######uncomment this block to turn on job launching
    
    ####Job launching code
    #si = schedInterface()
    #newSchedulerJobID = si.submitJob(request.user.id,job_folder,scriptlist,exedict)
    #dsflog.write("schedulerid: %s\n" % (newSchedulerJobID))
    
    NewJob=jobs()
    NewJob.owner=request.user
    NewJob.job_owner_index=job_owner_id
    NewJob.description=request.POST['job_description']
    NewJob.job_start_time=datetime.datetime.now()
    jobtype=job_types.objects.get(job_type_name='DAIM/SEED/FFLD Docking')
    NewJob.job_type=jobtype
    NewJob.save()
    
    ######### Create downloadable job folder
    os.chdir(charmming_config.user_dd_jobs_home + '/' + username)
    os.system("tar -zcvf %s.tar.gz %s" % (job_folder,job_basename))
    dsflog.write("tar -zcvf %s.tar.gz %s" % (job_folder,job_basename))
    #####create associate tar file with the job                                                                                       
    tar_file=files()
    tar_file.owner=u
    tar_file.file_name=job_basename+".tar.gz"
    tar_file.file_location=charmming_config.user_dd_jobs_home + '/' + username
    tar_file.description = "Archived Downloadable Drug Design Job " + str(job_owner_id) 
    tar_file.save()

    tarfile=common.objFile()
    tarfile.fullpath=tar_file.file_location
    tarfile.name=tar_file.file_name
    tarfile.id=tar_file.id
    #tarfile.tag=nativefilename[0+len("a-native-"):nativefilename.index(".pdb",0+len("a-native-"))]
    tar_files.append(tarfile)
    dsflog.write("appending tarfile %s\n" % (tarfile.name))
    
    new_job_tar_file_object=files_objects()
    new_job_tar_file_object.owner=u
    new_job_tar_file_object.file_id=tar_file.id
    new_job_tar_file_object.object_table_name='dd_infrastructure_jobs'
    new_job_tar_file_object.object_id=NewJob.id
    new_job_tar_file_object.save()
    ######### End of downloadable job folder creation    

    ####CHARMMING integrated association of job with task of the target structure
    dsftask.parent=pTask
    dsftask.start(job_folder,scriptlist,exedict)
    dsftask.save()
    dsftask.workstruct.save()
 
    NewJob.job_scheduler_id=dsftask.jobID
    NewJob.save()
    
    new_job = jobs.objects.filter(owner=u).order_by("-id")[0]
    
    
    #####associate conformation files with this job in the database
    #for job_conformation_file_id in job_conformation_file_ids:
    #    new_job_conformation_file_object=files_objects()
    #    new_job_conformation_file_object.owner=u
    #    new_job_conformation_file_object.file_id=job_conformation_file_id
    #    new_job_conformation_file_object.object_table_name='dd_infrastructure_jobs'
    #    new_job_conformation_file_object.object_id=new_job.id
    #    new_job_conformation_file_object.save()
        
    

    #for job_conformation_file_id in job_conformation_file_ids:
    #    new_job_conformation_file_object=files_objects()
    #    new_job_conformation_file_object.owner=u
    #    new_job_conformation_file_object.file_id=job_conformation_file_id
    #    new_job_conformation_file_object.object_table_name='dd_infrastructure_jobs'
    #    new_job_conformation_file_object.object_id=new_job.id
    #    new_job_conformation_file_object.save()
    

    ######

    #####associate ligand files with this job in the database
    for job_ligand_file_id in job_ligand_file_ids:
        new_job_ligand_file_object=files_objects()
        new_job_ligand_file_object.owner=u
        new_job_ligand_file_object.file_id=job_ligand_file_id
        new_job_ligand_file_object.object_table_name='dd_infrastructure_jobs'
        new_job_ligand_file_object.object_id=new_job.id
        new_job_ligand_file_object.save()
        
    #####create associate run.out file with the job
    output_file=files()
    output_file.owner=u
    output_file.file_name="run.out"
    output_file.file_location=job_folder
    output_file.description = "Job execution output for job dd_job_" + str(job_owner_id)
    output_file.save()
    
    new_output_file = files.objects.filter(owner=u).order_by("-id")[0]
    
    new_job_output_file_object=files_objects()
    new_job_output_file_object.owner=u
    new_job_output_file_object.file_id=new_output_file.id
    new_job_output_file_object.object_table_name='dd_infrastructure_jobs'
    new_job_output_file_object.object_id=new_job.id
    new_job_output_file_object.save()
    #####End of job launching code
    
    
    
    message="New DAIM/SEED/FFLD Docking job (Job Code: dd_job_%s)\n " \
            "was successuflly submitted! Use the link below to download the job files\n" % (job_owner_id)
    message_color = 'green'
    
    dsflog.close()

    lesson_ok, dd_ok = checkPermissions(request)    
    return render_to_response('html/dddsfform.html',{'conformations':conformation_info_list, 'message':message, 'message_color':message_color,'lesson_ok': lesson_ok, 'dd_ok': dd_ok, 'tarfiles':tar_files})


def DSFFormDisplay(request):
    dsflog=open("/tmp/dsf.log",'w')
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    dsflog.write("authenticated\n") 
    message=""
    username=request.user.username
    u = User.objects.get(username=username)   

    lesson_ok, dd_ok = checkPermissions(request)
    try:
        conformation_info_list=[]
        project = projects.objects.filter(owner=request.user,selected='y')[0]
        project_conformations = projects_protein_conformations.objects.filter(owner=request.user,project=project).select_related()
        for project_conformation in project_conformations:
            conformation_info=target_common.TargetConformationInfo()
            conformation_info.id=project_conformation.protein_conformation.id
            conformation_info.name=project_conformation.protein_conformation.conformation_name
            conformation_info.description=project_conformation.protein_conformation.description
            conformation_info.target_pdb_code=project_conformation.protein_conformation.protein.pdb_code
            
            conformation_info.file_objects_list = files_objects.objects.filter\
                                        (owner=request.user,object_table_name="dd_target_protein_conformations", \
                                        object_id=project_conformation.protein_conformation.id).select_related()
                                        
            conformation_info_list.append(conformation_info)
    except:
        pass
        
    #dsflog.write("got conformations\n")
    try:
        dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)

    except:
        pass
       
    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    public_user=User.objects.get(username='public')
    setssql="select '00' as setid, 'All User Uploaded' as setname, 'All Ligands available to the user' as setdescription, " \
            "(select count(distinct id) from dd_substrate_ligands where owner_id in (%s)) as ligandcount " \
            " union " \
            "select '0' as setid, 'All Available' as setname, 'All Ligands available to the user' as setdescription, " \
            "(select count(distinct id) from dd_substrate_ligands where owner_id in (%s,%s)) as ligandcount " \
	        " union " \
	        "select ls.id as setid, ls.ligand_set_name as setname, ls.description as setdescription, " \
            "(select count(distinct ligands_id) from dd_substrate_ligand_sets_ligands where " \
            "ligands_set_id=ls.id and owner_id in (%s,%s)) as ligandcount " \
	        "from dd_substrate_ligand_sets ls where ls.owner_id in (%s,%s) and ls.ligand_set_name<>'All Available'" \
	        % (request.user.id,request.user.id,public_user.id,request.user.id,public_user.id,request.user.id,public_user.id)
    #setssql="select '0' as setid, gs.name as setname from vcs_gridsets gs"
    
    cursor.execute(setssql)
    rows = cursor.fetchall()
    sets=[]

    for row in rows:
        set=substrate_common.LigandSetObj()
        set.id=row['setid']
        set.name=row['setname']
        set.ligand_count=row['ligandcount']
        sets.append(set)
        
    if not (request.POST):

        try:
            struct = Structure.objects.get(owner=request.user,selected='y')
        except:
            return HttpResponse("Please submit a structure first.")
        
        try:
            ws = WorkingStructure.objects.get(structure=struct,selected='y')
        except:
            return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

        #try:
        abadfiles=[]
        os.chdir(struct.location)
        start=""
        end=""
        nativefilename=""
        os.system("rm *-native-*.pdb")
        for file in glob.glob(struct.location + '/*-bad*.pdb'):
            dsflog.write("file is %s\n" % file)
            if not "badres" in file and not "segment" in file:
                
                os.system("awk '{if($4!=$res && $1==\"ATOM\") {res=$4; print >> \"%s-native-\"res\".pdb\";}}' %s" % (os.path.basename(file)[0:1],os.path.basename(file)))
        for native_lig_file in glob.glob(struct.location + '/*-native-*.pdb'):
            nativefilename=os.path.basename(native_lig_file)
            abadfile=common.objFile()
            abadfile.fullpath=native_lig_file
            abadfile.name=nativefilename
            abadfile.tag=nativefilename[0+len("a-native-"):nativefilename.index(".pdb",0+len("a-native-"))]
            #abadfile.tag=nativefilename[nativefilename.index("a-native-")+len("a-native-"):nativefilename.index(".pdb",nativefilename.index("a-native-")+len("a-native-"))]
            dsflog.write("adding native ligand: %s\n" % abadfile.tag)
            abadfiles.append(abadfile)
            
        dsflog.write("native ligands: %s\n" % str(abadfiles[0].name))
        #except:
        #    return HttpResponse("No native ligands found in the structure. Cannot identify binding site")

        tasks = Task.objects.filter(workstruct=ws,status='C',active='y').exclude(action='energy')
        return render_to_response('html/dddsfform.html', {'abadfiles':abadfiles,'ws_identifier':ws.identifier,'tasks':tasks, 'conformations':conformation_info_list,'existingsets':sets,'message':message,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
    
    
    #####this block contains the actions taking place when the form gets submitted
    ##### i.e. the job gets launched
    tar_files=[]
    job_owner_id=common.getNewJobOwnerIndex(u)
    dsflog.write("job_owner_id: %s\n" % (job_owner_id) )
    job_basename='dd_job_' + str(job_owner_id)
    job_folder = charmming_config.user_dd_jobs_home + '/' + username + '/' + job_basename
    os.system("mkdir %s" % (job_folder))
    os.system("chmod -R g+w %s" % (job_folder))
    os.chdir(job_folder)
    #os.system("pwd > pwd") 
    os.system("touch %s/ligands" % (job_folder))
    os.system("mkdir %s/Inputs" % (job_folder))
    os.system("chmod -R g+w %s/Inputs" % (job_folder))
    os.system("chown schedd %s/Inputs" % (job_folder))
    os.system("mkdir %s/Results" % (job_folder))
    os.system("chmod -R g+w %s/Results" % (job_folder))
    os.system("chown schedd %s/Results" % (job_folder))
    os.system("mkdir %s/ligands_to_dock" % (job_folder))
    os.system("chmod -R g+w %s/ligands_to_dock" % (job_folder))
    os.system("chown schedd %s/ligands_to_dock" % (job_folder))
    os.system("cp %s/charmm_prot_convert.sh %s\n" % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/cgenff_convert.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/targetprep.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/typeassign.py %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/ffld_paramconvert.pl %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/cgenff.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/cgenff_lig.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/add_substructure.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/add_lig_substructure.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/Run_FFLD_FLEA_yp.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/num_substructures.pl %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/pdb2mol.in %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/run.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    #os.system("cp %s/dodsfc.sh %s" % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/split_seed_input.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/run_seed.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/dock_and_score.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/submit_dock_and_score.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    os.system('cp %s/submit_dock_iteration.sh %s\n' % (charmming_config.dd_scripts_home,job_folder))
    
    #####charmm stuff
    os.system("mkdir %s/charmm_mini" % (job_folder))
    os.system("mkdir %s/charmm_mini/lig" % (job_folder))
    os.system("mkdir %s/charmm_mini/combined" % (job_folder))
    os.system("mkdir %s/charmm_mini/target" % (job_folder))
    os.system("mkdir %s/charmm_mini/poses" % (job_folder))
    os.system("mkdir %s/charmm_mini/minimized" % (job_folder))
    os.system("cp %s* %s/charmm_mini/" % (charmming_config.charmm_files,job_folder))
    os.system("cp %s/*.inp %s/charmm_mini/" % (charmming_config.dd_scripts_home,job_folder))
    os.system("cp %s/mol2crd %s/charmm_mini/" % (charmming_config.dd_scripts_home,job_folder))
    os.system("cp %s/run_mini_ligand.sh %s/charmm_mini/" % (charmming_config.dd_scripts_home,job_folder))
    os.system("cp %s/pdb2mol.sh  %s/charmm_mini/" % (charmming_config.dd_scripts_home,job_folder))
    
    ####end charmm stuff

    ####ffld scoring stuff
    os.system("mkdir %s/ffld_eval" % (job_folder))
    os.system("mkdir %s/ffld_eval/poses" % (job_folder))
    os.system("mkdir %s/ffld_eval/ligands" % (job_folder))
    os.system("cp %s/run_ffld_eval_ligand.sh  %s/ffld_eval/\n" % (charmming_config.dd_scripts_home,job_folder))
    ####

    ####seed scoring stuff
    os.system("mkdir %s/seed_eval" % (job_folder))
    os.system("mkdir %s/seed_eval/poses" % (job_folder))
    os.system("cp %s/run_seed_eval_ligand.sh  %s/seed_eval/\n" % (charmming_config.dd_scripts_home,job_folder))
    os.system("cp %s/make_seed_inp_ligand.sh  %s/seed_eval/\n" % (charmming_config.dd_scripts_home,job_folder))
    ####

    os.system("chmod -R g+w %s" % (job_folder))
    os.system("chmod -R g+x %s" % (job_folder))
    #copy native ligand into the job folder
    #native_ligand_destination = open(job_folder + "/native_ligand.pdb", 'w')
    try:
        struct = Structure.objects.get(owner=request.user,selected='y')
    except:
        return HttpResponse('No structure')
    
    filename = struct.location + '/' + request.POST['native_ligand']
    shutil.copy2(filename,job_folder + "/native_ligand.pdb")
    
    #for fchunk in request.FILES['native_ligand_file'].chunks():
    #for fchunk in request.POST['native_ligand'].chunks():
    #native_ligand_destination.write(fchunk)
    
    daim_ligands_list_file=open(job_folder + "/ligands",'w')

    user_ligands=ligands.objects.filter(Q(owner=request.user) | Q(owner=public_user))
    user_ligand_ids=[]
    wnp_ligprep_string=""
    for user_ligand in user_ligands:
        user_ligand_ids.append(user_ligand.id)
        

    dsflog.write("userligands: %s\n" % user_ligand_ids)
    ligandcount=0
    ligandfile_list=[]
    try:
        user_ligand_file_objects=files_objects.objects.filter(Q(owner=request.user) | Q(owner=public_user), \
                                                     object_table_name='dd_substrate_ligands', \
                                                     object_id__in=user_ligand_ids)
        dsflog.write("ligand file objects: %s\n" % user_ligand_file_objects)                                             
        job_ligand_file_ids=[]                                      
        for user_ligand_file_object in user_ligand_file_objects:
            try:
                dsflog.write("dbligandfile: %s\n" % str(user_ligand_file_object.file_id))
                dsflog.write("request ligand %s: %s\n" % (str(user_ligand_file_object.file_id),request.POST["ligand_" + str(user_ligand_file_object.file_id)]))
                if (request.POST["ligand_" + str(user_ligand_file_object.file_id)]=='on'):
                    ligand_file=files.objects.get(id=user_ligand_file_object.file_id)
                    job_ligand_file_ids.append(ligand_file.id)
                    dsflog.write("cp %s%s %s/%s\n" % (ligand_file.file_location, ligand_file.file_name,\
                                            job_folder, ligand_file.file_name))
                    os.system("cp %s%s %s/%s_%s" % (ligand_file.file_location, ligand_file.file_name,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name))
                    strfile=("%s%s" % (ligand_file.file_location, ligand_file.file_name.replace(".mol2",".str")))
                    psffile=("%s%s" % (ligand_file.file_location, ligand_file.file_name.replace(".mol2",".psf")))
                    rtffile=("%s%s" % (ligand_file.file_location, ligand_file.file_name.replace(".mol2",".rtf")))
                    dsflog.write("checking for str: %s\n" % (strfile))
                    if (os.path.isfile(strfile)):
                        os.system("cp %s %s/%s_%s" % (strfile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".str")))
                        dsflog.write("cp %s %s/%s_%s" % (strfile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".str")))
                        os.system("cp %s %s/%s_%s" % (rtffile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".rtf")))
                        dsflog.write("cp %s %s/%s_%s" % (rtffile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".rtf")))
                        os.system("cp %s %s/charmm_mini/lig/%s_%s" % (strfile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".str")))
                        dsflog.write("cp %s %s/charmm_mini/lig/%s_%s" % (strfile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".str")))
                        os.system("cp %s %s/charmm_mini/lig/%s_%s" % (psffile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".psf")))
                        dsflog.write("cp %s %s/charmm_mini/lig/%s_%s" % (psffile,\
                                            job_folder, ligand_file.owner_id,ligand_file.file_name.replace(".mol2",".psf")))
                       

                    ligandcount=ligandcount+1
                    ligandfile_list.append(str(ligand_file.owner_id)+"_"+ligand_file.file_name)
                    daim_ligands_list_file.write(str(ligand_file.owner_id)+"_"+ligand_file.file_name + '\n')
                    # """
                    # wnp_ligprep_string=wnp_ligprep_string + "read mol2 /var/tmp/dd/jobs/root/dd_job_1/%s\n" % (ligand_file.file_name.replace(".mol2",""))
                    # wnp_ligprep_string=wnp_ligprep_string + "modi hydrogen *
                    # wnp_ligprep_string=wnp_ligprep_string + "atom lig current *
                    # wnp_ligprep_string=wnp_ligprep_string + "atom charge auto *
                    # wnp_ligprep_string=wnp_ligprep_string + "atom type auto charmm *
                    # wnp_ligprep_string=wnp_ligprep_string + "done
                    # wnp_ligprep_string=wnp_ligprep_string + "atom name frag -done -auto seq
                    # wnp_ligprep_string=wnp_ligprep_string + "atom q * -done mpeoe
                    # wnp_ligprep_string=wnp_ligprep_string + "done
                    # wnp_ligprep_string=wnp_ligprep_string + "write mol2 /var/tmp/dd/jobs/root/dd_job_1/3_ligand_1_1 ZINC02796227
                    # """
            except:
                pass
    except:
        pass
          
    daim_ligands_list_file.close()                   
    settings_file=open(job_folder + "/settings.cfg",'w')
    settings_file.write("#! /bin/bash\n")
                                  
    current_project=projects.objects.get(owner=request.user,selected='y')
    #project_conformations=projects_protein_conformations.objects.filter(owner=request.user, \
    #                                                          project=current_project)
    conformation_ids=[]
    for conformation in project_conformations:
        conformation_ids.append(conformation.protein_conformation_id)
                                                                  
    dsflog.write("confids: " + str(conformation_ids) + "\n")
    #DAIM
    job_conformation_file_ids=[]                                                 
    #try:
    '''
    conformations_file_objects=files_objects.objects.filter(owner=request.user, \
                                                 object_table_name='dd_target_protein_conformations', \
                                                 object_id__in=conformation_ids)
    for conformation_file_object in conformations_file_objects:
        try:
            dsflog.write("dbconffile: %s\n" % str(conformation_file_object.file_id))
            if (request.POST["conformation_" + str(conformation_file_object.file_id)]=='on'):
                conformation_file=files.objects.get(id=conformation_file_object.file_id)
                job_conformation_file_ids.append(conformation_file.id)
                dsflog.write("cp %s%s %s/%s" % (conformation_file.file_location, conformation_file.file_name,\
                                            job_folder, conformation_file.file_name))
                dsflog.write("cp %s%s %s/%s" % (conformation_file.file_location, conformation_file.file_name.replace(".pdb",".psf"),\
                                            job_folder, conformation_file.file_name.replace(".pdb",".psf")))
                os.system("cp %s%s %s/%s" % (conformation_file.file_location, conformation_file.file_name,\
                                            job_folder, conformation_file.file_name))
                os.system("cp %s%s %s/%s" % (conformation_file.file_location, conformation_file.file_name.replace(".pdb",".psf"),\
                                            job_folder, conformation_file.file_name.replace(".pdb",".psf")))
                #SEED wants the receptor file to be in ./Inputs
                #os.system("cp %s%s %s/Inputs/%s" % (conformation_file.file_location, conformation_file.file_name,\
                #                            job_folder, conformation_file.file_name))
                #####os.system("cp %s%s %s/Inputs/%s\n" % (conformation_file.file_location, conformation_file.file_name,\
                #####                   job_folder, conformation_file.file_name))                            
                ###submitscript.write(charmming_config.daim_exe + " --seed --receptor " + \
                ###                  conformation_file.file_name + " --native " + str(ligand_file.owner_id)+"_"+ligand_file.file_name + \
                ###                  " -c " + str(ligandcount) + " < ligands\n")
                submitscript.write("%s/targetprep.sh %s/%s\n" % (job_folder,job_folder,os.path.splitext(conformation_file.file_name)[0]))
                submitscript.write("cd %s\n" % (job_folder))
                submitscript.write("while [ ! -s %s/%s.mol2 ]\n" % (job_folder,os.path.splitext(conformation_file.file_name)[0]))
                submitscript.write("do\n")
                submitscript.write("sleep 5\n")
                submitscript.write("echo 'waiting for target'\n") 
                submitscript.write("done\n")
                submitscript.write("echo 'found %s.mol2'\n"  % (os.path.splitext(conformation_file.file_name)[0]))
                submitscript.write(charmming_config.daim_exe + " --seed --receptor " + \
                                  os.path.splitext(conformation_file.file_name)[0] + ".mol2 --native " + "native_ligand.mol2" + \
                                  " -c " + str(ligandcount) + " < " + job_folder + "/ligands\n")
                submitscript.write ("sed -i 's/\\sp\\s/ b /g' targ_seed.inp\n")
                submitscript.write("cp %s/%s.mol2 %s/Inputs/\n" % (job_folder,os.path.splitext(conformation_file.file_name)[0],job_folder))
                submitscript.write("%s/add_lig_substructure.sh\n" % (job_folder))  
                               
        except:
            pass
    
    '''
    ##################
    #CHARMMING integrated target


    try:
        struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    try:
        ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
       return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    
    try:
        oldtsk = dsfdockingtask.objects.filter(workstruct=ws,active='y')[0]
        oldtsk.active = 'n'
        oldtsk.save()
    except:
         pass

    dsftask = dsfdockingtask()
    dsftask.setup(ws)
    dsftask.active = 'y'
    dsftask.action = 'dsfdocking'
    dsftask.save()

    if ws.isBuilt != 't':
        isBuilt = False
        pTask = ws.build(dsftask)
        pTaskID = pTask.id
    else:
        isBuilt = True
        pTaskID = int(request.POST['ptask'])


    pTask=Task.objects.get(id=pTaskID)
    dsflog.write("location : %s\n" % (dsftask.workstruct.structure.location))
    dsflog.write("identifier: %s\n" % (dsftask.workstruct.identifier))
    dsflog.write("identifier action: %s\n" % (dsftask.workstruct.identifier + '-' + pTask.action))
    dsflog.write("job folder: %s\n" % (job_folder))
    #struct=Structure.objects.get(id=t.id)
    
    pro_segment=Segment.objects.filter(structure=dsftask.workstruct.structure,type="pro",is_working="y")[0] #filter[0] is on purpose
    pro_name=pro_segment.name + '-' + str(pro_segment.id)
    #pro_name=dsftask.workstruct.identifier + '-' + pTask.action
    
    settings_file.write('protein="%s"\n' % (pro_name))
    settings_file.write('vmd_exe="%s"\n' % (charmming_config.vmd_exe))
    settings_file.write('daim_exe="%s"\n' % (charmming_config.daim_exe))
    settings_file.write('seed_exe="%s"\n' % (charmming_config.seed_exe))
    settings_file.write('ffld_exe="%s"\n' % (charmming_config.ffld_exe))
    settings_file.write('flea_exe="%s"\n' % (charmming_config.flea_exe))
    settings_file.write('charmm_exe="%s"\n' % (charmming_config.charmm_exe))
    settings_file.write('submit="%s" #serial execution of docking commands\n' % (charmming_config.dd_submit_serial))
    settings_file.write('#submit="%s" #parallel execution of docking commands\n' % (charmming_config.dd_submit_parallel))
    settings_file.write('docking_iterations=%s\n' % (charmming_config.docking_iterations))
    settings_file.write('clustering_energy_cutoff=%s\n' % (charmming_config.clustering_energy_cutoff))
    settings_file.write('num_energy_evaluations=%s\n' % (charmming_config.ffld_energy_evaluations))
    settings_file.write('num_generations=%s\n' % (charmming_config.ffld_generations))
    settings_file.close()

    dsflog.write("cp %s/%s.pdb %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    dsflog.write("cp %s/%s.psf %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    dsflog.write("cp %s/%s.crd %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    dsflog.write("cp %s/%s.psf %s/charmm_mini/target/" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    dsflog.write("cp %s/%s.crd %s/charmm_mini/target/" % (dsftask.workstruct.structure.location,pro_name,job_folder))

    os.system("cp %s/%s.pdb %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    os.system("cp %s/%s.psf %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    os.system("cp %s/%s.crd %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    os.system("cp %s/%s.psf %s/charmm_mini/target/" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    os.system("cp %s/%s.crd %s/charmm_mini/target" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    #os.system("cp %s/%s_vmd.mol2 %s" % (dsftask.workstruct.structure.location,pro_name,job_folder))
    #os.system("cp %s/%s.pdb %s" % (dsftask.workstruct.structure.location,pro_segment.name + '-' + str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.psf %s" % (dsftask.workstruct.structure.location,pro_segment.name + '-'+ str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.crd %s" % (dsftask.workstruct.structure.location,pro_segment.name + '-'+ str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.psf %s/charmm_mini/target/" % (dsftask.workstruct.structure.location,pro_segment.name + '-'+ str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.crd %s/charmm_mini/target" % (dsftask.workstruct.structure.location,pro_segment.name + '-'+ str(pro_segment.id),job_folder))
    #os.system("cp %s/%s.pdb %s" % (dsftask.workstruct.structure.location, dsftask.workstruct.identifier + '-' + pTask.action,job_folder))
    #os.system("cp %s/%s.psf %s" % (dsftask.workstruct.structure.location, dsftask.workstruct.identifier + '-' + pTask.action,job_folder))

    #end Charmming integrated target
    #################


    #except:
    #    pass
        
    #os.system("cp %s %s" % (charmming_config.daim_param,job_folder))
    #os.system("cp %s %s" % (charmming_config.daim_prop,job_folder))
    
    os.system("cp %s %s\n" % (charmming_config.daim_param,job_folder))
    os.system("cp %s %s\n" % (charmming_config.daim_prop,job_folder))
    
    
    #os.system("cp %s %s/Inputs/seed.par" % (charmming_config.seed_param, job_folder))
    os.system("cp %s %s/FFLD_param_cgenff_base" % (charmming_config.ffld_param, job_folder))
    os.system("cp %s %s/Inputs/seed.par\n" % (charmming_config.seed_param, job_folder))
    #os.system("cp %s %s/FFLD_param_cgenff\n" % (charmming_config.ffld_param, job_folder))
    os.system("cp %s %s\n" % (charmming_config.flea_param, job_folder))
    os.system("cp %s %s/\n" % (charmming_config.charmm_param, job_folder))    
    #os.system("chmod -R g+w %s/Inputs" % (job_folder))
        
    #for ligandfile in ligandfile_list:
        #submitscript.write('cp %s/INFO/%s %s/%s\n' % (job_folder,ligandfile.replace(".mol2","_mis.info"), \
        #                                           job_folder,ligandfile.replace(".mol2",".info")))
        #submitscript.write('%s --outputname ffld_out_%s.pdb %s\n' % (charmming_config.ffld_exe,ligandfile.replace(".mol2",""),ligandfile))
        
        #submitscript.write("count=`tail -1 %s/ffld_out_%s.pdb | awk '{print $5}'`\n" % (job_folder,ligandfile.replace(".mol2","")))
        #submitscript.write("mkdir %s/Results/%s\n" % (job_folder,ligandfile.replace(".mol2","")))
        #submitscript.write('for i in `seq 0 1 $count`; do grep "     *$i      " %s/ffld_out_%s.pdb > %s/Results/%s/ffld_out_$i.pdb; done\n' \
        #                % (job_folder,ligandfile.replace(".mol2",""),job_folder,ligandfile.replace(".mol2","")))
        
    os.system("chmod g+w -R %s" % (job_folder))
    #os.system("chmod g+x %s" % (job_folder + "/run.sh"))

    
    scriptlist=[]
    os.system("touch %s/run.inp" % (job_folder))
    scriptlist.append(job_folder + "/run.inp")
    #exedict={job_folder + "/run.inp": 'bash'}
    exedict={job_folder + "/run.inp": job_folder + '/run.sh'}
    #exedict={"": "bash " + job_folder + "/run.inp >& " + job_folder + "/run.out"}
    ######uncomment this block to turn on job launching
    
    ####Job launching code
    #si = schedInterface()
    #newSchedulerJobID = si.submitJob(request.user.id,job_folder,scriptlist,exedict)
    #dsflog.write("schedulerid: %s\n" % (newSchedulerJobID))
    
    NewJob=jobs()
    NewJob.owner=request.user
    NewJob.job_owner_index=job_owner_id
    NewJob.description=request.POST['job_description']
    NewJob.job_start_time=datetime.datetime.now()
    jobtype=job_types.objects.get(job_type_name='DAIM/SEED/FFLD Docking')
    NewJob.job_type=jobtype
    NewJob.save()
    
    ######### Create downloadable job folder
    os.chdir(charmming_config.user_dd_jobs_home + '/' + username)
    os.system("tar -zcvf %s.tar.gz %s" % (job_folder,job_basename))
    dsflog.write("tar -zcvf %s.tar.gz %s" % (job_folder,job_basename))
    #####create associate tar file with the job                                                                                       
    tar_file=files()
    tar_file.owner=u
    tar_file.file_name=job_basename+".tar.gz"
    tar_file.file_location=charmming_config.user_dd_jobs_home + '/' + username
    tar_file.description = "Archived Downloadable Drug Design Job " + str(job_owner_id) 
    tar_file.save()

    tarfile=common.objFile()
    tarfile.fullpath=tar_file.file_location
    tarfile.name=tar_file.file_name
    tarfile.id=tar_file.id
    #tarfile.tag=nativefilename[0+len("a-native-"):nativefilename.index(".pdb",0+len("a-native-"))]
    tar_files.append(tarfile)
    dsflog.write("appending tarfile %s\n" % (tarfile.name))
    
    new_job_tar_file_object=files_objects()
    new_job_tar_file_object.owner=u
    new_job_tar_file_object.file_id=tar_file.id
    new_job_tar_file_object.object_table_name='dd_infrastructure_jobs'
    new_job_tar_file_object.object_id=NewJob.id
    new_job_tar_file_object.save()
    ######### End of downloadable job folder creation    

    ####CHARMMING integrated association of job with task of the target structure
    dsftask.parent=pTask
    dsftask.start(job_folder,scriptlist,exedict)
    dsftask.save()
    dsftask.workstruct.save()
 
    NewJob.job_scheduler_id=dsftask.jobID
    NewJob.save()
    
    new_job = jobs.objects.filter(owner=u).order_by("-id")[0]
    
    
    #####associate conformation files with this job in the database
    #for job_conformation_file_id in job_conformation_file_ids:
    #    new_job_conformation_file_object=files_objects()
    #    new_job_conformation_file_object.owner=u
    #    new_job_conformation_file_object.file_id=job_conformation_file_id
    #    new_job_conformation_file_object.object_table_name='dd_infrastructure_jobs'
    #    new_job_conformation_file_object.object_id=new_job.id
    #    new_job_conformation_file_object.save()
        
    

    #for job_conformation_file_id in job_conformation_file_ids:
    #    new_job_conformation_file_object=files_objects()
    #    new_job_conformation_file_object.owner=u
    #    new_job_conformation_file_object.file_id=job_conformation_file_id
    #    new_job_conformation_file_object.object_table_name='dd_infrastructure_jobs'
    #    new_job_conformation_file_object.object_id=new_job.id
    #    new_job_conformation_file_object.save()
    

    ######

    #####associate ligand files with this job in the database
    for job_ligand_file_id in job_ligand_file_ids:
        new_job_ligand_file_object=files_objects()
        new_job_ligand_file_object.owner=u
        new_job_ligand_file_object.file_id=job_ligand_file_id
        new_job_ligand_file_object.object_table_name='dd_infrastructure_jobs'
        new_job_ligand_file_object.object_id=new_job.id
        new_job_ligand_file_object.save()
        
    #####create associate submitscript.out file with the job
    output_file=files()
    output_file.owner=u
    output_file.file_name="run.out"
    output_file.file_location=job_folder
    output_file.description = "Job execution output for job dd_job_" + str(job_owner_id)
    output_file.save()
    
    new_output_file = files.objects.filter(owner=u).order_by("-id")[0]
    
    new_job_output_file_object=files_objects()
    new_job_output_file_object.owner=u
    new_job_output_file_object.file_id=new_output_file.id
    new_job_output_file_object.object_table_name='dd_infrastructure_jobs'
    new_job_output_file_object.object_id=new_job.id
    new_job_output_file_object.save()
    #####End of job launching code
    
    
    
    message="New DAIM/SEED/FFLD Docking job (Job Code: dd_job_%s)\n " \
            "was successuflly submitted! Use the link below to download the job files\n" % (job_owner_id)
    message_color = 'green'
    
    dsflog.close()

    lesson_ok, dd_ok = checkPermissions(request)    
    return render_to_response('html/dddsfform.html',{'conformations':conformation_info_list, 'message':message, 'message_color':message_color,'lesson_ok': lesson_ok, 'dd_ok': dd_ok, 'tarfiles':tar_files})



def updateDSFLigands(request,selected_set_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    log=open("/tmp/dsfligands.txt",'w')
    log.write("selecte_set_id: %s\n" % (selected_set_id))
    ligand_info_list=[]
    #project = projects.objects.filter(owner=request.user,selected='y')[0]
    user_ligands = ligands.objects.filter(owner=request.user).select_related()
    log.write("user_ligands: %s" % (str(user_ligands)))
    public_user=User.objects.get(username='public')
    public_ligands = ligands.objects.filter(owner=public_user).select_related()

    if (selected_set_id=='0'): #show all ligands
        for ligand in user_ligands:
            ligand_info=substrate_common.LigandInfo()
            ligand_info.id=ligand.id
            ligand_info.name=ligand.ligand_name
            ligand_info.description=ligand.description
            ligand_info.file_objects_list = files_objects.objects.filter\
                                        (owner=request.user,object_table_name="dd_substrate_ligands", \
                                        object_id=ligand.id).select_related()
                                        
            ligand_info_list.append(ligand_info)
            
        for ligand in public_ligands:
            ligand_info=substrate_common.LigandInfo()
            ligand_info.id=ligand.id
            ligand_info.name=ligand.ligand_name
            ligand_info.description=ligand.description
            ligand_info.file_objects_list = files_objects.objects.filter\
                                        (owner=public_user,object_table_name="dd_substrate_ligands", \
                                        object_id=ligand.id).select_related()
                                        
            ligand_info_list.append(ligand_info)

    elif (selected_set_id=='00'): #user ligands
        for ligand in user_ligands:
            ligand_info=substrate_common.LigandInfo()
            ligand_info.id=ligand.id
            ligand_info.name=ligand.ligand_name
            ligand_info.description=ligand.description
            ligand_info.file_objects_list = files_objects.objects.filter\
                                        (owner=request.user,object_table_name="dd_substrate_ligands", \
                                        object_id=ligand.id).select_related()
                                        
            ligand_info_list.append(ligand_info)
            
    else:
        selected_set_ligands=ligands.objects.filter().extra \
        (where=["dd_substrate_ligands.id IN (select ligands_id from dd_substrate_ligand_sets_ligands where owner_id in (%s,%s) " \
        " and ligands_set_id=%s)" % (request.user.id,public_user.id,selected_set_id)]).select_related()
        log.write("id IN (select ligands_id from dd_substrate_ligand_sets_ligands where owner_id in (%s,%s) " \
        " and ligands_set_id=%s)" % (request.user.id,public_user.id,selected_set_id))
        log.write("selected_set_ligands: %s" % (str(user_ligands)))
        for ligand in selected_set_ligands:
            log.write("ligand id: %s" % (ligand.id))
            ligand_info=substrate_common.LigandInfo()
            ligand_info.id=ligand.id
            ligand_info.name=ligand.ligand_name
            ligand_info.description=ligand.description
            ligand_info.file_objects_list = files_objects.objects.filter\
                                        (object_table_name="dd_substrate_ligands", \
                                        object_id=ligand.id).select_related()
            #ligand_info.file_objects_list = files_objects.objects.filter\
            #                            (object_table_name="dd_substrate_ligands", \
            #                            object_id=ligand.id).select_related()                                
            ligand_info_list.append(ligand_info)

    lesson_ok, dd_ok = checkPermissions(request)        
    return render_to_response('html/dddsfligands.html', {'ligands':ligand_info_list,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})



#Lets user view file in browser, puts it in iframe
#def viewJobsContainer(request):
#    if not request.user.is_authenticated():
#        return render_to_response('html/loggedout.html')
#    return render_to_response('html/ddviewjobscontainer.html')

#Lets user see the list of Drug Design Jobs history
def viewJobs(request):
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

    try:
        struct = Structure.objects.get(owner=request.user,selected='y')
    except:
        return HttpResponse("Please submit a structure first.")

    try:
        ws = WorkingStructure.objects.get(structure=struct,selected='y')
    except:
       return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    jobs=[]   

    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    jobssql = "SELECT js.id, dj.job_scheduler_id, dj.job_owner_index, dj.description, js.queued, " \
              "if (state=1,timediff(now(),js.queued),timediff(started,js.queued)) 'Waited', " \
              "if (state=2,timediff(now(),js.started),timediff(ended,started)) 'Runtime', " \
              "if (state=2,timediff(now(),js.queued),timediff(ended,js.queued)) 'Total Time', " \
              "(select case state when 1 then 'Queued' when 2 then 'Running' when 3 then 'Done' when 4 then 'Failed' end) as 'Status', " \
              "(select case state when 1 then 'Yellow' when 2 then 'Blue' when 3 then 'Green' when 4 then 'Red' end) as 'Status Color' " \
              "FROM job_scheduler js, dd_infrastructure_jobs dj WHERE dj.job_scheduler_id=js.id and " \
              "js.id in (select jobID from structure_task where workstruct_id = %s) " \
              "and js.userid=%d order by dj.job_owner_index desc" % (ws.id,request.user.id)
    cursor.execute(jobssql)
    rows = cursor.fetchall()
    cursor.close()
    dba.close()

    for row in rows:
        job=common.objDDJob()    
        job.id=row['id']
        job.scheduler_id=row['job_scheduler_id']
        job.owner_index=row['job_owner_index']
        job.description=row['description']
        job.status="<font color='%s'>%s</font>" % (row['Status Color'],row['Status'])
        job.queued=row['queued']
        job.waited=row['Waited']
        job.runtime=row['Runtime']
        job.totaltime=row['Total Time']
        jobs.append(job)
        

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddviewjobs.html', {'jobs': jobs,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})


#Lets user view file in browser, puts it in iframe
#def viewJobDetailsContainer(request,job_id):
#    if not request.user.is_authenticated():
#        return render_to_response('html/loggedout.html')
#    return render_to_response('html/viewjobdetailscontainer.html', {'job_id':job_id.replace("/","")})


def viewJobDetails(request,job_id):
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
    rows=[]
    jobs=[]   
    message=""
       

    cursor = dba.cursor(MySQLdb.cursors.DictCursor)
    jobssql = "SELECT dj.id dj_id, js.id, dj.job_scheduler_id, dj.job_owner_index, dj.description, js.queued, " \
              "if (state=1,timediff(now(),js.queued),timediff(started,js.queued)) 'Waited', " \
              "if (state=2,timediff(now(),js.started),timediff(ended,started)) 'Runtime', " \
              "(select case state when 1 then 'Queued' when 2 then 'Running' when 3 then 'Done' when 4 then 'Failed' end) as 'Status', " \
              "(select case state when 1 then 'Yellow' when 2 then 'Blue' when 3 then 'Green' when 4 then 'Red' end) as 'Status Color' " \
              "FROM job_scheduler js, dd_infrastructure_jobs dj WHERE dj.job_scheduler_id=js.id and " \
              "js.userid=%d and js.id=%s order by dj.job_owner_index desc" % (request.user.id,job_id.replace("/",""))
    #jobssql = "SELECT vj.id as vjid, js.id, vj.scheduler_id, vj.owner_job_id, vj.job_description, s.name 'Status', s.color 'Status Color', js.queued, if (s.name='Queued',timediff(now(),js.queued),timediff(started,js.queued)) 'Waited', if (s.name='Running',timediff(now(),js.started),timediff(ended,started)) 'Runtime' FROM job_scheduler js, vcs_jobs vj, job_state s WHERE vj.scheduler_id=js.id and js.userid=%d and s.id=js.state and js.id=%s" % (request.user.id,job_id.replace("/",""))
    cursor.execute(jobssql)
    rows = cursor.fetchall()
    cursor.close()

    jobdetaillog=open("/tmp/jobdetail.log",'w')
    jobdetaillog.write("job detail sql: %s\n" % (jobssql))
    for row in rows:
        job=common.objDDJob()
        job.dj_id=row['dj_id']    
        job.id=row['id']
        job.scheduler_id=row['job_scheduler_id']
        job.owner_index=row['job_owner_index']
        job.description=row['description']
        job.status="<font color='%s'>%s</font>" % (row['Status Color'],row['Status'])
        job.queued=row['queued']
        job.waited=row['Waited']
        job.runtime=row['Runtime']
        
        jobdetaillog.write("checking path: %s\n" % (charmming_config.user_dd_jobs_home + '/' + request.user.username + '/' + 'dd_job_' + str(job.owner_index)))
        jobdetaillog.write("jobpath: %s\n" % (os.path.exists(charmming_config.user_dd_jobs_home + '/' + request.user.username + '/' + 'dd_job_' + str(job.owner_index))))
        if not os.path.exists(charmming_config.user_dd_jobs_home + '/' + request.user.username + '/' + 'dd_job_' + str(job.owner_index)):
            message="This is an old job whose files have been deleted from the server. The information shown on this page is incomplete."
            jobdetaillog.write(message)

#    conformation_files=files.objects.filter(owner=request.user).extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_infrastructure_jobs' and object_id=%s)" % (request.user.id,job.dj_id)]).extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_target_protein_conformations')" % (request.user.id)])
    jobdetaillog.write("initiationg conformations list\n")

    try:
        dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)

    except:
        pass

    rows=[]
    jobs=[]

    cursor = dba.cursor(MySQLdb.cursors.DictCursor)
    conformations_sql = "select distinct p.pdb_code, p.protein_name,c.conformation_name, f.file_id from " \
                        " dd_target_proteins p left join dd_target_protein_conformations c on p.id=c.protein_id, dd_infrastructure_files_objects f " \
                        " where c.id = (select object_id from dd_infrastructure_files_objects where file_id=f.file_id and object_table_name='dd_target_protein_conformations') " \
                        " and f.file_id in (select id from dd_infrastructure_files where id IN " \
                        " (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_infrastructure_jobs' and object_id=%s) and id IN " \
                        " (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_target_protein_conformations'))" \
                        % (request.user.id,job.dj_id,request.user.id)

    jobdetaillog.write("conformations_sql: %s\n" % (conformations_sql))
    cursor.execute(conformations_sql)
    job_conformations = cursor.fetchall()
    cursor.close()

    #ligand_files=files.objects.filter(owner=request.user).extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_infrastructure_jobs' and object_id=%s)" % (request.user.id,job.dj_id)]).extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_substrate_ligands')" % (request.user.id)])
    
    cursor = dba.cursor(MySQLdb.cursors.DictCursor)
    ligands_sql = "select distinct l.ligand_name,l.description, f.file_id from " \
                        " dd_substrate_ligands l, dd_infrastructure_files_objects f " \
                        " where l.id = (select object_id from dd_infrastructure_files_objects where file_id=f.file_id and object_table_name='dd_substrate_ligands') " \
                        " and f.file_id in (select id from dd_infrastructure_files where id IN " \
                        " (select file_id from dd_infrastructure_files_objects where owner_id in (%s,%s) and object_table_name='dd_infrastructure_jobs' and object_id=%s) and id IN " \
                        " (select file_id from dd_infrastructure_files_objects where owner_id in (%s,%s) and object_table_name='dd_substrate_ligands')) " \
                        % (request.user.id,public_user.id,job.dj_id,request.user.id,public_user.id)

    jobdetaillog.write("ligands_sql: %s\n" % (ligands_sql))
    cursor.execute(ligands_sql)
    job_ligands = cursor.fetchall()
    cursor.close()

    job_files=[]
    job_files=files.objects.filter(Q(file_name__iendswith='.out') | Q(file_name__iendswith='.tar.gz'), owner=request.user).extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_infrastructure_jobs' and object_id=%s)" % (request.user.id,job.dj_id)])


    lesson_ok, dd_ok = checkPermissions(request)    
    return render_to_response('html/ddviewjobdetails.html', {'job': job, 'job_conformations':job_conformations, 'job_ligands':job_ligands, 'job_files':job_files, 'message':message, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
    
def viewJobInfo(request,job_id):
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

    rows=[]
    jobs=[]   

    cursor = dba.cursor(MySQLdb.cursors.DictCursor)
    jobssql = "SELECT js.id, dj.job_scheduler_id, dj.job_owner_index, dj.description, js.queued, " \
              "if (state=1,timediff(now(),js.queued),timediff(started,js.queued)) 'Waited', " \
              "if (state=2,timediff(now(),js.started),timediff(ended,started)) 'Runtime', " \
              "if (state=2,timediff(now(),js.queued),timediff(ended,js.queued)) 'Total Time', " \
              "(select case state when 1 then 'Queued' when 2 then 'Running' when 3 then 'Done' when 4 then 'Failed' end) as 'Status', " \
              "(select case state when 1 then 'Yellow' when 2 then 'Blue' when 3 then 'Green' when 4 then 'Red' end) as 'Status Color' " \
              "FROM job_scheduler js, dd_infrastructure_jobs dj WHERE dj.job_scheduler_id=js.id and " \
              "js.userid=%d and js.id=%s order by dj.job_owner_index desc" % (request.user.id,job_id.replace("/",""))
    #jobssql = "SELECT vj.id as vjid, js.id, vj.scheduler_id, vj.owner_job_id, vj.job_description, s.name 'Status', s.color 'Status Color', js.queued, if (s.name='Queued',timediff(now(),js.queued),timediff(started,js.queued)) 'Waited', if (s.name='Running',timediff(now(),js.started),timediff(ended,started)) 'Runtime' FROM job_scheduler js, vcs_jobs vj, job_state s WHERE vj.scheduler_id=js.id and js.userid=%d and s.id=js.state and js.id=%s" % (request.user.id,job_id.replace("/",""))
    cursor.execute(jobssql)
    rows = cursor.fetchall()
    cursor.close()

    jobdetaillog=open("/tmp/jobinfo.log",'w')
    jobdetaillog.write("jobdetail:\n")
    for row in rows:
        job=common.objDDJob()    
        job.id=row['id']
        job.scheduler_id=row['job_scheduler_id']
        job.owner_index=row['job_owner_index']
        job.description=row['description']
        job.status="<font color='%s'>%s</font>" % (row['Status Color'],row['Status'])
        job.queued=row['queued']
        job.waited=row['Waited']
        job.runtime=row['Runtime']
        job.totaltime=row['Total Time']
        
    results=0
    

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddviewjobinfo.html', {'job':job, 'results':results,'lesson_ok': lesson_ok, 'dd_ok': dd_ok })

#@transaction.commit_manually    
def viewJobResults(request,job_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    username=request.user.username
    u = User.objects.get(username=username)
    #transaction.commit()
    result_poses=[]
    resultslog=open("/tmp/results.log","w")
    ###check if result files are in the database
    ###if not in db then see if the result files have been created and if so add them to the database
    #
    resultslog.write("job id: %s\n" % job_id)
    job=jobs.objects.get(owner=request.user,job_scheduler_id=job_id.replace("/",""))
    #transaction.commit()
    try:
        new_source=sources.objects.get(source_object_table_name="dd_infrastructure_jobs",source_object_id=job.id)
        #transaction.commit()
    except:
        new_source=sources()
        new_source.source_name="DSF Docking Job"
        new_source.description="Generated by DSF Docking Protocol"
        new_source.source_object_table_name = "dd_infrastructure_jobs"
        new_source.source_object_id = job.id
        new_source.save()
        #transaction.commit()

    path = charmming_config.user_dd_jobs_home + '/' + username + '/' + 'dd_job_' + \
           str(job.job_owner_index) + "/clustering/clusters/" 
    
    charmm_mini_path = charmming_config.user_dd_jobs_home + '/' + username + '/' + 'dd_job_' + \
           str(job.job_owner_index) + "/charmm_mini/"

    seed_eval_path = charmming_config.user_dd_jobs_home + '/' + username + '/' + 'dd_job_' + \
           str(job.job_owner_index) + "/seed_eval/"

    ffld_eval_path = charmming_config.user_dd_jobs_home + '/' + username + '/' + 'dd_job_' + \
           str(job.job_owner_index) + "/ffld_eval/"
    
    #
    ##try:
        #result_files=files.objects.filter(description='Docking By Decomposition Ligand Pose File') \
        #            .extra(where=["id IN (select file_id from dd_infrastructure_files_objects " + \
        #                         "where owner_id=%s and object_table_name='dd_infrastructure_jobs' and object_id=%s)" \
        #                         % (request.user.id,job.id)]).order_by('file_name')[:10]
    ##    try:
    ##        dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
    ##                                  db=settings.DATABASE_NAME, compress=1)
    ##    except:
    ##        pass

    ##    cursor = dba.cursor(MySQLdb.cursors.DictCursor)
    ##    poses_sql="SELECT @count := @count + 1 as rank, p.id, f.id as 'file_id',f.file_name, f.file_location, " \
    ##              " (select value from dd_analysis_object_attributes where object_table_name='dd_substrate_poses' and object_id=p.id) as 'score', " \
    ##              " (select ligand_name from dd_substrate_ligands where id=p.pose_object_id) as ligand_name " \
    ##              " from dd_substrate_poses p, dd_infrastructure_files f, (select @count := 0) c " \
    ##              " where p.source_id = " \
    ##              " (select id from dd_infrastructure_sources where source_object_table_name='dd_infrastructure_jobs' and source_object_id= %s) " \
    ##              " and f.id = (select file_id from dd_infrastructure_files_objects where object_table_name ='dd_substrate_poses' and object_id=p.id) " \
    ##              " order by (select value from dd_analysis_object_attributes where object_table_name='dd_substrate_poses' and object_id=p.id)" % (str(job.id))
        
    ##    resultslog.write("poses_sql: %s\n" % (poses_sql))
    ##    cursor.execute(poses_sql)
    ##    result_poses = cursor.fetchall()
    ##    cursor.close()


    ##    transaction.commit()
    ##except:
    ##    resultslog.write("poses_sql failed\n")
    ##    transaction.rollback()
    ##    pass
    #filename=result_files[0].file_name
    #user_id=filename[:filename.find("_")]
    #ligand_id_range=str(filename.find("_",filename.find("_")+1)+1) + " " + str(filename.find("_",filename.find("_",filename.find("_")+1)+1))
    #resultslog.write("ligand_id_range: %s\n" % (ligand_id_range))
    #ligand_id=filename[filename.find("_",filename.find("_")+1)+1:filename.find("_",filename.find("_",filename.find("_")+1)+1)]
    #resultslog.write("user_id: %s, ligand_id: %s\n" % (user_id,ligand_id))
    ##resultslog.write("len(poses): %s\n" % (str(len(result_poses))))
    ##resultslog.write("pose files: %s\n" % (str(len(glob.glob( os.path.join(path, '*_clus0.mol2') )))))
    #if len(result_poses)<len(glob.glob( os.path.join(path, '*_clus0.mol2') )):
    if 1==1:
        resultslog.write("result files\n")
        resultslog.write("ospathjoin: " + os.path.join(path, '*_clus*.mol2') + "\n")
        for file in glob.glob( os.path.join(path, '*_clus*.mol2') ):
            try:
                
                posefiles=files.objects.filter(owner=request.user,file_name=os.path.basename(file),file_location=os.path.dirname(file))
                resultslog.write("for file %s posefile: %s\n" % (os.path.basename(file),str(posefiles)))
                if len(posefiles)!=0:
                    continue
                    #i=1
            except:
                #resultslog.write("continuing\n")
                #continue
                pass
            #fileobject=files_objects.objects.get
           
            
            posefile=os.path.basename(file)
            posename=posefile.replace(".mol2","")
            ligname=posename[:posename.index("_clus")]
            #ligname=os.path.basename(file).replace("_clus0.mol2","")
              
            resultslog.write("file: " + file + " ligname: " + ligname + "\n")
            resultslog.write("posefile: " + file + " posename: " + posename + "\n")
            new_file=files()
            new_file.owner=u
            new_file.file_name=os.path.basename(file)
            new_file.file_location=os.path.dirname(file)
            new_file.description="Docking By Decomposition Ligand Pose File"
            #new_file.save()
            
            resultslog.write("new_file saved\n")
            #associate result file with the job
            #added_file = files.objects.filter(owner=request.user).order_by('-id')[0]
            new_file_object=files_objects()
            new_file_object.owner=u
            new_file_object.file_id=new_file.id
            new_file_object.object_table_name="dd_infrastructure_jobs"
            new_file_object.object_id=job.id
            #new_file_object.save()
            
            resultslog.write("new_file_object saved\n")
            #transaction.commit()
            #create a pose object
            new_pose_object=poses()
            new_pose_object.owner=u
            new_pose_object.pose_object_table_name="dd_substrate_ligands"
            filename=new_file.file_name
            file_owner_id=int(filename[:filename.find("_")])
            file_ligand_index=int(filename[filename.find("_",filename.find("_")+1)+1:filename.find("_",filename.find("_",filename.find("_")+1)+1)])
            ligand=ligands.objects.get(owner=file_owner_id,ligand_owner_index=file_ligand_index)
            #transaction.commit()
            new_pose_object.pose_object_id=ligand.id
            new_pose_object.source=new_source
            #new_pose_object.save()

            ##transaction.commit()
            resultslog.write("new_pose_object saved\n")
            #create an attribute for the pose energy and link it with the pose object
            attribute=attributes.objects.get(attribute_short_name="DSF Docking Score")
            new_object_attribute=object_attributes()
            new_object_attribute.owner=u
            new_object_attribute.attribute_id=attribute.id
            new_object_attribute.object_table_name="dd_substrate_poses"
            new_object_attribute.object_id=new_pose_object.id
            new_object_attribute.value=common.GetEnergyFromFLEAMol(new_file.file_location + "/" + new_file.file_name)
            #new_object_attribute.save()

            resultslog.write("new_object_attribute saved\n")
            #link result file with the pose object
            new_file_object1=files_objects()
            new_file_object1.owner=u
            new_file_object1.file_id=new_file.id
            new_file_object1.object_table_name="dd_substrate_poses"
            new_file_object1.object_id=new_pose_object.id
            #new_file_object1.save()
            ##transaction.commit()

            
            ######if charmm minimized pose exists
            if 1==1:
            #try:
                min_posefiles=files.objects.filter(owner=request.user,file_name=posename+"_min.mol2",file_location=charmm_mini_path + "minimized/")
                resultslog.write("for charmm mini file %s \n" % str(min_posefiles))
                if len(min_posefiles)!=0:
                    continue
                    #i=1
            #except:
            #     pass
 
            ###save minimized pose file and create all the objects for it
            
            new_min_file=files()
            new_min_file.owner=u
            new_min_file.file_name=posename+"_min.mol2"
            new_min_file.file_location=charmm_mini_path + "minimized/"
            new_min_file.description="Docking By Decomposition CHARMM Minimized Ligand Pose File"
            new_min_file.save()

            resultslog.write("new_min_file with id: %s saved\n" % (new_min_file.id))

            new_min_file_object=files_objects()
            new_min_file_object.owner=u
            new_min_file_object.file_id=new_min_file.id
            new_min_file_object.object_table_name="dd_infrastructure_jobs"
            new_min_file_object.object_id=job.id
            new_min_file_object.save()

            resultslog.write("new_min_file_object saved\n")
            #transaction.commit()

            #create a pose object                                                                                                                                             
            new_min_pose_object=poses()
            new_min_pose_object.owner=u
            new_min_pose_object.pose_object_table_name="dd_substrate_ligands"
            filename=new_min_file.file_name
            file_owner_id=int(filename[:filename.find("_")])
            file_ligand_index=int(filename[filename.find("_",filename.find("_")+1)+1:filename.find("_",filename.find("_",filename.find("_")+1)+1)])
            ligand=ligands.objects.get(owner=file_owner_id,ligand_owner_index=file_ligand_index)
            #transaction.commit()
            new_min_pose_object.pose_object_id=ligand.id
            new_min_pose_object.source=new_source
            new_min_pose_object.save()
            resultslog.write("new_min_pose_object for ligand %s saved\n" % (ligand.ligand_name))

            #link result file with the pose object                                                                                                                                             
            new_min_file_object1=files_objects()
            new_min_file_object1.owner=u
            new_min_file_object1.file_id=new_min_file.id
            new_min_file_object1.object_table_name="dd_substrate_poses"
            new_min_file_object1.object_id=new_min_pose_object.id
            new_min_file_object1.save()
            resultslog.write("new_min_file_object1 for file %s representing min_pose saved\n" % (new_min_file.file_name))

            ##transaction.commit()                                                                                                                                                             
            #create an attribute for the pose energy and link it with the pose object
            #CHARMM minimization
            charmm_vdw,charmm_elec,charmm_tot=common.GetEnergyFromCHARMMMini(charmm_mini_path,posename)
            resultslog.write("charmm_vdw: %s, charmm_elec: %s, charmm_tot: %s\n" % (charmm_vdw,charmm_elec,charmm_tot))

            attribute=attributes.objects.get(attribute_short_name="CHARMM VDW Energy")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=charmm_vdw
            new_min_object_attribute.save()
            resultslog.write("new object attribute: %s with value: %s for pose of the ligand: %s\n" % (attribute.attribute_short_name,new_min_object_attribute.value,ligand.ligand_name))

            attribute=attributes.objects.get(attribute_short_name="CHARMM Electrostatics")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=charmm_elec
            new_min_object_attribute.save()                                                                                                                        

            resultslog.write("new object attribute: %s with value: %s for ligand: %s\n" % (attribute.attribute_short_name,new_min_object_attribute.value,ligand.ligand_name))

            attribute=attributes.objects.get(attribute_short_name="CHARMM Total Energy")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=charmm_tot
            new_min_object_attribute.save()                                                                                                                       

            resultslog.write("new object attribute: %s with value: %s for ligand: %s\n" % (attribute.attribute_short_name,new_min_object_attribute.value,ligand.ligand_name))


            #CHARMM minimization end

            #SEED evaluation                                                                                              
            seed_vdw,seed_elec,seed_tot=common.GetEnergyFromSeedEval(seed_eval_path,posename)
            resultslog.write("seed_vdw: %s, seed_elec: %s, seed_tot: %s\n" % (seed_vdw,seed_elec,seed_tot))

            attribute=attributes.objects.get(attribute_short_name="SEED VDW Score")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=seed_vdw
            new_min_object_attribute.save()
            resultslog.write("new object attribute: %s with value: %s for pose of the ligand: %s\n" % (attribute.attribute_short_name,new_min_object_attribute.value,ligand.ligand_name))
            
            attribute=attributes.objects.get(attribute_short_name="SEED Electrostatics")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=seed_elec
            new_min_object_attribute.save()                                                                                                                                                        
            resultslog.write("new object attribute: %s with value: %s for ligand: %s\n" % (attribute.attribute_short_name,new_min_object_attribute.value,ligand.ligand_name))

            attribute=attributes.objects.get(attribute_short_name="SEED Total Score")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=seed_tot
            new_min_object_attribute.save()                                                                                                                                                       
                                                                                                                      
            resultslog.write("new object attribute: %s with value: %s for ligand: %s\n" % (attribute.attribute_short_name,new_min_object_attribute.value,ligand.ligand_name))
            #SEED evaluation end

            #FFLD evaluation

            ffld_vdw,ffld_elec,ffld_tot=common.GetEnergyFromFFLDEval(ffld_eval_path,posename)
            resultslog.write("ffld_vdw: %s, ffld_elec: %s, ffld_tot: %s\n" % (ffld_vdw,ffld_elec,ffld_tot))

            attribute=attributes.objects.get(attribute_short_name="FFLD VDW Score")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=ffld_vdw
            new_min_object_attribute.save()
            resultslog.write("new object attribute: %s with value: %s for pose of the ligand: %s\n" % (attribute.attribute_short_name,new_min_object_attribute.value,ligand.ligand_name))

            attribute=attributes.objects.get(attribute_short_name="FFLD Electrostatics")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=ffld_elec
            new_min_object_attribute.save()                                                                                                                        

            resultslog.write("new object attribute: %s with value: %s for ligand: %s\n" % (attribute.attribute_short_name,new_min_object_attribute.value,ligand.ligand_name))

            attribute=attributes.objects.get(attribute_short_name="FFLD Total Score")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=ffld_tot
            new_min_object_attribute.save()                                                                                                                        

            resultslog.write("new object attribute: %s with value: %s for ligand: %s\n" % (attribute.attribute_short_name,new_min_object_attribute.value,ligand.ligand_name))

            #FFLD evaluation end

            #create unified scoring attribute
            scores=[]
            scores.extend([float(seed_vdw),float(seed_elec),float(seed_tot),float(ffld_vdw),float(ffld_elec),float(ffld_tot),float(charmm_vdw),float(charmm_elec),float(charmm_tot)])
            resultslog.write("Scores: %s\n" % (scores))
            scores.sort(key=float)
            median_score=scores[4]
            average_score=sum(scores)/len(scores)
            resultslog.write("Median is: %s and Average is: %s\n" % (median_score,average_score))           
 
            attribute=attributes.objects.get(attribute_short_name="DSF Median Score")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=median_score
            new_min_object_attribute.save() 

            attribute=attributes.objects.get(attribute_short_name="DSF Average Score")
            new_min_object_attribute=object_attributes()
            new_min_object_attribute.owner=u
            new_min_object_attribute.attribute_id=attribute.id
            new_min_object_attribute.object_table_name="dd_substrate_poses"
            new_min_object_attribute.object_id=new_min_pose_object.id
            new_min_object_attribute.value=average_score
            new_min_object_attribute.save()           
           
            #create unified scoring attribute end
           
            ######end charmm minimized pose


        #transaction.commit()
        resultslog.write("files processed: %s\n" % (len(glob.glob( os.path.join(path, '*_clus0.mol2')))))
        if len(glob.glob( os.path.join(path, '*_clus0.mol2')))!=0:
            #transaction.commit()
            resultslog.write("committing\n")
            #result_poses = 
            try:
                dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                                      db=settings.DATABASE_NAME, compress=1)
            except:
                pass


            cursor = dba.cursor(MySQLdb.cursors.DictCursor)
            poses_sql="SELECT @count := @count + 1 as rank, selection.* FROM (SELECT p.id, f.id as 'file_id',f.file_name, f.file_location, " \
                      " (select 0+value from dd_analysis_object_attributes where object_table_name='dd_substrate_poses' and object_id=p.id and attribute_id = " \
                      "    (select id from dd_analysis_attributes where attribute_short_name='SEED Total Score')) as seed_score, " \
                      " (select 0+value from dd_analysis_object_attributes where object_table_name='dd_substrate_poses' and object_id=p.id and attribute_id = " \
                      "    (select id from dd_analysis_attributes where attribute_short_name='FFLD Total Score')) as ffld_score, " \
                      " (select 0+value from dd_analysis_object_attributes where object_table_name='dd_substrate_poses' and object_id=p.id and attribute_id = " \
                      "    (select id from dd_analysis_attributes where attribute_short_name='CHARMM Total energy')) as charmm_energy, " \
                      " (select 0+value from dd_analysis_object_attributes where object_table_name='dd_substrate_poses' and object_id=p.id and attribute_id = " \
                      "    (select id from dd_analysis_attributes where attribute_short_name='DSF Median Score')) as median_score, " \
                      " (select 0+value from dd_analysis_object_attributes where object_table_name='dd_substrate_poses' and object_id=p.id and attribute_id = " \
                      "    (select id from dd_analysis_attributes where attribute_short_name='DSF Average Score')) as average_score, " \
                      " (select ligand_name from dd_substrate_ligands where id=p.pose_object_id) as ligand_name " \
                      " from dd_substrate_poses p, dd_infrastructure_files f, (select @count := 0) c " \
                      " where p.source_id = " \
                      " (select id from dd_infrastructure_sources where source_object_table_name='dd_infrastructure_jobs' and source_object_id= %s) " \
                      " and f.id = (select file_id from dd_infrastructure_files_objects where object_table_name ='dd_substrate_poses' and object_id=p.id) " \
                      " order by 0+average_score asc) as selection where selection.file_name like '%%_min.mol2'" % (str(job.id))
        
            resultslog.write("poses_sql2: %s\n" % (poses_sql))
            cursor.execute(poses_sql)
            result_poses = cursor.fetchall()
            cursor.close()
            
             
            #result_files=files.objects.filter(description='Docking By Decomposition Ligand Pose File') \
            #    .extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=" + \
            #    str(request.user.id) + " and object_table_name='dd_infrastructure_jobs' and object_id=" + \
            #    str(job.id) + ")"]).order_by('file_name')
            #transaction.commit()
        else:
            #transaction.rollback()
            resultslog.write("no results\n")

        #result_files=files.objects.filter(description='Docking By Decomposition Ligand Pose File') \
        #        .extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=" + \
        #        str(request.user.id) + " and object_table_name='dd_infrastructure_jobs' and object_id=" + \
        #        str(job.id) + ")"]).order_by('file_name')
    
    #transaction.rollback()
    #########
    

    lesson_ok, dd_ok = checkPermissions(request)
    return render_to_response('html/ddviewjobresults.html', {'result_poses':result_poses, 'job':job,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

"""
#Lets user view file in browser, puts it in iframe
def viewvcsjobfilescontainer(request,fileid):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    #jobfile=Files()
    #fileid=fileid.replace("/","")
    #jobfile = Files.objects.filter(owner=u).order_by("-id")[0]
    #filename=fileid
    #filename=jobfile.location + jobfile.filename
    #command_list=[]
    #command_list.append(filename)
    #return render_to_response('html/debugpageform.html', {'commands':command_list})
    return render_to_response('html/viewvcsjobfilescontainer.html', {'fileid':fileid })
"""
#Lets user view file in browser
def viewDDFile(request,file_id, mimetype = None):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    username = request.user.username

    
    file_id=file_id.replace("/","")
    u = User.objects.get(username=username)
    dfile = files.objects.get(id=file_id)
    filename=dfile.file_location + "/" +dfile.file_name

    try:
       os.stat("%s" % (filename))
    except:
       return HttpResponse("That file doesn't seem to exist. Maybe you deleted it?")

    mimetype = "Content-Type: text/richtext"
    response = HttpResponse(mimetype=mimetype)
    response.write(file("%s" % (filename), "rb").read())
    return response


#Lets user download the actual file
def downloadDDFile(request,file_id, mimetype = None):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    username = request.user.username

    file_id=file_id.replace("/","")
    u = User.objects.get(username=username)
    dfile = files.objects.get(id=file_id)
    filename=dfile.file_location + "/" +dfile.file_name

    if mimetype is None:
        mimetype,encoding = mimetypes.guess_type("%s" % (filename))
    try:
        os.stat("%s" % (filename))
    except:
        return HttpResponse('Oops ... that file no longer exists.')
    response = HttpResponse(mimetype=mimetype)
    response['Content-Disposition'] = 'attachment; filename=%s' % dfile.file_name
    response.write(file(filename, "rb").read())
    return response
    
    
def viewResultPose(request,result_file_id,job_id):

    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')

    try:
        struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")

    try:
        ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
       return HttpResponse("Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")


    
    username = request.user.username
    u = User.objects.get(username=username)
    #ligand_filename="vcsdrawn30"
    job=jobs.objects.get(owner=request.user,job_scheduler_id=job_id.replace("/",""))
    job_folder = charmming_config.user_dd_jobs_home + '/' + username + '/' + 'dd_job_' + str(job.job_owner_index)
    resultfile = files.objects.get(owner=u,id=result_file_id)
    filename=resultfile.file_location + "/" + resultfile.file_name


    viewpose=open("/tmp/viewpose.log",'w')
    viewpose.write("jobid:%s\n" % job_id)
    #Getting the pdb
    
    #conformation_file=files.objects.filter(owner=request.user).extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_infrastructure_jobs' and object_id=%s)" % (request.user.id,job.id)]).extra(where=["id IN (select file_id from dd_infrastructure_files_objects where owner_id=%s and object_table_name='dd_target_protein_conformations')" % (request.user.id)])[0]
    
    pro_segment=Segment.objects.filter(structure=ws.structure,type="pro",is_working="y")[0] #filter[0] is on purpose                                         
    pro_name=pro_segment.name + '-' + str(pro_segment.id)


    proteinfile= "%s/%s.pdb" % (ws.structure.location,pro_name)
    
    ligandfile=filename.replace("/home/schedd/","/charmming/pdbuploads/")
    #proteinfile=conformation_file.file_location+conformation_file.file_name
    proteinfile=proteinfile.replace("/home/schedd/","/charmming/pdbuploads/")
    #viewpose.write("babel -imol2 %s -opdb %s" % (proteinfile,proteinfile.replace(".mol2",".pdb")))
    #os.system("babel -imol2 %s -opdb %s" % (proteinfile,proteinfile.replace(".mol2",".pdb")))
    #proteinfile=conformation_file.file_location.replace("/var/tmp/","/charmming/")+conformation_file.file_name
    proteinfile=proteinfile.replace("/var/tmp/","/charmming/")
    viewpose.write("ligandfile: %s\nproteinfile: %s\n" % (ligandfile,proteinfile))
    

    
    lesson_ok, dd_ok = checkPermissions(request)
    try:
        username = request.user.username
        return render_to_response('html/ddviewresultposejmol.html', {'ligandfile': ligandfile,'proteinfile': proteinfile,'lesson_ok': lesson_ok, 'dd_ok': dd_ok })
    except:
        return HttpResponse("No PDB Uploaded")

def updateDSFLigandsChemSpider(request,selected_set_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    lesson_ok, dd_ok = checkPermissions(request)
    if not dd_ok:
        return render_to_response('html/unauthorized.html')


    log=open("/tmp/dsfligandscs.txt",'w')
    log.write("selecte_set_id: %s\n" % (selected_set_id))
    ligand_info_list=[]
    
    
    #log.write("user_ligands: %s" % (str(user_ligands)))
    public_user=User.objects.get(username='public')
    
    compounds=[]
    rid=chemspipy.asyncfind('chol')
    #compounds.append(chemspipy.find('cholesterol')[0])                                        
    #ligand_info_list.append(ligand_info)
    compounds=chemspipy.getAsyncSearchResult(rid)
    count=len(compounds)

    lesson_ok, dd_ok = checkPermissions(request)        
    return render_to_response('html/dddsfligandschemspi.html', {'ligands':compounds,'count':count,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})
