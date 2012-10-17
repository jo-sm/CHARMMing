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
from pdbinfo.models import PDBFile, PDBFileForm, ParseException, energyParams
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
from pdbinfo.qmmm import makeQChem_tpl, handleLinkAtoms, writeQMheader
from django.contrib.auth.models import User
from django.core import validators 
from django.core.mail import mail_admins
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from account.views import isUserTrustworthy
from pdbinfo.editscripts import generateHTMLScriptEdit
from pdbinfo.aux import checkNterPatch
#dd stuff
from dd_infrastructure.models import projects

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
import common


#Lets user view file in browser, puts it in iframe
def viewProjectsContainer(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    return render_to_response('html/ddviewprojectscontainer.html')

#Lets user see the list of projects they have created
def viewProjects(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    user_projects = projects.objects.filter(owner=request.user)
    return render_to_response('html/ddviewprojects.html', {'user_projects': user_projects})

def viewUpdateProjectInfo(request,projectid):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    username = request.user.username
    owner = User.objects.get(username=username)

    #update or add is taking place
    if ((request.POST['action']=='update') or (request.POST['action']=='addnew')):

	if (len(request.POST['name'].strip())==0): 
	    message="Operation Failed!<br>Project Name Cannot be Blank"
	    if (request.POST['action']=='update'):
                project = projects.objects.filter(owner=request.user,id=request.POST['projectid'])
	        return render_to_response('html/ddprojectinfo.html', {'project':project, 'message':message, 'messagecolor':'Red'})
            else:
		return render_to_response('html/ddnewprojectinfo.html', {'description':request.POST['description'], 'name':request.POST['name'], 'message':message, 'messagecolor':'Red'})
        else:
	    if (request.POST['action']=='update'):
                project = projects.objects.filter(owner=request.user,id=request.POST['projectid'])
                project.project_name=request.POST['name'].strip()
                project.description=request.POST['description'].strip()
                project.save()
	        project_id=request.POST['projectid']
		message="Project Info was successfully updated"
	    elif (request.POST['action']=='addnew'):

                if common.ProjectExists(request.POST['name'],request.user)==True:
		    return render_to_response('html/ddnewprojectinfo.html', {'description':request.POST['description'], 'name':request.POST['name'], 'message':'Operation Failed!<br>Duplicate Project Found.', 'messagecolor':'Red'})

	        project=projects()
                project.owner=owner
                project.date_created=datetime.datetime.now().isoformat(' ')
	        project.project_name=request.POST['name'].strip()
                project.description=request.POST['description'].strip()
	        project.save()
	        
	        added_project = projects.objects.filter(owner=request.user).order_by('-id')[0]
                project_id = added_project.id
		message="Project was successfully created"

	project = projects.objects.filter(owner=request.user,id=project_id)
	return render_to_response('html/ddprojectinfo.html', {'project':project, 'message':message, 'messagecolor':'Green'})

    #no action, just display blank or populated form
    else: 
	if request.POST['projectid']!='0':
	    project = projects.objects.filter(owner=request.user,id=request.POST['projectid'])
            return render_to_response('html/ddprojectinfo.html', {'set':set, 'message':'', 'messagecolor':'Red'})
	else:
	    return render_to_response('html/ddnewprojectinfo.html', {'description':'', 'message':'', 'messagecolor':'Red'})

def projectsaddnew(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    
    return render_to_response('html/ddprojectsaddnew.html', {})

#This is for changing the currently selected project
#on the manage projects page
def switchProjects(request,switch_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    try:
        oldproject =  projects.objects.filter(owner=request.user,selected='y')[0]
        oldproject.selected = ''
        oldproject.save()
    except:
        pass
    newproject = projects.objects.filter(owner=request.user,id=switch_id)[0] 
    newproject.selected = 'y'
    newproject.save()
    return render_to_response('html/ddswitchproject.html',{'oldproject':oldproject,'newproject':newproject})

def projectDetails(request,project_id):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    
    try:
       dba = MySQLdb.connect("localhost", user="charmming", passwd="charmming", db="charmming", compress=1)

    except:
       pass

    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    conformations_sql = "select fo.file_id as file_id, pc.id as conformation_id, pc.conformation_name as conformation_name, s.id as structureid, p.name as protein_name,  "\
   	       "from dd_infrastructure_files_objects fo left join dd_target_protein_conformations pc on fo.object_id=pc.id "\
	       "left join dd_target_proteins p on p.id=pc.protein_id " \
	       "where fo.object_table_name = 'dd_target_protein_conformations' and pc.owner_id=%s and fo.owner_id=%s" % (request.user.id,request.user.id)
    cursor.execute(conformations_sql)
    rows = cursor.fetchall()
    grids=[]
    for row in rows:
	structuregrid=common.StructureGridFileObj()
	structuregrid.grid_id=row['gridid']
	structuregrid.file_id=row['fileid']
        structuregrid.grid_code=row['gridcode']      
        structuregrid.structure_name=row['structurename']
        structuregrid.structure_id=row['structureid']
	structuregrid.protein_type=row['shortproteintype']
        grids.append(structuregrid)   
  
    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    #setssql = "select '0' as setid, gs.name as setname, (select count(distinct id) from vcs_structuregrids where owner_id=%s) as gridcount, " \
    #	      "(select count(distinct id) from vcs_structures where owner_id=%s) as structurecount from vcs_gridsets gs where gs.name='All Available' union " \
    setssql = "select '0' as setid, 'All Available' as setname, (select count(distinct id) from vcs_structuregrids where owner_id=%s) as gridcount, " \
	      "(select count(distinct id) from vcs_structures where owner_id=%s) as structurecount union " \
	      "select gs.id as setid, gs.name as setname, (select count(distinct grid_id) from vcs_gridsetgrids where gridset_id=gs.id and owner_id=%s) as gridcount, " \
	      "(select count(distinct structure_id) from vcs_structureconformations where id in (select structure_conformation_id from vcs_structuregrids where owner_id=%s and id in " \
	      "(select grid_id from vcs_gridsetgrids where owner_id=%s and gridset_id=gs.id))) as structurecount " \
	      "from vcs_gridsets gs where gs.owner_id=%s and gs.name<>'All Available'" % (request.user.id,request.user.id,request.user.id,request.user.id,request.user.id,request.user.id)
    
    cursor.execute(setssql)
    rows = cursor.fetchall()
    sets=[]

    for row in rows:
	set=common.GridSetObj()
	set.id=row['setid']
	set.name=row['setname']
        set.grid_count=row['gridcount']
	set.structure_count=row['structurecount']
        sets.append(set) 

    set_id=set_id.replace("/","")
    

    proteintypes=[]
    proteintypes=ProteinTypes.objects.all().order_by('short_name')

    if set_id=="0":
	return render_to_response('html/ddprojectdetails.html', {'grids':grids, 'sql':gridssql, 'sets':sets, 'proteintypes':proteintypes})

    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    editsetsql = "select id as setid, name as setname, description, (select count(id) from vcs_gridsetgrids where gridset_id=%s) as gridcount, "\
		"(select count(distinct id) from vcs_structures where id in "\
		"(select distinct structure_id from vcs_structureconformations where id in "\
		"(select distinct structure_conformation_id from vcs_structuregrids where id in "\
		"(select grid_id from vcs_gridsetgrids where gridset_id=%s)))) as structurecount "\
	        "from vcs_gridsets where id=%s " % (set_id,set_id,set_id)

    cursor.execute(editsetsql)
    rows = cursor.fetchall()

    for row in rows:
	editset=common.GridSetObj()
	editset.id=row['setid']
	editset.name=row['setname']
	editset.description=row['description']
   
    cursor.close()
    dba.close()

    return render_to_response('html/ddprojectdetails.html', {'grids':grids, 'sql':gridssql, 'sets':sets, 'proteintypes':proteintypes, 'editset':editset})

def projectAvailableConformations(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    
    try:
       dba = MySQLdb.connect("localhost", user="charmming", passwd="charmming", db="charmming", compress=1)

    except:
       pass

    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    if (request.POST['ptid']!='0') and (request.POST['ptid']!=''):
	gridssql = "select fo.file_id as fileid, sg.id as gridid, sg.code as gridcode, s.id as structureid, s.name as structurename, pt.name as shortproteintype "\
   	       "from vcs_filesobjects fo left join vcs_structuregrids sg on fo.object_id=sg.id left join vcs_structureconformations sc on sc.id=sg.structure_conformation_id "\
	       "left join vcs_structures s on s.id=sc.structure_id " \
	       "left join vcs_proteintypes pt on s.protein_type_id=pt.id where fo.object_table_id in (select id from vcs_objecttables where name = 'vcs_structuregrids') "\
	       "and s.protein_type_id=%s and sg.owner_id=%s" % (request.POST['ptid'],request.user.id)
    	
    else:
	
	gridssql = "select fo.file_id as fileid, sg.id as gridid, sg.code as gridcode, s.id as structureid, s.name as structurename, pt.name as shortproteintype "\
   	       "from vcs_filesobjects fo left join vcs_structuregrids sg on fo.object_id=sg.id left join vcs_structureconformations sc on sc.id=sg.structure_conformation_id "\
	       "left join vcs_structures s on s.id=sc.structure_id " \
	       "left join vcs_proteintypes pt on s.protein_type_id=pt.id where fo.object_table_id in (select id from vcs_objecttables where name = 'vcs_structuregrids') "\
	       "and sg.owner_id=%s" % (request.user.id)

    if (request.POST['setid']!='0') and (request.POST['setid']!=''):
	setcriteria=" and sg.id in (select grid_id from vcs_gridsetgrids where owner_id=%s and gridset_id=%s)" % (request.user.id,request.POST['setid'])
	gridssql=gridssql + setcriteria

    cursor.execute(gridssql)
    rows = cursor.fetchall()
    grids=[]
    for row in rows:
	structuregrid=common.StructureGridFileObj()
	structuregrid.grid_id=row['gridid']
	structuregrid.file_id=row['fileid']
        structuregrid.grid_code=row['gridcode']      
        structuregrid.structure_name=row['structurename']
        structuregrid.structure_id=row['structureid']
	structuregrid.protein_type=row['shortproteintype']
        grids.append(structuregrid)     
    
    cursor.close()
    dba.close()

    proteintypes=[]
    proteintypes=ProteinTypes.objects.all().order_by('short_name')
    
    return render_to_response('html/dockgrids.html', {'grids':grids, 'proteintypes':proteintypes})

def projectConformations(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    
    try:
       dba = MySQLdb.connect("localhost", user="charmming", passwd="charmming", db="charmming", compress=1)

    except:
       pass


    set_id=request.POST['setid']
    managesetlog = open("/tmp/manageset.log",'w')
    managesetlog.write(request.POST['removedids'])
    addedids=request.POST['addedids']
    removedids=request.POST['removedids']
    if (addedids!=""):
	cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
	grididssql="select distinct(object_id) as gridid from vcs_filesobjects where file_id in (%s) and object_table_id = (select id from vcs_objecttables where name = 'vcs_structuregrids')" % (addedids)
        cursor.execute(grididssql)
        rows = cursor.fetchall()
	for row in rows:
	    gridsetgrid=GridSetGrids()
	    gridsetgrid.owner_id=request.user.id
	    gridsetgrid.gridset_id=set_id
	    gridsetgrid.grid_id=row['gridid']
	    gridsetgrid.save()
        
    
    if (removedids!=""):
	cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
	removegridssql="delete from vcs_gridsetgrids where gridset_id=%s and grid_id in " \
	               "(select object_id from vcs_filesobjects where file_id in (%s) and object_table_id = (select id from vcs_objecttables where name = 'vcs_structuregrids'))" % (set_id,removedids)
	managesetlog.write(removegridssql)
        cursor.execute(removegridssql)


    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 
    if request.POST['setid']!='0':
	gridssql = "select fo.file_id as fileid, sg.id as gridid, sg.code as gridcode, s.id as structureid, s.name as structurename, pt.name as shortproteintype "\
   	       "from vcs_filesobjects fo left join vcs_structuregrids sg on fo.object_id=sg.id left join vcs_structureconformations sc on sc.id=sg.structure_conformation_id "\
	       "left join vcs_structures s on s.id=sc.structure_id " \
	       "left join vcs_proteintypes pt on s.protein_type_id=pt.id where fo.object_table_id in (select id from vcs_objecttables where name = 'vcs_structuregrids') "\
	       " and sg.owner_id=%s " % (request.user.id)
    	
        editsetsql = gridssql + " and sg.id in (select grid_id from vcs_gridsetgrids where owner_id=%s and gridset_id=%s)" % (request.user.id,set_id)

    cursor.execute(editsetsql)
    rows = cursor.fetchall()
    grids=[]
    for row in rows:
	structuregrid=common.StructureGridFileObj()
	structuregrid.grid_id=row['gridid']
	structuregrid.file_id=row['fileid']
        structuregrid.grid_code=row['gridcode']      
        structuregrid.structure_name=row['structurename']
        structuregrid.structure_id=row['structureid']
	structuregrid.protein_type=row['shortproteintype']
        grids.append(structuregrid)     
    
    cursor.close()
    cursor = dba.cursor(MySQLdb.cursors.DictCursor) 

    editsetsql = "select id as setid, name as setname, description, (select count(id) from vcs_gridsetgrids where gridset_id=%s) as gridcount, "\
		"(select count(distinct id) from vcs_structures where id in "\
		"(select distinct structure_id from vcs_structureconformations where id in "\
		"(select distinct structure_conformation_id from vcs_structuregrids where id in "\
		"(select grid_id from vcs_gridsetgrids where gridset_id=%s)))) as structurecount "\
	        "from vcs_gridsets where id=%s " % (set_id,set_id,set_id)
    #command_list=[]
    #command_list.append(editsetsql)
    #return render_to_response('html/debugpageform.html', {'commands':command_list})
    cursor.execute(editsetsql)
    rows = cursor.fetchall()

    for row in rows:
	editset=common.GridSetObj()
	editset.id=row['setid']
	editset.name=row['setname']
	editset.grid_count=row['gridcount']
	editset.structure_count=row['structurecount']
   

    cursor.close()
    dba.close()

    
    return render_to_response('html/managesets_setgrids.html', {'grids':grids, 'set':editset, 'addedids':addedids, 'setid':set_id})

