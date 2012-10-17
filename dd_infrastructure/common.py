import smtplib  
import os
from django.contrib.auth.models import User  
from email.mime.multipart import MIMEMultipart   
from email.mime.base import MIMEBase  
from email.mime.text import MIMEText  
from email.utils import COMMASPACE, formatdate  
from email import encoders  
#dd stuff
from dd_infrastructure.models import projects, files, files_objects, jobs
from dd_target.models import protein_conformations, proteins
###
import time
import datetime
import charmming_config
import MySQLdb
import MySQLdb.cursors

def GetEnergyFromFLEAMol(filename):
    mol2file = open(filename, 'r')
    energy = "9999"
    for line in mol2file:
        if line.find(":")>0:
            energy=line[line.find(":")+1:].strip()
            return energy

    return energy

def ProjectExists(projectname, user):
    try:
        projects=projects.objects.filter(name=projectname,owner=user)
        if projects:
	        return True
        else:
            return False
    except:
        return False

def getNewJobOwnerIndex(user):

    try:
        last_job = jobs.objects.filter(owner=user).order_by('-id')[0]
        next_job_owner_index=int(last_job.job_owner_index)+1

    except:
        next_job_owner_index=1

    return next_job_owner_index
       
class objProjectProteinConformation():

    #file_id=0
    conformation_id=0
    project_conformation_id=0
    conformation_name=""
    protein_pdb_code=""
    #file=files()
    #conformation=protein_conformations()

class objDDJob():
    id=0
    scheduler_id=0
    owner_index=0
    description=""
    status=""
    queued=datetime.datetime.now()
    waited=datetime.datetime.now()
    runtime=datetime.datetime.now()
    totaltime=datetime.datetime.now() 
    output_file=""
    target_conformations_file_list=[]
    ligands_file_list=[]
    results_file_list=[]
    
    
def mergeLists(list1,list2):

    for item in list2:
        list1.append(item)
      
    return list1

class objFile():
    fullpath=""
    name=""
