import smtplib  
import os  
from email.mime.multipart import MIMEMultipart   
from email.mime.base import MIMEBase  
from email.mime.text import MIMEText  
from email.utils import COMMASPACE, formatdate  
from email import encoders  
#dd stuff
from dd_infrastructure.models import projects, files, files_objects
from dd_target.models import proteins, protein_conformations
##
import time
import datetime
import charmming_config
import MySQLdb
import MySQLdb.cursors

def getNewProteinOwnerIndex(user):
    next_target_id=1
    next_conformation_id=1
    next_target_owner_index=1
    #uploadtargetlog1=open("/tmp/targetupload1.log",'w')
     
    try:

        last_target = proteins.objects.filter(owner=user).order_by('-id')[0]
        next_target_owner_index=int(last_target.protein_owner_index)+1

    except:
	next_target_owner_index=1

    return next_target_owner_index


def getNewProteinFilename(user):

    return "target_%s" % (str(getNewProteinOwnerIndex(user)))

def getNewConformationProteinIndex(protein):
    next_target_id=1
    next_conformation_id=1
    #uploadtargetlog1=open("/tmp/targetupload1.log",'w')
     
    try:

        last_conformation = protein_conformations.objects.filter(protein=protein).order_by('-id')[0]
        next_conformation_protein_index=int(last_conformation.conformation_protein_index)+1

    except:
        next_conformation_protein_index=1

    return next_conformation_protein_index


def getNewConformationFilename(protein):

    return "target_%s_%s" % (str(getNewProteinOwnerIndex(user)),str(getNewConformationProteinIndex(protein)))

class objProteinConformation():

    #file_id=0
    conformation_id=0
    conformation_name=""
    protein_pdb_code=""
    files=[]
    #conformation=protein_conformations()
   
  
class TargetConformationInfo():
    
    id=0
    name=""
    description=""
    file_objects_list=[]
    filecount=0
    target_pdb_code=""
    conformations=[]    


