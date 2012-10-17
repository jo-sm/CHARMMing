import smtplib  
import os  
from email.mime.multipart import MIMEMultipart   
from email.mime.base import MIMEBase  
from email.mime.text import MIMEText  
from email.utils import COMMASPACE, formatdate  
from email import encoders  
#dd stuff
from dd_infrastructure.models import projects, files, files_objects
from dd_substrate.models import ligands, poses
##
from django.contrib.auth.models import User
from django.db import transaction
import time
import datetime
#from charmming import charmming_config
import MySQLdb
import MySQLdb.cursors

def getNewLigandOwnerIndex(user):
    next_ligand_id=1
    next_conformation_id=1
    next_ligand_owner_index=1
    rllog=open("/tmp/rl.log",'w')
     
    try:
        
        last_ligand = ligands.objects.filter(owner=user).order_by('-id')[0]
        next_ligand_owner_index=int(last_ligand.ligand_owner_index)+1

    except:
    	next_ligand_owner_index=1

    return next_ligand_owner_index


def getNewLigandFilename(user):

    return "ligand_%s" % (str(getNewligandOwnerIndex(user)))
    
    
        
class LigandInfo():
    
    id=0
    name=""
    description=""
    file_objects_list=[]
    filecount=0
    files=[]
    file=""

class LigandSetObj():
    id=0
    name=""
    description=""
    ligand_count=0
    datecreated=""
    
def LigandSetExists(setname, user):
    try:
        sets=ligand_sets.objects.filter(name=setname,owner=user)
        if sets:
	        return True
        else:
            return False
    except:
        return False
    
@transaction.commit_manually
def DeletePose(user,delete_id):
    public_user=User.objects.get(username='public')
    transaction.commit()
    deleteposelog=open("/tmp/deletelig.log",'a')
    deleteposelog.write("starting pose deletion\n")
    #make sure request data isn't malicious                                                                                                                                                    
    #input.checkRequestData(request)                                                                                                                                                              
    if delete_id<>0:
        deleteposelog.write("deleting pose with id %s\n" % (delete_id))
        try:
            delete_pose=poses.objects.get(id=delete_id)
            deleteposelog.write("found pose: %s, %s\n" % (delete_pose.id,delete_pose.pose_name))
            pose_files_objects=files_objects.objects.filter(owner__in=[user,public_user],object_table_name='dd_substrate_poses',object_id=delete_pose.id)
            for pose_file_object in pose_files_objects:
                deleteposelog.write("    found pose file object: %s\n" % (pose_file_object.id))
                try:
                    pose_files=files.objects.filter(owner__in=[user,public_user],id=pose_file_object.file_id)
                    for pose_file in pose_files:
                        deleteposelog.write("        found pose file: %s, %s\n" % (pose_file.id, pose_file.file_name))
                        try:
                            pose_associations=files_objects.objects.filter(owner__in=[user,public_user], file=pose_file).exclude(id=pose_file_object.id).delete()
                            #transaction.commit()
                        except:
                            transaction.rollback()
                            return False
                    pose_files.delete()
                    #transaction.commit()
                except:
                    transaction.rollback()
                    return False
            pose_files_objects.delete()
            delete_pose.delete()
            transaction.commit()
        except:
            transaction.rollback()   
            return False
                #pose_associations=files_objects.objects.filter(owner__in=[user,public_user], file=pose_file).exclude(id=pose_file_object.id)
                #for pose_association in pose_associations:
                #    deleteposelog.write("            found pose association: %s, %s, %s\n" % (pose_association.id, pose_association.object_table_name, pose_association.object_id))
    else:
        return False

    return True
