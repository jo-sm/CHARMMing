import commands, datetime, sys, re, os, glob, shutil
import charmming_config, input
import smtplib  
import cPickle, copy, traceback, socket
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
    ligandsetid=0
    name=""
    description=""
    file_objects_list=[]
    filecount=0
    files=[]
    file=""
    file_id=0

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


def makeMatch(tempname):
    """                                                                                                                                                   
    Uses the match program from the Charlie Brooks group to                                                                                               
    try to build topology and parameter files.                                                                                                            
    """

    logfp = open('/tmp/dd_match.txt', 'w')

    os.putenv("PerlChemistry","%s/MATCH_RELEASE/PerlChemistry" % charmming_config.data_home)
    os.putenv("MATCH","%s/MATCH_RELEASE/MATCH" % charmming_config.data_home)
    os.chdir("/tmp")

    #self.rtf_list = '%s/toppar/top_all27_prot_na.rtf' % charmming_config.data_home
    #self.prm_list = '%s/toppar/par_all27_prot_na.prm' % charmming_config.data_home
    #for badRes in badResList:
    exe_line = '%s/MATCH_RELEASE/MATCH/scripts/MATCH.pl %s.mol2 > /tmp/%s.txt' % (charmming_config.data_home,tempname,tempname)
    logfp.write('Xcute: %s\n' % exe_line)
    status, output = commands.getstatusoutput(exe_line)
    if status != 0:
        logfp.write("sorry ... nonzero status :-(\n")
        logfp.write(output + '\n')
        logfp.close()
        return -1
    else:
        os.chdir("/tmp")
        #os.system("mv %s.rtf %s_temp.rtf" % (tempname,tempname))
        os.system("awk '/RESI/{$2=\"ULIG\"; print;next}{print}' %s.rtf > %s_temp.rtf" % (tempname,tempname))
        os.system("mv %s.prm %s_temp.prm" % (tempname,tempname)) 
        os.system("grep -v \"MASS\" %s_temp.rtf > %s.rtf" % (tempname,tempname))
        os.system("awk '/NONBONDED/{exit;}{print;}' %s_temp.prm > %s.prm" % (tempname,tempname))
        os.system("echo 'read rtf card append' > %s.str" % (tempname))
        os.system("cat %s.rtf >> %s.str" % (tempname,tempname))
        os.system("echo 'read para card flex append' >> %s.str" % (tempname))
        os.system("cat %s.prm >> %s.str" % (tempname,tempname))
        os.system("printf 'END\nRETURN' >> %s.str" % (tempname))
 
    logfp.write("OK!\n")

        # ToDo: add the rtf/param files that have been generated to a master topology                                                                     
        # and parameter files for the entire segment.                                                                                                     
     #self.rtf_list += ' %s-badres-h-%s.rtf' % (self.name,badRes)
     #self.prm_list += ' %s-badres-h-%s.prm' % (self.name,badRes)

    logfp.close()
    return 0


def makeAntechamber(tempname):

    logfp = open('/tmp/dd_ac.log', 'w')
    logfp.write('In makeAntechamber.\n')

    os.putenv("AMBERHOME", charmming_config.amber_home)
    os.chdir("/tmp")

    for badRes in badResList:

        basefile = '%s' % (tempname)
        mol2file = basefile + '.mol2'
        cinpfile = basefile + '.inp'
        rtffile = basefile + '.rtf'
        prmfile = basefile + '.prm'

        cmd = charmming_config.amber_home + "/bin/antechamber -an n -rn " + badRes + \
              " -fi mol2 -i " + mol2file + " -fo charmm -o " + basefile
        logfp.write(cmd + '\n')
        status = os.system(cmd)
        if status != 0:
            # debug purposes only                                                                                                                         
            # raise("Antechamber screwed up!")                                                                                                            
            return -1

        try:
            # see if output file exists                                                                                                                   
            os.stat(cinpfile)
            os.stat(rtffile)
            os.stat(prmfile)
        except:
            logfp.write('No CHARMM input file produced!\n')
            return -1

        # set the newly generated RTF/PRM to be used.                                                                                                     
        if self.rtf_list:
            self.rtf_list += " " + rtffile
        else:
            self.rtf_list = rtffile
        if self.prm_list:
           self.prm_list += " " + prmfile
        else:
           self.prm_list = prmfile

    self.save()
    logfp.close()
    return 0

def makeCGenFF(ws,tempname):
    """                                                                                                                                                   
    Connects to dogmans.umaryland.edu to build topology and                                                                                               
    parameter files using CGenFF                                                                                                                          
    """
    logfp = open('/tmp/dd_makecgenff.txt', 'w')
    logfp.write('Try to use CGenFF\n')
    #self.rtf_list = '%s/toppar/top_all36_cgenff.rtf' % charmming_config.data_home
    #self.prm_list = '%s/toppar/par_all36_cgenff.prm' % charmming_config.data_home

    header = '# USER_IP 165.112.184.52 USER_LOGIN %s\n\n' % ws.structure.owner.username
    #for badRes in badResList:
    
    fp = open('/tmp/%s.mol2' % (tempname), 'r')
    payload = fp.read()
    content = header + payload
    fp.close()

    recvbuff = ''
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((charmming_config.cgenff_host, charmming_config.cgenff_port))
        s.sendall(content)
        s.shutdown(socket.SHUT_WR)
        while True:
            data = s.recv(1024)
            if not data:
                break
            recvbuff += data

    except Exception, e:
        logfp.write('crap something went wrong::\n')
        traceback.print_exc(file=logfp)
        logfp.close()
        return -1

    fp = open('/tmp/%s-dogmans.txt' % (tempname), 'w')
    fp.write(recvbuff)
    fp.close()
    lines=open('/tmp/%s-dogmans.txt' % (tempname)).readlines()
    #open('/tmp/%s.str' % (tempname),'w').writelines(lines[1:])
    logfp.write('Dumped CGenFF result for analysis.\n')

    # parse file returned from dogmans, and make sure that it has both topology and                                                                   
    # parameters.                                                                                                                                     
    rtffound = False
    prmfound = False
    inrtf = False
    inprm = False
    rtffp = open('/tmp/%s-dogmans.rtf' % (tempname), 'w')
    prmfp = open('/tmp/%s.prm' % (tempname), 'w')
    for line in recvbuff.split('\n'):
        if inrtf:
            rtffp.write(line + '\n')
        if inprm:
            prmfp.write(line + '\n')

        line = line.strip()
        if line == 'end' or line == 'END':
            inrtf = False
            inprm = False
        if line.startswith('read rtf'):
            rtffound = True
            inrtf = True
        if line.startswith('read param'):
            prmfound = True
            inprm = True

    rtffp.close()
    prmfp.close()

    if not (rtffound and prmfound):
        logfp.write('Aw crap ... did not find both topology and parameters.\n')
        logfp.close()
        return -1
    else:
        os.chdir("/tmp")
        os.system("awk '/RESI/{$2=\"ULIG\"; print;next}{print}' dd_user_ligand-dogmans.rtf > dd_user_ligand.rtf")
        logfp.write("awk '/RESI/{$2=\"ULIG\"; print;next}{print}' dd_user_ligand-dogmans.rtf > dd_user_ligand.rtf\n")
        os.system("awk '/Toppar/,0' dd_user_ligand-dogmans.txt > dd_user_ligand.str")



    #self.rtf_list += ' %s-%s-dogmans.rtf' % (self.name,badRes)
    #self.prm_list += ' %s-%s-dogmans.prm' % (self.name,badRes)


    logfp.write('All looks OK\n')
    logfp.close()
    return 0

def checkDAIM(tempname):

    logfp = open('/tmp/dd_daim.txt', 'w')

    os.chdir("/tmp")

    #self.rtf_list = '%s/toppar/top_all27_prot_na.rtf' % charmming_config.data_home                                                                                                                          
    #self.prm_list = '%s/toppar/par_all27_prot_na.prm' % charmming_config.data_home                                                                                                                          
    #for badRes in badResList:                                                                                                                                                                               
    os.system("cp %s ." % (charmming_config.daim_param))
    os.system("cp %s ." % (charmming_config.daim_prop))
    os.system("cp %s ." % (charmming_config.daim_weight))

    exe_line = '%s --paramfile %s  %s.mol2 > /tmp/%s_daim.txt' % (charmming_config.daim_exe,charmming_config.daim_param,tempname,tempname)
    logfp.write('Xcute: %s\n' % exe_line)
    status, output = commands.getstatusoutput(exe_line)
    if status != 0:
        logfp.write("sorry ... nonzero status :-(\n")
        logfp.write(output + '\n')
        logfp.close()
        return -1
    else:
        os.chdir("/tmp")
        if "WARNING" in open("/tmp/daim.log").read():
            return -1

    logfp.write("OK!\n")
    logfp.close()
    return 0
