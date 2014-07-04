#!/usr/bin/python
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

from django.core.management import setup_environ
#from charmming import settings

from httplib import HTTPConnection
from django import forms
from django.template.loader import get_template




import sys, os, string, getopt, commands, glob, common
import common
#dd stuff
#
from models import ligands 
from charmming.dd_infrastructure.models import files,files_objects,file_types, sources
###


# main routine
if __name__ == '__main__':

    log=open("/tmp/ligand_batch_upload.log",'w')
    location="/var/tmp/dd/ligands/public/zinc_initial_test_compounds/"
    u = User.objects.get(id=2)
    for file in glob.glob(os.path.join(location, '*.mol2')):
        file_obj=open(file,'r')
        mol_title=getMOleculeTitle(file_obj)
        newligand = ligands()
        newfile = files()
        newfileobject= files_objects()
        try: 
            newligand=ligands.objects.filter(owner=2,name=mol_title)[0]
        except:
            newligand.owner=u
            newligand.ligand_name=mol_title
            newligand.ligand_owner_index=str(common.getNewLigandOwnerIndex(u))
            newligand.description="System preloaded drug-like compound"
            source=sources.objects.get(source_name='Zinc Database')
            newligand.source=source
            print("saving ligand with owner id %s and name %s\n" % (newligand.ligand_owner_index,newligand.ligand_name))
            #newligand.save()
       
            new_location=charmming_config.user_dd_ligands_home + '/' + 'public' + '/' + 'ligand_' + str(newligand.ligand_owner_index) + '/'
            filename='ligand_' + str(newligand.ligand_owner_index) + '_' + str(newligand.ligand_owner_index)+'.mol2'
            ligand_file_type=file_types.objects.filter(file_type_name="Ligand Structure File")[0]
            newfile.owner=u
            newfile.file_name=filename
            newfile.file_location=location
            newfile.file_type=ligand_file_type
            newfile.description="Structure File for Ligand %s, %s" % (newligand.ligand_name,newligand.description)
            print("saving file with location %s and name %s\n" % (newfile.file_location,newfile.file_name))
            #newfile.save()
            
            newfileobject.owner=u
            newfileobject.file=newfile
            newfileobject.object_table_name='dd_substrate_ligands'
            newfileobject.object_id=newligand.id
            print("saving new fileobject\n")
            #newfileobject.save()
            
            os.system("mkdir " + charmming_config.user_dd_ligands_home + '/' + "public" + '/')
            os.system("cp %s %s" % (file,new_location))
            print("mkdir " + charmming_config.user_dd_ligands_home + '/' + "public" + '/')
            print("cp %s %s" % (file,new_location))
        
    log=close()


def getMoleculeTitle(file):
    while 1:
        line = file.readline()
        if not line:
            break
        if line.strip() == "@<TRIPOS>MOLECULE":
            titleline.readline()
            title=titleline.strip()
            break
    return title
