"""
Small Python script to get all residues we have top/par
files for and map them to the database as a key/value pair
in the following format:
    'residue_name':residue_description
This makes error checking in ligand design page easy
and also allows for checking against this list at any
other point (a-bad/a-good generation).
Since this is a database object it's very mobile and doesn't
require any path information beyond wher ethe toppar files
will be stored.


Vinushka Schalk
Bernard R. Brooks Group
NIH, NHLBI, LCB

07/05/2013
v 0.1 (+/- 0.1)
"""
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
import os

try:
    import settings
except ImportError:
    import sys
    sys.stderr.write("Error: settings.py not in root charmming directory. It appears to be causing an import error.")
    sys.exit(1)

import django.db
from django.db import models
from toppar.models import Residue
import charmming_config

if __name__ == "__main__":
    write_toppar_info()

"""
Uses the path defined in charmming_config.data_home,
adds toppar to the path, then reads in all the top/par files from there
and creates a database object.
"""
def write_toppar_info():
#Since half the files will be ignored most likely, we run a grep first
    logfp = open("/tmp/toppar_load_test.txt", "w")
    toppar_set = set() #So that we don't repeat residues
    try:
        os.chdir(charmming_config.data_home+"/toppar/")
#        p = Popen(['grep','-r','"RESI"','./'],stdout=f) #This will greatly reduce the files we have to analyze
# subprocess hates grep.
        os.system("grep -r RESI ./ > /tmp/toppar.txt")
        logfp.write("grep successful.\n")
        f = open("/tmp/toppar.txt","r") #Open the same file for reading
        for line in f.readlines(): #Read each line (we could do a regex for this I guess?)
        #RESI is in each line so just split by : and see what happens
            if line.startswith("Binary"):
                continue
            filename_end = line.split(":")[0].split("/")[-1].split(".")[-1]
            #This gets the first half of the string, chops the / of the filepath, 
            # splits at the filename dot, then gets the file type
            if filename_end == "out" or filename_end == "history": #Just skip these.
                continue
            if line.split(":")[1].strip().startswith("!"): #Commented line. Skip.
                continue
            tmpstr = line.split(":")[1].split(" ",2) #Grep is always formatted file_path:match
            resname = tmpstr[1] #file_path:RESI CLA for example
            desc = tmpstr[2].split("!")
            if len(desc) < 2:
                desc = ""
            else:
                desc = desc[1].strip() #This strips by comment. If there's nothing, print nothing.
            if resname not in toppar_set:
                toppar_set = toppar_set.union(set([resname]))
                residue = Residue(residue_name=resname, residue_desc=desc) #Database
                residue.save()
                logfp.write(str(toppar_set) + "\n")
        logfp.write("Database saving successful.")
        #Since grep gets every match, we should be ok.
        residue = Residue(residue_name="TIP",residue_desc="") #This is for JSmol since it can't read more than 3-char residue names and it gets confused about what's a solvent otherwise
        f.close()
        os.unlink("/tmp/toppar.txt")
    except Exception as e:
        logfp.write("Something went wrong along the line. This is the exception:\n" + str(e))
    finally:
        logfp.close()

