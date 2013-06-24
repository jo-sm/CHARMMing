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
from django.db import models
from structure.models import Task, WorkingFile, WorkingStructure
from pychm.io.pdb import PDBFile
from django.http import HttpResponse
import os, copy, shutil
import subprocess
from subprocess import Popen
import cPickle
import charmming_config

class mutateTask(Task):
    MutResi = models.PositiveIntegerField(default=0)
    MutResName = models.CharField(max_length=3,null=True)
    NewResn = models.CharField(max_length=3,null=True) #ALA, SYN, GLY, &c.
    MutSegi = models.CharField(max_length=50,null=True)
    MutFile = models.CharField(max_length=500,null=True)
    localpickle = models.CharField(max_length=500,null=True)

    def finish(self):
        logfp = open("/tmp/logmut.txt","a+")
        #overrides basic finish, tests if the job succeeded
        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier

        #create WorkingFile for the input file...
        wfinp = WorkingFile()
        wfinp.task = self
        wfinp.path = loc + '/' + bnm + '-mutate.inp'
        wfinp.canonPath = wfinp.path
        wfinp.type = 'inp'
        wfinp.description = 'mutation script input'

        # Check if an output file was created and if so create
        # a WorkingFile for it.
        try:
            os.stat(loc + '/' + bnm + '-mutate.out')
        except:
            self.status = 'F'
            logfp.write("Could not create inp/out files.\n")
            return
        logfp.write("Created inp/out files.\n")

        wfout = WorkingFile()
        wfout.task = self
        wfout.path = loc + '/' + bnm + '-mutate.out'
        wfout.canonPath = wfout.path
        wfout.type = 'out'
        wfout.description = 'mutation script output'
        wfout.save()

        if self.status == 'F':
            logfp.write("Did not write output files.\n")
            return
        logfp.write("Wrote output files.\n")

        # check and make sure that the output PSF/CRD were 
        # created
        try:
            os.stat(loc + '/' + bnm + '-mutation.crd')
        except:
            self.status = 'F'
            self.save()
            logfp.write("PSF/CRD were not created.\n")
            return
        logfp.write("PSF/CRD were created.\n")

        # create Working files for PDB, CRD, and PSF.
        wf = WorkingFile()
        wf.task = self
        wf.path = loc + '/' + bnm + '-mutation.crd'
        wf.canonPath = wf.path
        wf.type = 'crd'
        wf.description = 'mutated structure'
        wf.pdbkey = 'mut_' + self.workstruct.identifier
        wf.save()

        wfpsf = WorkingFile()
        wfpsf.task = self
        wfpsf.path = loc + '/' + bnm + '-mutation.psf'
        wfpsf.canonPath = wfpsf.path
        wfpsf.type = 'psf'
        wfpsf.description = 'mutated structure'
        wfpsf.save()

        wfpdb = WorkingFile()
        wfpdb.task = self
        wfpdb.path = loc + '/' + bnm + '-mutation.pdb'
        wfpdb.canonPath = wfpdb.path
        wfpdb.type = 'pdb'
        wfpdb.description = 'mutated structure'
        wfpdb.save()
        #At this point we can assume the script ran successfully so we can wipe tmp.crd
        try:
            os.remove(loc + "/tmp.crd")
        except:
            pass 
        #OR at least attempt to.
        #This is where this diverges from normal task finishes
        #We need to reprocess the PDBFile object, create a new pickle...
        pdbloc = bnm + "mutation.pdb"
        logfp.write(pdbloc + "\n")
        try:
            os.chdir(loc) #We use Popen because we need to use > and we don't want to use shell=True (prevents need for input sanitization)
            f = open("tmp.tmp","w")
            p = Popen(["stride","-o",pdbloc],stdout=f) #Identifier (bnm) assumed to already be "1vom-mt1" for example.
            f.close()
#            os.system("stride -o " + pdbloc + " > tmp.tmp")
            logfp.write("STRIDE ran successfully.\n")
            f = open("tmp.ss","w")
            p = Popen([charmming_config.data_home + "/stride2pdb","tmp.tmp"],stdout=f)
            f.close()
#            os.system("stride2pdb tmp.tmp > tmp.ss")
            logfp.write("stride2pdb conversion complete.\n")
            f = open("temp-pdb.pdb","w")
            p = Popen(["cat","tmp.ss",pdbloc],stdout=f)
            f.close()
#            os.system("cat tmp.ss " + pdbloc + " > temp-pdb.pdb")
            logfp.write("PDB/SS concatenation complete.\n")
            os.remove("tmp.ss")
            os.remove("tmp.tmp")
            logfp.write("temp files cleared.\n")
        except Exception as ex:
            logfp.write(str(ex) + "\n")
            self.status = "F"
            self.save()
            return
        try:
            oldpickle = open(self.localpickle,'r')
        except:
            oldpickle = open(self.workstruct.structure.pickle,'r')
        pdb = cPickle.load(oldpickle)
        metadata = pdb.get_metaData() #This is where things get fun
        #Now the PDB file has all the atoms stored in it, what we need now is to use the tmp PDB with the HELIX/SHEET info
        #take the old metadata and write it in, create a PDBFile object to store the new metadata, create the new pickle,
        #then wipe that PDB file since we only use the one that has only atoms.
        #TODO: Add SEQRES update function
        #TODO: Add CONECT update function
        #the metadata's REMARK lines and stuff might remain the same, the issue comes in SEQRES, HELIX, SHEET and a few others.
        newPDB = open("temp-pdb.pdb","a+")    #Write them at the end because PDBFile doesn't care where they end up.
        for key in metadata.keys():
            if key not in (['helix', 'sheet', 'seqres']):
                #Write them to your new PDB file...
                linelength = len(metadata[key][0]) #This way we can even out the spacing
                for item in metadata[key]:
                    newPDB.write(key.upper() + (((linelength - len(item)) + 4) * " ") + item + "\n") #4 spaces or less depending on length of the record
        ##TODO: Ask Frank about which of these keys are dangerous for mutation
        newPDB.close()
        logfp.write("PDB updated.\n")
        newPDBFile = PDBFile("temp-pdb.pdb")
        outputstring = self.workstruct.structure.location + "/" + bnm + "-pickle.dat"
        newpickle = open(outputstring,"w") #we're still os.chdir'd to loc
        cPickle.dump(newPDBFile, newpickle)
        logfp.write("New pickle file created.\n")
        os.remove("temp-pdb.pdb")
        logfp.write("Oldtmp file removed.\n")
        logfp.write(outputstring + "\n")
        randomws = WorkingStructure.objects.get(id=self.workstruct.id)
        randomws.localpickle = outputstring #Make sure localpickle points to the right place...
        randomws.finalTopology = 'alfalfa'
#        self.workstruct.addCRDToPickle(wf.path,'mut_' + self.workstruct.identifier) #THis seems like a bad holdover from -minimization...
        logfp.write("LP1: " + randomws.localpickle + "\n")
        randomws.save() #Make sure localpickle is saved right...
        logfp.write("LP2: " + randomws.localpickle + "\n")
        logfp.close()
        self.status = 'C'
        self.localpickle = outputstring
        self.save()
