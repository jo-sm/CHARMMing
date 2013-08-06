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
from django.contrib.auth.models import User
import django.forms, os, re
import structure
from structure.models import WorkingStructure, WorkingFile, Task, CGWorkingStructure

#Replica exchange parameters
class rexParams(models.Model):

    firstStep = models.PositiveIntegerField(default=1)
    lastStep = models.PositiveIntegerField(default=2)
    cleanStep = models.PositiveIntegerField(default=50)
    npref = models.DecimalField(default=0,null=False,max_digits=8,decimal_places=3)
    nbath = models.PositiveIntegerField(default=4)
    temperatures = models.CharField(max_length=250)


class dynamicsTask(Task):

    nstep = models.PositiveIntegerField(default=1000)
    make_movie = models.BooleanField(default=False)
    movie_status = models.CharField(max_length=250,null=True)
    replica_exchange = models.ForeignKey(rexParams,null=True)
    scpism = models.BooleanField(default=False)

    #pre: Requires a Structure object
    #Combines all the smaller PDBs make in the above method into one large PDB that
    #jmol understands
    #type is md,ld,or sgld 
    def combinePDBsForMovie(self, inp_file, cgws):
        ter = re.compile('TER') 
        remark = re.compile('REMARK')

        #movie_handle will be the final movie created
        #frame movie is the PDBs which each time step sepearated into a new PDB
        movie_handle = open(self.workstruct.structure.location + '/' + self.workstruct.identifier + "-" + self.action + '-movie.pdb','a')
        for i in range(1,11):
            frame_fname = self.workstruct.structure.location + '/' + self.workstruct.identifier + "-" + self.action + "-movie" + str(i) + ".pdb"
            if inp_file != False:
                frame_status = cgws.addBondsToPDB(inp_file, frame_fname)
                if frame_status != True:
                    self.movie_status = 'Failed'
                    self.save()
                    return 'Failed'
            frame_handle = open(frame_fname,'r')
            movie_handle.write('MODEL ' + str(i) + "\n")
            for line in frame_handle:
                if not remark.search(line) and not ter.search(line): movie_handle.write(line)
            movie_handle.write('ENDMDL\n')
        movie_handle.close()

        self.movie_status = 'Done'
        self.save()
        return "Done"


    def finish(self):
        """test if the job suceeded, create entries for output"""

        logfp = open('/tmp/finishdyn.txt', 'w')

        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier

        # There's always an input file, so create a WorkingFile
        # for it.
        basepath= loc + "/" + bnm + "-" + self.action

        logfp.write('Create input WF\n')
        path = basepath + ".inp"
        wfinp = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfinp.task = self
            wfinp.path = path
            wfinp.canonPath = wfinp.path
            wfinp.type = 'inp'
            wfinp.description = self.action.upper() + ' input script'
            wfinp.save()

        # Check if an output file was created and if so create
        # a WorkingFile for it.
        logfp.write('Check output...\n')
        path = basepath + ".out"
        try:
            os.stat(path)
        except:
            logfp.write("Flunked on " + path + '\n')
            self.status = 'F'
            self.save()
            return
        logfp.write('Done!\n')
        logfp.close()

        wfout = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfout.task = self
            wfout.path = path
            wfout.canonPath = wfout.path
            wfout.type = 'out'
            wfout.description = 'output from ' + self.action.upper()
            wfout.save()

        # check if the final coordinates were created, if so then
        # we can also add the PSF, PDB, trajectory, and restart files
        path = basepath + ".crd"
        try:
            os.stat(path)
        except:
            self.status = 'F'
            self.save()
            return

        wfcrd = WorkingFile()
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfcrd.task = self
            wfcrd.path = path
            wfcrd.canonPath = wfcrd.path
            wfcrd.type = 'crd'
            wfcrd.description = 'coordinates from ' + self.action.upper()
            wfcrd.save()
            self.workstruct.addCRDToPickle(wfcrd.path, self.action + '_' + self.workstruct.identifier)

        wfpsf = WorkingFile()
        path = basepath + ".psf"
        inp_file = path
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfpsf.task = self
            wfpsf.path = path
            wfpsf.canonPath = wfpsf.path
            wfpsf.type = 'psf'
            wfpsf.description = 'PSF from ' + self.action.upper()
            wfpsf.save()

        wfpdb = WorkingFile()
        path = basepath + ".pdb"
        out_file = path
        try:
            wftest = WorkingFile.objects.get(task=self,path=path)
        except:
            wfpdb.task = self
            wfpdb.path = path
            wfpdb.canonPath = wfpdb.path
            wfpdb.type = 'pdb'
            wfpdb.description = 'PDB coordinates from ' + self.action.upper()
            wfpdb.save()
        cgws = False
        #Generic coarse-grain code goes here.
        logfp = open("/tmp/speed_test.txt","w")
        try:
            cgws = CGWorkingStructure.objects.get(workingstructure_ptr=self.workstruct.id)
            logfp.write("Found a CGWorkingStructure.")
        except CGWorkingStructure.MultipleObjectsReturned: #Uh oh. This MAY be alright if AA/CG...
            self.status = "F"
            logfp.write("Found too many CGWorkingStructures.")
            return
        except Exception as ex: #Catch everything else
            logfp.write(str(ex))
            logfp.flush()
#            pass
        if cgws != False:
            cgws.addBondsToPDB(inp_file,out_file)
        logfp.close()
        if self.make_movie:
            # go to hollywood
            if cgws:
                result = self.combinePDBsForMovie(inp_file, cgws)
            else:
                result = self.combinePDBsForMovie(False, False)
            if result != "Done":
                self.status = "F"
                return
            wfmoviePDB = WorkingFile()
            path = basepath + "-movie.pdb"
            try:
                wftest = WorkingFile.objects.get(task=self,path=path)
            except:
                wfmoviePDB.task = self
                wfmoviePDB.path = path
                wfmoviePDB.canonPath = wfmoviePDB.path
                wfmoviePDB.type = 'pdb'
                wfmoviePDB.description = 'Movie PDB for ' + self.action.upper()
                wfmoviePDB.save()

class mdTask(dynamicsTask):

    ensemble = models.CharField(max_length=50)

    #temp will represent the temperature Kelvin if "type" is heat
    temp = models.FloatField(null=True)
    firstt = models.FloatField(null=True)
    finalt = models.FloatField(null=True)
    teminc = models.FloatField(null=True)
    ihtfrq = models.FloatField(null=True)
    tbath = models.FloatField(null=True)


class ldTask(dynamicsTask):

    fbeta = models.FloatField(default=60.0,null=True)
    sgld = models.BooleanField(default=False)


class sgldTask(ldTask):

    tsgavg = models.FloatField(default=0.5,null=True)
    tempsg = models.FloatField(default=1.0,null=True)
    make_sgld_movie = models.BooleanField(default=False)
    sgld_movie_status = models.CharField(max_length=250,null=True)
    sgld_movie_req = models.CharField(max_length=250,null=True)
    sgld_movie_req = models.BooleanField(default=False)

# Create your models here.
