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
from structure.models import WorkingStructure, WorkingFile, Task

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
    def combinePDBsForMovie(self):
        ter = re.compile('TER') 
        remark = re.compile('REMARK')

        #movie_handle will be the final movie created
        #frame movie is the PDBs which each time step sepearated into a new PDB
        movie_handle = open(self.workstruct.structure.location + '/' + self.workstruct.identifier + "-" + self.action + '-movie.pdb','a')
        for i in range(1,11):
            frame_handle = open(self.workstruct.structure.location +  "/" + self.workstruct.identifier + "-" + self.action + "-movie" + str(i) + ".pdb",'r')
            movie_handle.write('MODEL ' + str(i) + "\n")
            for line in frame_handle:
                if not remark.search(line) and not ter.search(line): movie_handle.write(line)
            movie_handle.write('ENDMDL\n')
        movie_handle.close()

        self.movie_status = 'Done'
        self.save()
        return "Done."


    def finish(self):
        """test if the job suceeded, create entries for output"""

        logfp = open('/tmp/finishdyn.txt', 'w')

        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier

        # There's always an input file, so create a WorkingFile
        # for it.
        logfp.write('Create input WF\n')
        wfinp = WorkingFile()
        wfinp.task = self
        wfinp.path = loc + '/' + bnm + '-' + self.action + '.inp'
        wfinp.canonPath = wfinp.path
        wfinp.type = 'inp'
        wfinp.description = self.action.upper() + ' input script'
        wfinp.save()

        # Check if an output file was created and if so create
        # a WorkingFile for it.
        logfp.write('Check output...\n')
        try:
            os.stat(loc + '/' + bnm + '-' + self.action + '.out')
        except:
            logfp.write("Flunked on " + loc + '/' + bnm + '-' + self.action + '.out\n')
            self.status = 'F'
            self.save()
            return
        logfp.write('Done!\n')
        logfp.close()

        wfout = WorkingFile()
        wfout.task = self
        wfout.path = loc + '/' + bnm + '-' + self.action + '.out'
        wfout.canonPath = wfout.path
        wfout.type = 'out'
        wfout.description = 'output from ' + self.action.upper()
        wfout.save()

        # check if the final coordinates were created, if so then
        # we can also add the PSF, PDB, trajectory, and restart files
        try:
            os.stat(loc + '/' + bnm + '-' + self.action + '.crd')
        except:
            self.status = 'F'
            self.save()
            return

        wfcrd = WorkingFile()
        wfcrd.task = self
        wfcrd.path = loc + '/' + bnm + '-' + self.action + '.crd'
        wfcrd.canonPath = wfcrd.path
        wfcrd.type = 'crd'
        wfcrd.description = 'coordinates from ' + self.action.upper()
        wfcrd.save()
        self.workstruct.addCRDToPickle(wfcrd.path, self.action + '_' + self.workstruct.identifier)

        wfpsf = WorkingFile()
        wfpsf.task = self
        wfpsf.path = loc + '/' + bnm + '-' + self.action + '.psf'
        wfpsf.canonPath = wfpsf.path
        wfpsf.type = 'psf'
        wfpsf.description = 'PSF from ' + self.action.upper()
        wfpsf.save()

        wfpdb = WorkingFile()
        wfpdb.task = self
        wfpdb.path = loc + '/' + bnm + '-' + self.action + '.pdb'
        wfpdb.canonPath = wfpdb.path
        wfpdb.type = 'pdb'
        wfpdb.description = 'PDB coordinates from ' + self.action.upper()
        wfpdb.save()

        if self.make_movie:
            # go to hollywood
            self.combinePDBsForMovie()
            wfmoviePDB = WorkingFile()
            wfmoviePDB.path = loc + '/' + bnm + '-' + self.action + '-movie.pdb'
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
