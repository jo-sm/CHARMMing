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
from scheduler.schedInterface import schedInterface
from structure.models import WorkingFile, Task, saveQCWorkingFiles
import os,re

class nmodeTask(Task):
    # type 1 = all atom, 2 - ENM
    type = models.PositiveIntegerField(default=0)

    nmodes = models.PositiveIntegerField(default=0)
    rcut = models.FloatField(default=0.0)
    kshort = models.FloatField(default=0.0)
    klong = models.FloatField(default=0.0)
    selected = models.CharField(max_length=1)
    nma_movie_status = models.CharField(max_length=250,default=None,null=True)
    make_nma_movie = models.BooleanField(default=False)
    nma_movie_req = models.BooleanField(default=False)
    num_trjs = models.PositiveIntegerField(default=5) #Holds how many trajectories are generated...
    useqmmm= models.CharField(max_length=1,null=True,default="n")
    modelType = models.CharField(max_length=30,null=True,default=None)

    #Combines all the smaller PDBs make in the above method into one large PDB that
    #jmol understands
    #type is nmode 
    def combinePDBsForMovie(self):
        ter = re.compile('TER')
        remark = re.compile('REMARK')
        #movie_handle will be the final movie created
        #frame movie is the PDBs which each time step sepearated into a new PDB
        for curr_traj in range(1,self.num_trjs+1):
            movie_handle = open(self.workstruct.structure.location + '/' + self.workstruct.identifier + "-" + self.action + '-' + str(curr_traj) + '-movie.pdb','a')
            for curr_frame in range(1,13):
                frame_handle = open(self.workstruct.structure.location +  "/" + self.workstruct.identifier + "-" + self.action + "-movie" + str(curr_frame) + '-' + str(curr_traj) + ".pdb",'r')
                movie_handle.write('MODEL ' + str(curr_frame) + "\n")
                for line in frame_handle:
                    if not remark.search(line) and not ter.search(line): movie_handle.write(line)
                movie_handle.write('ENDMDL\n')
                frame_handle.close()
            movie_handle.close()

        self.nma_movie_status = 'Done'
        self.save()
        return "Done."

    def finish(self):
        """test if the job suceeded, create entries for output"""
        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier

        # There's always an input file, so create a WorkingFile
        # for it.
        basepath = loc + '/' + bnm + '-' + self.action

        path = basepath + ".inp"
        wfinp = WorkingFile()
        try:
            WorkingFile.objects.get(task=self,path=path)
        except:
            wfinp.task = self
            wfinp.path = path
            wfinp.canonPath = wfinp.path
            wfinp.type = 'inp'
            wfinp.description = 'normal mode input script'
            wfinp.save()


        # Check if an output file was created and if so create
        # a WorkingFile for it.
        try:
            os.stat(loc + '/' + bnm + '-nmode.out')
        except:
            self.status = 'F'
            return

        path = basepath + ".out"
        wfout = WorkingFile()
        try:
            WorkingFile.objects.get(task=self,path=path)
        except:
            wfout.task = self
            wfout.path = path
            wfout.canonPath = wfout.path
            wfout.type = 'out'
            wfout.description = 'normal modes script output'
            wfout.save()

        if self.status == 'F':
            return

        #Generic Qchem inp/out function
        saveQCWorkingFiles(self,basepath)

        # if there is a movie, make the files.
        if self.make_nma_movie:
            self.combinePDBsForMovie()
            try:
                os.stat(basepath + "-1-movie.pdb")
            except:
                self.status = "F"
                return
            for curr_traj in range(1,self.num_trjs + 1):
                wfmovie = WorkingFile()
                path = basepath + "-" + str(curr_traj) + "-movie.pdb"
                try:
                    wftest = WorkingFile.objects.get(task=self,path=path)
                except:
                    wfmovie.task = self
                    wfmovie.path = path
                    wfmovie.canonPath = wfmovie.path
                    wfmovie.type = "pdb"
                    wfmovie.description = 'normal modes movie for trajectory ' + str(curr_traj)
                    wfmovie.save()
            for curr_traj in range(1,self.num_trjs+1):
                for curr_frame in range(1,13):
                    os.unlink(basepath + "-movie" + str(curr_frame) + "-" + str(curr_traj) + ".pdb")

        self.status = 'C'
        self.save()
