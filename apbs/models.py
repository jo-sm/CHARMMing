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
from structure.models import Task, WorkingFile

import os


class apbsParams(models.Model):
    gdim_x    = models.PositiveIntegerField(default=0)
    gdim_y    = models.PositiveIntegerField(default=0)
    gdim_z    = models.PositiveIntegerField(default=0)
    npoint_x  = models.PositiveIntegerField(default=0)
    npoint_y  = models.PositiveIntegerField(default=0)
    npoint_z  = models.PositiveIntegerField(default=0)

class redoxTask(Task):
    redoxsite  = models.CharField(max_length=250) # segid + 1-/2- or 2-/3-
    apbsparams = models.ForeignKey(apbsParams,null=True)
    #couple = models.CharField(max_length=250)
    #chargeko = models.CharField(max_length=250)

    def finish(self):
        # This is a pain, since we need to get all of the input, output, and
        # PDB files from the REDOX calculation.
        """test if the job suceeded, create entries for output"""
        
        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier

        # There's always an input file, so create a WorkingFile
        # for it.
        filenamesuffix=['-build-redall',
                        '-build-redsite',
                        '-build-oxiall',
                        '-build-oxisite',
                        '-combinegrid-site',
                        '-mkgrid-full',
                        '-mkgrid-site',
                        '-modpot',
                        '-modpotref',
                        '-oxipot',
                        '-oxipotref']

        for suffix in filenamesuffix:
            wfinp = WorkingFile()
            wfinp.task = self
            wfinp.path = loc + '/redox-' + bnm + suffix + '.inp'
            wfinp.canonPath = wfinp.path
            wfinp.type = 'inp'
            wfinp.description = 'a redox input script'
            wfinp.save()

        # Check if an output file was created and if so create
        # a WorkingFile for it.
        logfp = open('/tmp/redoxdone.txt', 'w')
        for suffix in filenamesuffix:
            logfp.write('looking for %s\n' % (loc + '/redox-' + bnm + suffix + '.out'))
            try:
                os.stat(loc + '/redox-' + bnm + suffix + '.out')
                logfp.write('OK!\n')
            except:
                logfp.write('Bad\n')
                self.status = 'C'
                self.finished = 'y'
                return
            wfout = WorkingFile()
            wfout.task = self
            wfout.path = loc + '/redox-' + bnm + suffix + '.out'
            wfout.canonPath = wfout.path
            wfout.type = 'out'
            wfout.description = 'a redox script output'
            wfout.save()
        self.status = 'C'
        logfp.close()

        if self.status == 'F':
            return

        # create Working files for PDB, CRD, and PSF.
        '''
        OH WOW! I just stopped caring!!!!
        '''

        self.status = 'C'
        #self.finished = 'y'
