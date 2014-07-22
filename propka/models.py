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
from pychm.io import pdb
from structure.models import Task, WorkingFile, CGWorkingStructure
import os, copy, re
import traceback
from lib import Atom

class propkaTask(Task):
    #doesn't actually have any fields, but we include it so we can get
    # a finish method and classification

    def finish(self):
        """test if the job suceeded, create entries for output"""

        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier
        basepath = loc + '/' + bnm + "-" + self.action
        #first see if we could write propka_input
        try:
            os.stat(basepath + ".propka_input")
        except:
            self.status = 'F'
            self.save()
        #might happen that the .pka isn't present. Unlikely, but possible.
        try:
            os.stat(basepath + '.pka')
        except:
            self.status = 'F'
            self.save()
            return
        if self.status == 'F':
            return
        self.status = 'C'
        self.createStatistics() #This only fails when the task succeeds, so we shouldn't worry.
        return
