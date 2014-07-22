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

import os, time, string, commands
import charmming_config
from tempfile import *

class jobScheduler:
   
   def getState(self,pbs_id):
      if self.debugMode:
         self.logfd.write("getState asked to look up %s\n" % pbs_id)
         self.logfd.flush()
      status, out = commands.getstatusoutput("%s %s" % (self.querycmd,pbs_id))
      if status != 0:
         # damn thing must not exist
         if self.debugMode:
            self.logfd.write("Did not find it -- return 3.\n")
         return 3
      for line in out.split('\n'):
         line = line.strip()
         if self.debugMode:
            self.logfd.write("line> %s\n" % line)
            self.logfd.flush()
         if line.startswith(pbs_id + ".ctb4"):
            self.logfd.write("OK found it\n")
            larr = string.splitfields(line)
            if larr[4] == 'R' or larr[4] == 'E' or larr[4] == 'H':
               self.logfd.write("Running return 2\n")
               self.logfd.flush()            
               return 2
            elif larr[4] == 'Q':
               self.logfd.write("Queued return 1\n")
               self.logfd.flush()            
               return 1
      # if the job exists and is not running, it must be done
      if self.debugMode:
         self.logfd.write("Did not find it, return 3.\n")
         self.logfd.flush()
      return 3

   def killJob(self,pbs_id):
      # returns 0 on success or an error code:
      # -1: error
      if self.debugMode:
         self.logfd.write("killJob asked to kill %s\n" % pbs_id)
         self.logfd.flush()
      cmd = "%s %s" % (self.killcmd,pbs_id)
      status, out = commands.getstatusoutput(cmd)
      if status != 0:
         if self.debugMode:
            self.logfd.write("Failed to kill!")
            self.logfd.flush()
         return -1

      self.logfd.write("Success!")
      self.logfd.flush()
      return 0

   def newJob(self,user,dir,scripts,exelist,nprlist):
      try:
         os.chdir(dir)
      except:
         return None
      if self.debugMode:
         self.logfd.write("newJob: Trying to start job %s %s %s\n" % (user,dir,scripts))
         self.logfd.flush()
      # create a temp file and put the PBS junk in it
      ids = []
      for i in range(len(scripts)):
         self.logfd.write("newJob: Trying script %s.\n" % (scripts[i]))
         self.logfd.write("newJob exelist: " + str(exelist) + "\n")
         self.logfd.flush()
         if not scripts[i].startswith('customscript'):
            self.logfd.write("Is NOT custom script!\n")
            self.logfd.flush()
            tfo = mkstemp(prefix="charmming-",dir="/tmp")
            tf = os.fdopen(tfo[0], 'w+b')
            tf.write("#!/bin/bash\n")
            tf.write("#PBS -l nodes=1:ppn=%d\n" % int(nprlist[i]))
            tf.write("#PBS -N charmming-job\n")
            tf.write("#PBS -j oe\n")
            if len(ids) > 0:
               tf.write("#PBS -W depend=afterok:%s\n" % ids[-1])
            tf.write("# PBS job file for CHARMMING\n")
            tf.write("# Created %s for charmming user id %d\n" % (time.ctime(), user))
            tf.write("# Number of processes requested: %d\n\n" % nprlist[i])
            tf.write("export GFORTRAN_UNBUFFERED_ALL=Y\n")
            tf.write("cd %s\n" % dir)
            tf.write("umask 0022\n")
#            if nprlist[i] == 1:
    #we can't rely on the number of processors being reliable. We make this check BEFOREHAND and use the right exec before passing it into here.
    #TODO: Update logic further.
            if exelist[i] == charmming_config.apps['charmm'] or exelist[i] == charmming_config.apps['charmm-mscale'] or exelist[i] == charmming_config.apps['charmm-apbs']:
                tf.write("%s < %s >& %s\n" % (exelist[i], scripts[i], scripts[i].replace("inp","out")))
            elif exelist[i] == charmming_config.apps['qchem']:
                tf.write("%s %s > %s\n" % (exelist[i], scripts[i], scripts[i].replace("inp","out")))
            elif exelist[i] == charmming_config.apps['propka']:
                tf.write("%s %s" %(exelist[i],scripts[i])) #in this case the "script' is really just a PDB file.
            else:
                tf.write("%s < %s >& %s\n" % (exelist[i], scripts[i], scripts[i].replace("inp","out")))
#            elif nprlist[i] > 1:
#               tf.write("%s %s < %s >& %s\n" % (charmming_config.mpirun_exe, exelist[i], scripts[i], scripts[i].replace("inp","out")))
            tf.close()
            cmd = "%s %s" % (self.subcmd,tfo[1])
         else:
            # the executable is something we can qsub directly
            self.logfd.write("Is INDEED a custom script!\n")
            self.logfd.flush()
            cmd = "%s %s" % (self.subcmd,exelist[i])

         if self.debugMode:
            self.logfd.write("newJob: Issuing %s\n" % cmd)
            self.logfd.flush()
         status, out = commands.getstatusoutput(cmd)
         if status != 0:
            if self.debugMode:
               self.logfd.write("Bad status from command.\n")
               self.logfd.flush()
            return None
         for line in out.split("\n"):
            try:
               line = line.strip()
               pbsID = line.split('.')[0]
               ids.append(pbsID)
            except:
               if self.debugMode:
                  self.logfd.write("Bad result from job submission.\n")
                  self.logfd.flush()
               return None

      return ids

   def __init__(self):
      self.subcmd = "qsub"
      self.querycmd = "qstat"
      self.killcmd = "qdel"
      self.debugMode = True

      if self.debugMode:
         self.logfd = open("pbs.log", "a+")
