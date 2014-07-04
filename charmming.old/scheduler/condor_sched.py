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
import os, time, string, commands, re
from tempfile import *

class jobScheduler:
   
   def getState(self,condor_id):
      #self.logfd.write("getState asked to look up %s\n" % condor_id)
      #self.logfd.flush()
      status, out = commands.getstatusoutput("%s %s" % (self.querycmd,condor_id))
      if status != 0:
         return -1
      for line in out.split('\n'):
         line = line.strip()
         #self.logfd.write("line> %s\n" % line)
         #self.logfd.flush()
         if line.startswith(condor_id):
            # job is queued or running, for now return running
            #self.logfd.write("OK found it\n")
            larr = string.splitfields(line)
            if larr[5] == 'R':
               #self.logfd.write("Running return 2\n")
               #self.logfd.flush()            
               return 2
            else:
               #self.logfd.write("Queued return 1\n")
               #self.logfd.flush()            
               return 1
      # if the job exists and is not running, it must be done
      #self.logfd.write("Did not find it, return 3\n")
      #self.logfd.flush()
      return 3

   def killJob(self,condor_id):
      status, out = commands.getstatusoutput("%s %s" % (self.killcmd,condor_id))
      if status != 0:
         return -1
      return 0

   def newJob(self,user,dir,exe,scripts):
      try:
         os.chdir(dir)
      except:
         self.logfd.write("Could not change to directory '%s'\n" % dir)
         return None
      self.logfd.write("Trying to start job %s %s %s %s\n" % (user,dir,exe,scripts))
      self.logfd.flush()
      
      if len(scripts) > 1:

         # this needs to be done in two steps -- first create the individual script files
         # and then create the DAG 

         # create individual script files
         subfiles = []
         for script in scripts:
            tfo = mkstemp(prefix="lcsub-dags-",dir="/tmp")
            tf = os.fdopen(tfo[0], 'w+b')
            tf.write("# Condor submit file for learn CHARMM\n")
            tf.write("# Created %s for charmming user id %d\n\n" % (time.ctime(), user))
            tf.write("Executable      = /usr/local/charmming/gfortran-xxlg.one\n")
            tf.write("Universe        = vanilla\n")
            tf.write("input           = %s\n" % script)
            tf.write("output          = %s\n" % script.replace("inp","out"))
            tf.write("error           = %s\n" % script.replace("inp","err"))
            tf.write("log             = %s\n" % script.replace("inp","log"))
            tf.write("Queue\n")

            tf.close()
            subfiles.append(tfo[1])

         # create the DAG file
         tfd = mkstemp(prefix="lcsub-dagm-",dir="/tmp")
         tf = os.fdopen(tfd[0], 'w+b')
         tf.write("# DAG job master for learn CHARMM\n\n")
         for i in range(len(scripts)):
            tf.write("Job STEP%d  %s\n" % (i+1,subfiles[i]))
         tf.write("\n")
         for j in range(len(scripts)-1):
            tf.write("PARENT STEP%d CHILD STEP%d\n" % (j+1,j+2))
         tf.close()
         cmd = "%s %s" % (self.dagcmd,tfd[1])

      else:
         script = scripts[0]
         tfo = mkstemp(prefix="lcsub-",dir="/tmp")
         tf = os.fdopen(tfo[0], 'w+b')
         tf.write("# Condor submit file for learn CHARMM\n")
         tf.write("# Created %s for charmming user id %d\n\n" % (time.ctime(), user))
         tf.write("Executable      = /usr/local/charmming/gfortran-xxlg.one\n")
         tf.write("Universe        = vanilla\n")
         tf.write("input           = %s\n" % script)
         tf.write("output          = %s\n" % script.replace("inp","out"))
         tf.write("error           = %s\n" % script.replace("inp","err"))
         tf.write("log             = %s\n" % script.replace("inp","log"))
         tf.write("Queue\n")
         tf.close()
         cmd = "%s %s" % (self.subcmd,tfo[1])

      self.logfd.write("Issuing %s\n" % cmd)
      self.logfd.flush()
      status, out = commands.getstatusoutput(cmd)
      self.logfd.write("RESULT:\n\n%s\n" % out)
      self.logfd.flush() 
      if status != 0:
         self.logfd.write("Bad status from command.\n")
         self.logfd.flush()
         return None
      if len(scripts) > 1:
         # DAG submission

         # hack to force rescheduling to avoid the 5 minute wait in DAG jobs...
         s2, o2 = commands.getstatusoutput("condor_reschedule -all")

         clmatch = re.compile("submitted to cluster ([0-9]+)")
         for line in out.split("\n"):
            x = clmatch.search(line)
            if x:
               condorID = x.group(1) + ".0"
               self.logfd.write("Got ID %s\n" % condorID)
               self.logfd.flush()

               # unlike PBS, we can't determine job IDs in advance so assign all scripts
               # the ID of the DAG manager.
               rlist = []
               for i in range(len(scripts)):
                  rlist.append(condorID)
               return rlist
      else:
         # regular condor submission
         for line in out.split("\n"):
            line = line.strip()
            if line.startswith("** Proc"):
               line = line.replace("** Proc ","")
               line = line.replace(":","")
               condorID = line
               self.logfd.write("Got ID %s\n" % condorID)
               self.logfd.flush()
               return [condorID]

      self.logfd.write("Did not find job number\n")
      self.logfd.flush()
      return None

   def __init__(self):
      # set up condor environment variables
      self.subcmd = "condor_submit -verbose"
      self.dagcmd = "condor_submit_dag"
      self.killcmd = "condor_rm"
      self.querycmd = "condor_q"
      os.environ['CONDOR_CONFIG'] = "/opt/condor-6.8.4/etc/condor_config"

      self.logfd = open("condor.log", "a+")
