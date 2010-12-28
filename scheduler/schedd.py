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

import sys, os, time, string, select, getopt, signal, MySQLdb
import commands
import errno
from socket import *
from pbs_sched import jobScheduler

## Global variables
port = 9995

## classes
class logger:
   loglevel = 4

   def write_log(self, msglvl, msg):
      if msglvl <= self.loglevel:
         self.__fd.write("%s <%d>: %s\n" % (time.ctime(), msglvl, msg))
         self.__fd.flush()
      if msglvl <= self.bomlev:
         self.__fd.write("Fatal error, shutting down %d.\n" % os.getpid())
         self.__fd.close()
         sys.exit(-1)

   def shutDown(self):
      self.write_log(2, "Closing log.")
      self.__fd.close()

   def __init__(self, fname, bomlev):
      try:
         self.__fd = open(fname, 'a+')
         self.bomlev = bomlev
      except:
         sys.stderr.write("Logger class could not open %s.\n" % fname)
         return None

class tracker:
   shutdown = 0

   def shutDown(self):
      self.logif.write_log(2, "Job tracker shutdown complete.")

   def updateStatus(self):
      # initialize MySQL
      try:
         dba = MySQLdb.connect("localhost", user="charmming", passwd="charmming", db="charmming")
      except:
         self.logif.write_log(0, "Authentication to MySQL db failed!")
         return None

      cra = dba.cursor()

      sqlstat = "SELECT id,sched_id,state,script,batchid FROM job_scheduler WHERE ended is null"
      try:
         cra.execute(sqlstat)
      except:
         self.logif.write_log(4, "Running '%s' failed!" % sqlstat)
         cra.close()
         dba.close()
         return None

      rset = cra.fetchall()
      for row in rset:
         try:
            j_id = int(row[0])
            j_schedid = row[1]
            j_state = int(row[2])
            j_script = row[3]
            j_batchid = row[4]
         except:
            self.logif.write_log(1, "Unconvertable value from db id %s state %s" % (row[0],row[2]))
            return

         currstate = self.js.getState(j_schedid)

         if currstate < 0:
            self.logif.write_log(1, "Error fetching current state information for job %d." % j_id)
            continue 

         #self.logif.write_log(4, "updateStatus ... checking currstate.")
         if currstate != j_state:

            if currstate == 3:
               # we need to figure out if the job is *really* done or failed
               j_outfile = j_script.replace(".inp",".out")
               self.logif.write_log(4, "Checking for normal termination in %s." % j_outfile)
               cmd = "grep -iw 'normal termination' %s" % j_outfile
               self.logif.write_log(4, "updateStatus ... running %s." % cmd)
               print cmd
               s, o = commands.getstatusoutput(cmd)
               if s != 0:
                  # couldn't find a normal termination line, have to assume that the job failed...
                  self.logif.write_log(2, "Uh oh ... looks like job %d failed -- grep returned %s." % (j_id,s))
                  currstate = 4
               else:
                  # looks OK, check all jobs to make sure no prior failures.
                  crx = dba.cursor()
                  sqlstat = "SELECT script FROM job_scheduler WHERE batchid='%s'" % j_batchid
                  crx.execute(sqlstat)
                  rset = crx.fetchall()
                  for row in rset:
                      j2_outfile = row[0].replace(".inp", ".out")
                      cmd = "grep -iw 'normal termination' %s" % j_outfile
                      self.logif.write_log(4, "updateStatus ... running %s." % cmd)
                      print cmd
                      s, o = commands.getstatusoutput(cmd)
                      if s != 0:
                         # couldn't find a normal termination line, have to assume that the job failed...
                         self.logif.write_log(2, "Uh oh ... looks like job %d failed -- grep returned %s." % (j_id,s))
                         currstate = 4
                         break
                  crx.close()

            ltime = time.localtime()
            self.logif.write_log(3, "Job %d went from state %d to %d." % (j_id,j_state,currstate))
            sstat2 = "UPDATE job_scheduler SET state=%d WHERE id=%d" % (currstate,j_id)
            crb = dba.cursor()
            crb.execute(sstat2)

            if currstate >= 2:
               if currstate == 2:
                  sstat2 = "UPDATE job_scheduler SET started='%s' WHERE id=%d" % (MySQLdb.Timestamp(ltime[0],ltime[1],ltime[2],ltime[3],ltime[4],ltime[5]),j_id)
               elif currstate == 3 or currstate == 4:            
                  sstat2 = "UPDATE job_scheduler SET ended='%s' WHERE id=%d" % (MySQLdb.Timestamp(ltime[0],ltime[1],ltime[2],ltime[3],ltime[4],ltime[5]),j_id)
               try:
                  crb.execute(sstat2)
                  crb.close()
               except:
                  self.logif.write_log(1, "Failed updating job time on '%s'" % sqlstat)


      cra.close()
      dba.close()

   def lookupJobStatus(self,jid):
      # initialize MySQL
      try:
         dbc = MySQLdb.connect("localhost", user="charmming", passwd="charmming", db="charmming")
         crc = dbc.cursor()
      except:
         self.logif.write_log(0, "Authentication to MySQL db failed!")
         return None

      # query DB
      sqlstat = "SELECT userid,state,started,ended,queued FROM job_scheduler WHERE id=%d" % jid
      try:
         crc.execute(sqlstat)
      except:
         self.logif.write_log(1, "Error executing SQL '%s'" % sqlstat)
         return None

      rarr = crc.fetchall()
      crc.close()
      dbc.close()

      if len(rarr) != 1:
         return None
      if rarr[0][1] == 0:
         sstr = 'submitted'
      elif rarr[0][1] == 1:
         sstr = "queued since %s" % rarr[0][4]
      elif rarr[0][1] == 2:
         sstr = "running since %s" % rarr[0][2]
      elif rarr[0][1] == 3:
         sstr = "complete since %s" % rarr[0][3]
      elif rarr[0][1] == 4:
         sstr = 'failed'
      else:
         self.logif.write_log(1, "Bad state received for job %d." % jid)
         return None
      rstr = "Job %d User %s %s" % (jid,rarr[0][0],sstr)
      self.logif.write_log(3, "Returning '%s' to CHARMMING." % rstr)
      return rstr

   def killJob(self,request):
      request = request.strip()
      try:
         self.logif.write_log(4, "killJob request: %s" % request)
         rarr = string.splitfields(request)
         user = int(rarr[0])
         jid = rarr[1]
      except:
         self.logif.write_log(1, "Malformed delete request string '%s'" % request)
         return -1

      # initialize MySQL
      try:
         dbc = MySQLdb.connect("localhost", user="charmming", passwd="charmming", db="charmming")
         crc = dbc.cursor()
      except:
         self.logif.write_log(0, "Authentication to MySQL db failed!")
         return -1

      # get the batch ID and kill all jobs in the batch
      sqlstat = "SELECT batchid FROM job_scheduler WHERE id=%s" % jid
      try:
         crc.execute(sqlstat)
      except:
         self.logif.write_log(1, "Error executing SQL '%s'" % sqlstat)
         return -1
      rarr = crc.fetchall()
      if len(rarr) < 1:
        self.logif.write_log(2, "Could not find any jobs with ID %s" % jid)
        return -3
      batchid = rarr[0][0]

      sqlstat = "SELECT userid,state,id,sched_id FROM job_scheduler WHERE batchid=%s" % batchid
      try:
         crc.execute(sqlstat)
      except:
         self.logif.write_log(1, "Error executing SQL '%s'" % sqlstat)
         return -2

      rarr = crc.fetchall()
      if len(rarr) < 1:
         self.logif.write_log(2, "Could not find any jobs with batch ID %s" % jid)
         return -3
      for rline in rarr:
         if int(rline[0]) != user:
             self.logif.write_log(1, "User IDs do not match on job delete request, deletion attempt by %d on job belonging to %d " % (int(rarr[0][0]),user))
             return -3
         # actually kill the job, killed jobs are marked as failed.
         sqlstat = "UPDATE job_scheduler SET state=4 where id=%d" % int(rline[2])
         try:
            crc.execute(sqlstat)
         except MySQLdb.Error, e:
            self.logif.write_log(1, "Error executing SQL '%s'" % sqlstat)
            self.logif.write_log(1, "Error %s" % (e.args))
            return -1
         except:
            self.logif.writelog(1, "Unknown error on SQL statement '%s'" % sqlstat)
            return -1

         self.js.killJob(rline[3])

      crc.close()
      dbc.close()
      return 0
        
   def startJob(self,request,execs,nprocs):
      request = request.strip()
      execs = execs.strip()
      nprocs = nprocs.strip()
      try:
         self.logif.write_log(4, "startJob request is %s." % request)
         self.logif.write_log(4, "startJob execs is %s." % execs)
         self.logif.write_log(4, "startJob nprocs is %s." % nprocs)
         rarr = string.splitfields(request)
         user = int(rarr[0])
         dir  = rarr[1]
         scripts = rarr[2:]
      except:
         self.logif.write_log(1, "Malformed job request string '%s'" % request)
         return None

      exelist = execs.split()
      nprlist = nprocs.split()
      if len(exelist) != len(scripts) or len(nprlist) != len(scripts):
         self.logif.write_log(1, "Mismatched list lengths! Lengths are %d (scripts) %d (execs) %d (nproc)." % (len(scripts),len(exelist),len(nprlist)))
         return None

      # initialize MySQL
      try:
         dbc = MySQLdb.connect("localhost", user="charmming", passwd="charmming", db="charmming")
         crc = dbc.cursor()
      except:
         self.logif.write_log(0, "Authentication to MySQL db failed!")
         return None

      # create the job with the real scheduling module to get the sched_id
      self.logif.write_log(4, "startJob: creating job u %d d %s s %s" % (user,dir,scripts))
      schedids = self.js.newJob(user,dir,scripts,exelist,[ int(x) for x in nprlist ])
      if not schedids:
         self.logif.write_log(2, "startJob: job failure user %d script %s" % (user,scripts))
         crc.close()
         dbc.close()
         return None
      batchid = schedids[-1]
      ltime = time.localtime()
      
      for i in range(len(scripts)):
         sqlstat = "INSERT INTO job_scheduler (userid,sched_id,batchid,state,exe,script,queued) VALUES "
         sqlstat += "(%d,'%s','%s',0,'%s','%s','%s')" % (user,schedids[i],batchid,exelist[i],scripts[i],MySQLdb.Timestamp(ltime[0],ltime[1],ltime[2],ltime[3],ltime[4],ltime[5]))
         crc.execute(sqlstat)
         self.logif.write_log(4, "startJob: created %s OK" % scripts[i])

      # we only return the ID of the last (significant) job in the batch.
      self.logif.write_log(4, "startJob: sched_id %s and script %s" % (schedids[-1],scripts[-1]))
      sqlstat = "SELECT id FROM job_scheduler WHERE sched_id='%s' AND script='%s'" % (schedids[-1],scripts[-1])
      crc.execute(sqlstat)
      rarr = crc.fetchall()
      jobid = rarr[0][0]

      # update status
      for i in range(len(scripts)):
         tvar = self.js.getState(schedids[i])
         if tvar != 0:
            sqlstat = "UPDATE job_scheduler set state=%d where sched_id='%s' AND script='%s'" % (tvar,schedids[i],scripts[i])
            crc.execute(sqlstat)
         self.logif.write_log(4, "Exiting startJob with job number %d for script %s and state at %d." % (jobid,scripts[i],tvar))

      # return ID to user
      crc.close()
      dbc.close()
      return "queued %d" % jobid

   def handle_conn(self,sock):
      line = sock.recv(256)
      line = line.strip()
      if line != "HELLO SCHEDD":
         sock.send("BAD PROTOCOL\n")
         return
      sock.send("HELLO CLIENT\n")
      while 1:
         request = sock.recv(256)
         request = request.strip()
         if request.startswith("STATUS"):
            sock.send("OK WHATNUM\n")
            line = sock.recv(256)
            line = line.strip()
            try:
               jobid = int(line)
            except:
               sock.send("BAD JOB ID\n")
               continue
            s = self.lookupJobStatus(jobid)
            if not s:
               self.logif.write_log(3, "STATUS ret bad job!")
               sock.send("BAD NO SUCH JOB\n")
            else:
               self.logif.write_log(3, "STATUS ret %s" % s)
               sock.send("OK %s\n" % s)
         elif request == "NEW":
            sock.send("OK GIVE USER DIR SCRIPTS\n")
            udsline = sock.recv(1024)
            udsline = udsline.strip()
            sock.send("OK GIVE EXELIST\n")
            exeline = sock.recv(1024)
            exeline = exeline.strip()
            sock.send("OK GIVE NPROCLIST\n")
            nprline = sock.recv(1024)
            nprline = nprline.strip()

            s = self.startJob(udsline,exeline,nprline)
            if not s:
               sock.send("BAD JOB SPEC\n")
            else:
               sock.send("OK %s\n" % s)
         elif request == "ABORT":
            sock.send("OK GIVE USER AND NUMBER\n")
            line = sock.recv(1024)
            line = line.strip()
            s = self.killJob(line)
            if s != 0:
               sock.send("BAD JOB: ERR %s\n" % s)
            else:
               sock.send("OK DONE\n")
         elif request == "END":
            sock.send("OK GOODBYE\n")
            return
         else:
            try:
               sock.send("BAD REQUEST\n")
            except:
               pass
            return

   def __init__(self, log):
      self.logif = log
      self.js = jobScheduler()
 

## functions

# variables global to functions but not classes
logif = None
jobTrack = None

def become_daemon():
   if os.fork() == 0:
      os.setsid()
      sys.stdout = open('/dev/null', 'w')
      sys.stdin = open('/dev/null', 'r')
      if os.fork() == 0:
         return

   sys.exit(0)

def shutdown_fn():
   jobTrack.shutDown()
   logif.shutDown()
   sys.exit(0)

def handle_term(signo, frame):
   logif.write_log(2, "Got shut down signal.")
   jobTrack.shutdown = 1

def handle_hup(signo, frame):
   logif.write_log(2, "Got SIGHUP ... reopening log.")
   logif.shutDown()
   logif = logger("schedd.log", 0)

## main program
if __name__ == "__main__":
   # process command line arguments
   host = '0.0.0.0'
   opts, args = getopt.getopt(sys.argv[1:], 'h:p:') 
   for o, a in opts:
      if o == '-h':
         host = a
      elif o == '-p':
         try:
            port = int(a)
         except:
            print "%s is not a valid port." % a
            sys.exit(-1)
      else:
         print "Unknown argument %s." % o
         sys.exit(-1)

   # set up handlers for SIGINT, SIGTERM and (todo) SIGHUP
   signal.signal(signal.SIGINT, handle_term)
   signal.signal(signal.SIGTERM, handle_term)
   signal.signal(signal.SIGHUP, handle_hup)

   # set up interfaces and go into main loop
   logif = logger("schedd.log", 0)
   if not logif:
      sys.stderr.write("Failed initializing logger.\n")
      sys.exit(-1)
   jobTrack = tracker(logif)
   if not jobTrack:
      sys.stderr.write("Could not start job tracker.\n")
      sys.exit(-1)
   logif.write_log(3,"Tracker initialized in main.")

   try:
      s = socket(AF_INET, SOCK_STREAM)
      s.bind(("0.0.0.0", port))
      s.listen(25)
   except:
      logif.write_log(0, "Could not set up the socket. Check nothing else is bound to port %d." % port)

   become_daemon()
   logif.write_log(2, "Daemonization complete, running PID is %d." % os.getpid())
   while jobTrack.shutdown == 0:
      isr = [s]
      isw = []
      isx = []

      try:
         r, w, e = select.select(isr, isw, isx, 1.0)
      except select.error, er:
         if er[0] != errno.EINTR:
            # interrupted system calls are harmless
            logif.write_log(1, "Got errno %d in select!" % er[0])

      if r:
         c_sock, c_addr = s.accept()

         frkres = os.fork()
         if frkres:
            # parent process
            c_sock.close()
         else:
           # child process -- close listening socket and handle connection
           jt2 = tracker(logif)
           if not jt2:
              sys.stderr.write("Could not start job tracker.\n")
              sys.exit(-1)
           jt2.handle_conn(c_sock)
           c_sock.close()
           sys.exit(0)
      
      # clean up any dead children that we've acquired
      while 1:
         try:
            deadkid, deadstat = os.waitpid(-1, os.WNOHANG)
            if deadkid == 0:
               break
         except:
            break

      # update status of jobs
      jobTrack.updateStatus()

   # loop breaks -- need to shut down
   logif.write_log(2, "Shutdown request received.")
   shutdown_fn()
   s.close()
