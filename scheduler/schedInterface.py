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
from django.http import HttpResponse
from socket import *
import string
import charmming_config

class schedInterface:

   def checkStatus(self,jobID):
      self.sock.send("STATUS\n")
      buf = self.sock.recv(512)
      if not buf.startswith("OK"):
         return -1
      self.sock.send("%d\n" % jobID)
      buf = self.sock.recv(512)
      if not buf.startswith("OK"):
         return -2
      return buf[3:]  

   def submitJob(self,user,dir,scriptlist,exedict={},nprocdict={}):

      # error out if the length of the executable list or the
      # length of the nproc list is different from the length of
      # the scriptlist (assuming they exist).
      le = len(exedict.keys())
      lp = len(nprocdict.keys())
      if le > len(scriptlist):
         return -3
      if lp > len(scriptlist):
         return -3

      # build exelist and nproclist from exedict and nprocdict
      exelist = []
      nproclist = []
      for i in range(len(scriptlist)):
         if scriptlist[i] in nprocdict.keys():
            nproclist.append(nprocdict[scriptlist[i]])
         else:
            nproclist.append(1)

         if scriptlist[i] in exedict.keys():
            exelist.append(exedict[scriptlist[i]])
         elif nproclist[-1] > 1:
            exelist.append(charmming_config.charmm_mpi_exe)
         else:
            exelist.append(charmming_config.charmm_exe)

      self.sock.send("NEW\n")
      buf = self.sock.recv(512)
      if not buf.startswith("OK"):
         return -1
      scriptstr = string.join(scriptlist)
      self.sock.send("%d %s %s\n" % (user,dir,scriptstr))
      buf = self.sock.recv(512)
      if not buf.startswith("OK"):
         return -2

      nprbuf = ""
      npl2 = []
      for i in range(len(nproclist)):
         nprbuf += "%s " % nproclist[i]

      exebuf = ""
      for i in range(len(scriptlist)):
         if exelist:
            exebuf += "%s " % exelist[i]
         else:
            if npl2[i] == 1:
               exebuf += "%s " % charmming_config.charmm_exe
            elif npl2[i] > 1:
               exebuf += "%s " % charmming_config.charmm_mpi_exe

      self.sock.send("%s\n" % exebuf)
      buf = self.sock.recv(512)
      if not buf.startswith("OK"):
         return -2

      self.sock.send("%s\n" % nprbuf)
      buf = self.sock.recv(512)
      if not buf.startswith("OK"):
         return -2

      return int(buf.split()[2])

   def doJobKill(self,user,jobid):
      self.sock.send("ABORT\n")
      buf = self.sock.recv(512)
      if not buf.startswith("OK"):
         return -1
      self.sock.send("%d %s\n" % (user,jobid))
      buf = self.sock.recv(512)
      if not buf.startswith("OK"):
         return -2
      return 0

   def submitOrderedList(self,user,dir,nscript,slist):
      pass

   def __del__(self):
      self.sock.send("END\n")
      self.sock.close()

   def __init__(self):
      self.sock = socket(AF_INET,SOCK_STREAM)
      self.sock.connect(('localhost',9995))
      self.sock.send("HELLO SCHEDD\n")
      buf = self.sock.recv(1024)
      if buf != "HELLO CLIENT\n":
         self.sock.close()
         raise "Protocol error from server."

# kills a job based on user request...
def killJob(request,jobid):
    if not request.user.is_authenticated():
       	return render_to_response('html/loggedout.html')
    uid = request.user.id
    jobid = jobid.replace('/','')

    si = schedInterface()
    r = si.doJobKill(uid,jobid)
    return HttpResponse("code %s" % r)
