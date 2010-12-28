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
#!/usr/bin/python
# charmm22_format: puts the PDB into a format that can safely be read in CHARMM, returns the following error codes:
# 0: success
# -1: I/O error
# -2: parse error
# -3: no atoms found
# Tim Miller, 2007.

import os
import pdbinfo
from pdbinfo.models import PDBFile,ParseException
import lesson1

# split a multi model PDB into its constituent parts and creates a PDBFile model for each
def getModels(filename,owner,location,lessontype,lessonid):
   stripped_name = filename.split('.pdb')[0].lower() 
   infp = open(location + filename, 'r')
   outfp = None
   if not infp:
      return None
   modnum = 0
   natom = 0
   curr_fname = ""
   headertxt = ""
   rlist = []
   errorList = []

   for line in infp:
      if line.startswith('HEADER') or line.startswith('TITLE') or line.startswith('COMPND') or \
           line.startswith('SOURCE') or line.startswith('JRNL') or line.startswith('AUTHOR'):
         # ignore header information if already in a model
         if natom == 0:
            headertxt += "%s" % line
      if line.startswith('SSBOND'):
         # ignore header information if already in a model
         if natom == 0:
            headertxt += line
      if line.startswith('ATOM') or line.startswith('HETATM'):
         if modnum == 0:
            # not in a model, we must create one (only one model in file)
            modnum = 1
            natom = 1
            curr_fname = stripped_name + "-1.pdb"
            try:
               outfp = open("%s/%s" % (location,curr_fname), 'w')
               outfp.write(headertxt)
               outfp.write(line)
            except:
               return None
         else:
            natom += 1
            try:
               outfp.write(line)
            except:
               return None
      if line.startswith("MODEL"):
         if natom > 0:
            # uh-oh, sttarting a new model with atoms from the last one still
            # floating about -- better barf up a parse error.
            return None
         modnum += 1
         curr_fname = stripped_name + "-" + str(modnum) + ".pdb"
         try:
            outfp = open(location + '/' + curr_fname, 'w+')
            outfp.write(headertxt)
         except:
            return None
      if line.startswith("ENDMDL"):
         natom = 0
         outfp.write("TER\nEND\n")
         outfp.close()

         # OK, created a file for the model, now let's create the requisite
         # object.
         nfile = PDBFile()
         nfile.owner = owner
         nfile.location = location
         nfile.filename = curr_fname
         if lessontype == 'lesson1':
             lesson = lesson1.models.Lesson1.objects.filter(id = lessonid)[0]
         else:
             lesson = None
         if lesson:
             nfile.lesson_type = lessontype
             nfile.lesson_id = lessonid
         nfile.save()
	 #Add an error to the errorList if the parser crashes
         try:
            pdbinfo.views.parser(nfile)
            errorList.append('')
	    rval = 1
         except ParseException, e:
            errorList.append("PDB model number: " + str(modnum) + ",\n " + e.reason)
	    rval = -1
         except:
            # unexpected
	    rval = -1
         
	 if rval == -1 or rval == -2:
            pass
            """
            try:
               nfile.delete()
               os.unlink(location + curr_fname)
            except:
               pass
            """
         else:
            rlist.append(nfile)
         curr_fname = ""

   if natom > 0:
      # must write out last model
      outfp.write("TER\nEND\n")
      outfp.close()
      nfile = PDBFile()
      nfile.owner = owner
      nfile.location = location
      nfile.filename = curr_fname
      if lessontype == 'lesson1':
         lesson = lesson1.models.Lesson1.objects.filter(id = lessonid)[0]
      else:
         lesson = None
      if lesson:
         nfile.lesson_type = lessontype
         nfile.lesson_id = lessonid
      nfile.save()

      try:
         rval = pdbinfo.views.parser(nfile)
         errorList.append('')
      except ParseException, e:
         if(modnum > 1):
             errorList.append("PDB model number: " + str(modnum) + ", " + e.reason)
	 else:
             errorList.append(e.reason)
	 rval = -1
      except:
         # unexpected
      	 rval = -2
      if rval == -1 or rval == -2:
         try:
            nfile.delete()
            os.unlink(location + curr_fname)
         except:
            pass
      else:
         rlist.append(nfile)

   infp.close()
   os.unlink(location + filename)
   modelinfo = {'models' : rlist, 'errors' : errorList}
   return modelinfo 

