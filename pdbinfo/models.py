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
from django.template import *
from django.template.loader import get_template
from django.db import models
from django import forms
from django.contrib.auth.models import User
from django.core.mail import mail_admins
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from normalmodes.aux import parseNormalModes, getNormalModeMovieNum
import charmming_config, pdbinfo
import commands, datetime, sys, re, os, glob, smtplib
import lesson1, normalmodes, dynamics, minimization
import solvation, lessonaux, apbs
import string, output, charmming_config
import toppar.Top, toppar.Par, lib.Etc

class PDBFile(models.Model):

    owner = models.ForeignKey(User)
    lesson_type = models.CharField(max_length=50,null=True)
    lesson_id = models.PositiveIntegerField(default=0,null=True)

    natom = models.PositiveIntegerField(default=0)
    pdbupload = models.FileField(upload_to='usr/')
    filename = models.CharField(max_length=100)
    segpatch_name = models.CharField(max_length=100)
    patch_name = models.CharField(max_length=100)
    pdb_disul = models.CharField(max_length=100)
    location = models.CharField(max_length=200) 
    title = models.CharField(max_length=250) 
    author = models.CharField(max_length=250) 
    journal = models.CharField(max_length=250)
    pub_date = models.DateTimeField(default=datetime.datetime.now)
    rtf = models.CharField(max_length=250)  
    rtf_append_replace = models.CharField(max_length=250,null=True)  
    prm = models.CharField(max_length=250)  
    prm_append_replace = models.CharField(max_length=250,null=True)  
    selected = models.CharField(max_length=1)  
    good_het = models.CharField(max_length=250)  
    nongood_het = models.CharField(max_length=250)  
    segids = models.CharField(max_length=250)
    pid = models.CharField(max_length=6)
    append_status = models.CharField(max_length=250) 
    solvation_structure = models.CharField(max_length=50) 
    crystal_x = models.DecimalField(max_digits=8,decimal_places=5,null=True)
    crystal_y = models.DecimalField(max_digits=8,decimal_places=5,null=True)
    crystal_z = models.DecimalField(max_digits=8,decimal_places=5,null=True)
    domains = models.CharField(max_length=250,default='')
    fes4list = models.CharField(max_length=250,default='')
    doblncharge = models.CharField(max_length=1,default='f')
    
    minimization_jobID = models.PositiveIntegerField(default=0)
    solvation_jobID = models.PositiveIntegerField(default=0)
    nma_jobID = models.PositiveIntegerField(default=0)
    ld_jobID = models.PositiveIntegerField(default=0)
    md_jobID = models.PositiveIntegerField(default=0)
    sgld_jobID = models.PositiveIntegerField(default=0)
    redox_jobID = models.PositiveIntegerField(default=0)

    def calcJobTime(self,time_filenames,nstep):
        time_filenames = time_filenames.rstrip(',')
        time_list = time_filenames.split(',')
	total_atom_num = 0
	a = .053
	b = 20
	for pdb in time_list:
            total_atom_num += self.countAtomsInPDB(pdb)
	x = total_atom_num 
	time = a*x + b 
        time = time * (int(nstep)/1000)
	time = time/60
        decimal = time - int(time)
        if decimal > 0.5:
          return str(int(time) + 1)
        return str(int(time))
		
    # get a list of disulfide bonds in the PDB (we parse these out in getHeader, this
    # function just returns a nice list for the template to use.
    # Note: this returns a list of tuples
    def getPDBDisulfides(self):
	if not self.pdb_disul:
            return None
        try:
            dfp = open(self.pdb_disul, 'r')
        except:
            return None
        if not dfp:
            return None

        rarr = []
        for line in dfp:
           line = line.strip()
           pres = line.split('-')
           if len(pres) != 2:
               # todo, give an error
               continue
           firstr = pres[0].split(':')
           secondr = pres[1].split(':')
           if len(firstr) != 2 or len(secondr) != 2:
               # todo, give an error
               continue
           rarr.append((firstr[0],firstr[1],secondr[0],secondr[1]))
        dfp.close()
        return rarr

    # This is a dumb way of counting the atoms, but what the hey
    def countAtomsInPDB(self,filename):
        ifp = open(filename, "r")
        if not ifp:
            return 0
        count = 0
        for line in ifp:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                count += 1
        ifp.close()
        return count

    def countAtomsInSeg(self,seg):
        pdbfilename = "new_" + self.stripDotPDB(self.filename) + "-" + seg + ".pdb"
        return self.countAtomsInPDB(pdbfilename)
	
    def countSolvatedAtoms(self):
        try:
            solvation_params = solvation.models.solvationParams.objects.filter(pdb=self,selected='y')[0]
	except:
	    return 0
        if not "Done" in solvation_params.statusHTML:
            return 0
        try:
           solvpdb = "new_" + self.stripDotPDB(self.filename) + "-solvated.pdb"
           os.stat(solvpdb)
        except:
           return 0
        return self.countAtomsInPDB(solvpdb)

    def countNeutralizedAtoms(self):
        solvation_params = solvation.models.solvationParams(pdb=self,selected='y')[0]
        if not "Done" in solvation_params.statusHTML:
            return 0
        try:
           solvpdb = "new_" + self.stripDotPDB(self.filename) + "-neutralized.pdb"
           os.stat(solvpdb)
        except:
           return 0
        return self.countAtomsInPDB(solvpdb)

    def countMinimizedAtoms(self):
        minimize_params = solvation.models.solvationParams(pdb=self,selected='y')[0]
        if not "Done" in minimize_params.statusHTML:
            return 0
        try:
           minpdb = "new_" + self.stripDotPDB(self.filename) + "-min.pdb"
           os.stat(solvpdb)
        except:
           return 0
        return self.countAtomsInPDB(minpdb)
   
    #pre: requires the residue number of the current PDB line
    #pre: requires a list of each PDB line separated by word(whitespace)
    #prints out the PDB in CHARMM-readable format
    def formatLine(self,res_number,linefields2):
        if(res_number < 10):
            try:
                line = " ".join(linefields2)
                printline = '%-6s%5s  %-4s%-4s%s   %-2s    %7s %7s %7s  %4s %4s\n'\
                %(linefields2[0],int(linefields2[1]),linefields2[2],\
                linefields2[3],linefields2[4],linefields2[5],linefields2[6],\
                linefields2[7],linefields2[8],linefields2[9],linefields2[10])
            except:
                 printline = line + '\n'
        elif(res_number < 1000):
            try:
                line = " ".join(linefields2)
                printline = '%-6s%5s  %-4s%-4s%s%4s     %7s %7s %7s  %4s %4s\n'\
                %(linefields2[0],int(linefields2[1]),linefields2[2],\
                linefields2[3],linefields2[4],linefields2[5],linefields2[6],\
                linefields2[7],linefields2[8],linefields2[9], linefields2[10])
            except:
                printline = line + '\n'
        else:
            try:
                line = " ".join(linefields2)
                printline = '%-6s%5s  %-4s%-4s%s %4s    %7s %7s %7s  %4s %4s\n'\
                %(linefields2[0],int(linefields2[1]),linefields2[2],\
                linefields2[3],linefields2[4],linefields2[5],linefields2[6],\
                linefields2[7],linefields2[8],linefields2[9], linefields2[10])
            except:
                printline = line + '\n'
        return printline

    #pre: options must be a string that is minimize or solvate
    #pre: this will remove the minimzie or solvated file
    #pre: appropriately
    #This is used so -f.pdb files aren't listed
    #and only files that solvation/minimization take
    #advantage of are returned

    #Options allow you to exclude items from the list
    def getLimitedFileList(self,options):
        #file_list = self.getFileList()
        dash_f = re.compile('-final')
	#solvated = re.compile('-solv')
	#md = re.compile('-md')
	#ld = re.compile('-ld')
	#sgld = re.compile('-sgld')
	minimized = re.compile('-min')
	seg_list = self.segids.split()
	het_list = self.good_het.split()
	tip_list = self.nongood_het.split()
	os.chdir(self.location)
	protein_list = []

	for item in seg_list + het_list + tip_list:
	    try:
	        if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb"):
		    protein_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
            except:
                try:
                    if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb"):
                        protein_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                except:
                    pass

	try:
	    if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-min" + ".pdb") and options != "minimization":
	        protein_list.append("new_" + self.stripDotPDB(self.filename) + "-min" +".pdb")
        except:
	    pass
	try:
	    if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-solv" + ".pdb") and options != "solvation":
	        protein_list.append("new_" + self.stripDotPDB(self.filename) + "-solv" +".pdb")
        except:
	    pass
	try:
	    if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-neutralized.pdb") and options != "solvation":
	        protein_list.append("new_" + self.stripDotPDB(self.filename) + "-neutralized.pdb")
        except:
	    pass
	try:
	    if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-md" + ".pdb"):
	        protein_list.append("new_" + self.stripDotPDB(self.filename) + "-md" +".pdb")
        except:
	    pass
	try:
	    if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-mdavg" + ".pdb"):
	        protein_list.append("new_" + self.stripDotPDB(self.filename) + "-mdavg" +".pdb")
        except:
	    pass
	try:
	    if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-ld" + ".pdb"):
	        protein_list.append("new_" + self.stripDotPDB(self.filename) + "-ld" + ".pdb")
        except:
	    pass
	try:
	    if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-sgld" + ".pdb"):
	        protein_list.append("new_" + self.stripDotPDB(self.filename) + "-sgld" + ".pdb")
        except:
	    pass
	return protein_list
   

   
    #This returns every single file that has been made through the interface
    def getAllFiles(self):
        file_list = []
        seg_list = self.segids.split(' ')
        tip_list = self.good_het.split(' ')
        het_list = self.nongood_het.split(' ')
        done = re.compile('Done')
        fail = re.compile('Failed')
        for item in seg_list:
            try:
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                temphandle.close()
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.psf")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.crd")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".out")
                temphandle.close()
            except:
                pass
        for item in tip_list:
            try:
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                temphandle.close()
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.psf")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.crd")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".out")
                temphandle.close()
            except:
                pass
        for item in het_list:
            try:
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.psf")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.crd")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".out")
                temphandle.close()
            except:
               pass
        if done.search(self.append_status):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-final.pdb")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-final.psf")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-final.crd")
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-append.inp")
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-append.out")
	try:
            minimize_params = solvation.models.solvationParams(pdb=self,selected='y')[0]
	except:
	    minimize_params = ''
        if done.search(minimize_params.statusHTML):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-min.pdb")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-min.psf")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-min.crd")
        if fail.search(minimize_params.statusHTML) or done.search(minimize_params.statusHTML):
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-min.inp")
            try:
                os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-min.out")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-min.out")
            except:
                pass
	try:
	    solvation_params = solvation.models.solvationParams.objects.filter(pdb = self, selected = 'y')[0]	
	except:
	    solvation_params = ''
        if done.search(solvation_params.statusHTML):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-solv.pdb")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-solv.psf")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-solv.crd")
            try:
                os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-neutralized.crd")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-solv.pdb")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-solv.psf")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-solv.crd")
            except:
                pass
        if fail.search(solvation_params.statusHTML) or done.search(solvation_params.statusHTML):
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-solv.inp")
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-solv.out")
            file_list.append("/solvation/water.crd")
            file_list.append("/scripts/savegv.py")
            file_list.append("/scripts/savechrg.py")
            file_list.append("/scripts/calcewald.pl")
            try:
                os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-neutralize.inp")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-neutralize.inp")
                os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-neutralized.pdb")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-neutralized.pdb")
                os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-neutralized.psf")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-neutralized.psf")
                os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-neutralized.crd")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-neutralized.crd")
                file_list.append("loop-" + self.stripDotPDB(self.filename) + ".log")
            except:
                pass
            try:
                os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-neutralize.out")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-neutralize.out")
            except:
                pass
	try:
	    md_status = dynamics.models.mdParams.objects.filter(pdb = self, selected='y')[0].statusHTML
	except:
	    md_status = ''
        if done.search(md_status):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-md.pdb")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-md.psf")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-md.crd")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-md.dcd")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-md.res")
        if fail.search(md_status) or done.search(md_status):
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-md.inp")
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-md.out")
	try:
	    ld_status = dynamics.models.ldParams.objects.filter(pdb = self, selected='y')[0].statusHTML
	except:
	    ld_status = ''
        if(ld_status):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-ld.pdb")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-ld.psf")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-ld.crd")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-ld.dcd")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-ld.res")
        if fail.search(ld_status) or done.search(ld_status):
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-ld.inp")
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-ld.out")
	try:
	    sgld_status = dynamics.models.sgldParams.objects.filter(pdb = self, selected='y')[0].statusHTML
	except:
	    sgld_status = ''
        if(done.search(sgld_status)):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-sgld.pdb")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-sgld.psf")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-sgld.crd")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-sgld.dcd")
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-sgld.res")
        if fail.search(sgld_status) or done.search(sgld_status):
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-sgld.inp")
            file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-sgld.out")
	try:
	    md_movie_status = dynamics.models.mdParams.objects.filter(pdb = self, selected='y')[0].md_movie_status
	except:
	    md_movie_status = ''
        if done.search(md_movie_status):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-md-mainmovie.pdb")
	try:
	    ld_movie_status = dynamics.models.mdParams.objects.filter(pdb = self, selected='y')[0].ld_movie_status
	except:
	    ld_movie_status = ''
        if done.search(ld_movie_status):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-ld-mainmovie.pdb")
	try:
	    sgld_movie_status = dynamics.models.mdParams.objects.filter(pdb = self, selected='y')[0].sgld_movie_status
	except:
	    sgld_movie_status = ''
        if done.search(sgld_movie_status):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-sgld-mainmovie.pdb")

        return file_list

    # Updates the status of in progress operations
    def updateActionStatus(self):
        done = re.compile('Done')
        fail = re.compile('Failed')
        si = schedInterface()        

        if self.minimization_jobID != 0:
            sstring = si.checkStatus(self.minimization_jobID)
	    miniparam_obj = minimization.models.minimizeParams.objects.filter(pdb=self,selected='y')[0]
            miniparam_obj.statusHTML = statsDisplay(sstring,self.minimization_jobID)
	    miniparam_obj.save()
            if done.search(miniparam_obj.statusHTML) or fail.search(miniparam_obj.statusHTML):
               self.minimization_jobID = 0
               if self.lesson_type:
                   lessonaux.doLessonAct(self,'onMinimizeDone')
        if self.solvation_jobID != 0:
            sstring = si.checkStatus(self.solvation_jobID)
	    solvparam_obj = solvation.models.solvationParams.objects.filter(pdb=self,selected='y')[0]
            solvparam_obj.statusHTML = statsDisplay(sstring,self.solvation_jobID)
	    solvparam_obj.save()
            if done.search(solvparam_obj.statusHTML):
               self.solvation_jobID = 0
               if self.lesson_type:
                   lessonaux.doLessonAct(self,'onSolvationDone')
        if self.nma_jobID != 0:
            sstring = si.checkStatus(self.nma_jobID)
	    nmaparam_obj = normalmodes.models.nmodeParams.objects.filter(pdb=self,selected='y')[0]
            nmaparam_obj.statusHTML = statsDisplay(sstring,self.nma_jobID)
	    nmaparam_obj.save()
            nma_status = statsDisplay(sstring,self.nma_jobID)
            if done.search(nma_status):
               self.nma_jobID = 0
	       parseNormalModes(self)
               if nmaparam_obj.nma_movie_req:
	           nmaparam_obj.make_nma_movie = True
                   nmaparam_obj.save()
        if self.md_jobID != 0:
            sstring = si.checkStatus(self.md_jobID)
	    mdparam_obj = dynamics.models.mdParams.objects.filter(pdb=self,selected='y')[0]
            mdparam_obj.statusHTML = statsDisplay(sstring,self.md_jobID)
	    mdparam_obj.save()
	    #self.md_status = statsDisplay(sstring)
            if done.search(mdparam_obj.statusHTML) and mdparam_obj.md_movie_req:
               self.md_jobID = 0
	       mdparam_obj.make_md_movie = True
               mdparam_obj.save()
	       self.save()
               if self.lesson_type:
                  lessonaux.doLessonAct(self,'onMDDone')
            elif done.search(mdparam_obj.statusHTML):
               if self.lesson_type:
                  lessonaux.doLessonAct(self,'onMDDone')
               self.md_jobID = 0
               self.save()
        if self.ld_jobID != 0:
            sstring = si.checkStatus(self.ld_jobID)
	    ldparam_obj = dynamics.models.ldParams.objects.filter(pdb=self,selected='y')[0]
            ldparam_obj.statusHTML = statsDisplay(sstring,self.ld_jobID)
	    ldparam_obj.save()
            if done.search(ldparam_obj.statusHTML) and ldparam_obj.ld_movie_req:
	       ldparam_obj.make_ld_movie = True
               ldparam_obj.save()
               self.ld_jobID = 0
               if self.lesson_type:
                  lessonaux.doLessonAct(self,'onLDDone')
            elif done.search(ldparam_obj.statusHTML):
               if self.lesson_type:
                  lessonaux.doLessonAct(self,'onLDDone')
               self.ld_jobID = 0
               self.save()
        if self.sgld_jobID != 0:
            sstring = si.checkStatus(self.sgld_jobID)
	    sgldparam_obj = dynamics.models.sgldParams.objects.filter(pdb=self,selected='y')[0]
            sgldparam_obj.statusHTML = statsDisplay(sstring,self.sgld_jobID)
	    sgldparam_obj.save()
            self.sgld_status = statsDisplay(sstring,self.sgld_jobID)
            if done.search(sgldparam_obj.statusHTML) and sgldparam_obj.sgld_movie_req:
	       sgldparam_obj.make_sgld_movie = True
               sgldparam_obj.save()
               if self.lesson_type:
                  lessonaux.doLessonAct(self,'onSGLDDone')
               self.sgld_jobID = 0
            elif done.search(sgldparam_obj.statusHTML):
               if self.lesson_type:
                  lessonaux.doLessonAct(self,'onSGLDDone')
               self.sgld_jobID = 0
               self.save()
        if self.redox_jobID != 0:
            sstring = si.checkStatus(self.redox_jobID)
            redox_obj = apbs.models.redoxParams.objects.filter(pdb=self,selected='y')[0]
            redox_obj.statusHTML = statsDisplay(sstring,self.redox_jobID)
            redox_obj.save()
            # ToDo, add lesson hooks
    
        self.save()

    # Returns a list of PDB files with -final
    def getDashFPDBList(self):
        file_list = []
	total_list = []
        seg_list = self.segids.split(' ')
	het_list = self.nongood_het.split(' ')
	tip_list = self.good_het.split(' ')
	big_list = seg_list + tip_list + het_list
        for item in big_list:
            try:
    	        os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
    	        file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
	    except:
	        pass
	return file_list
  
    #Returns a list of PDB protein segments
    def getProteinSegPDBList(self):
        file_list = []
        seg_list = self.segids.split(' ')
        for item in seg_list:
  	    fname = "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb"
            try:
    	        os.stat(self.location + fname)
    	        file_list.append(fname)
	    except:
                fname = "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb"
                try:
                    os.stat(self.location + fname)
                    file_list.append(fname)
                except:
                    pass

        return file_list
    
    #Returns a list of Goodhet PDBs
    def getGoodHetPDBList(self):
        file_list = []
        tip_list = self.good_het.split(' ')
        for item in tip_list:
            try:
    	        os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
    	        file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
	    except:
                try:
             	    os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                    file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                except:
	            pass
	return file_list
   
    #Returns a list of NonGoodhet PDBs
    def getNonGoodHetPDBList(self):
        file_list = []
        het_list = self.nongood_het.split(' ')
        for item in het_list:
            try:
    	        os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
    	        file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
	    except:
                try:
                    os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                    file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                except:
	            pass
	return file_list
  
    #Gets basic protein segment file list for the download files page
    def getDownloadSegmentList(self):
        file_list = []
        seg_list = self.segids.split(' ')
        tip_list = self.good_het.split(' ')
        het_list = self.nongood_het.split(' ')
        done = re.compile('Done')
        fail = re.compile('Failed')
        for item in seg_list:
            try:
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                temphandle.close()
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.psf")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.crd")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".out")
                temphandle.close()
            except:
                pass
        for item in tip_list:
            try:
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                temphandle.close()
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.psf")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.crd")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".out")
                temphandle.close()
            except:
                pass
        for item in het_list:
            try:
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
                temphandle = open(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-f.pdb",'r')
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.psf")
                file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.crd")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
                file_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".out")
                temphandle.close()
            except:
               pass
        return file_list
    
    #Returns a list of files not specifically associated with the structure
    def getNonProteinFiles(self):
        file_list = []
	file_list.append("/solvation/water.crd")
        file_list.append("/scripts/savegv.py")
        file_list.append("/scripts/savechrg.py")
        file_list.append("/scripts/calcewald.pl")
        return file_list

    #Gets list of PDBs that have been processed
    def getFileList(self):

        file_list = []
        seg_list = self.segids.split(' ')
        tip_list = self.good_het.split(' ')
        het_list = self.nongood_het.split(' ')

        done = re.compile('Done')
        for item in seg_list + het_list + tip_list:
            try:
                if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb"):
                    file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + "-final.pdb")
            except:
                try:
                    if os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb"):
                       	file_list.append("new_" + self.stripDotPDB(self.filename) + "-" + item + ".pdb")
               	except:
                    pass

	if(done.search(self.append_status)):
          file_list.append("new_" + self.stripDotPDB(self.filename) + "-final.pdb")
	try:
	    minimize_params = minimization.models.minimizeParams.objects.filter(pdb=self,selected='y')[0].statusHTML
	except:
	    minimize_params = ''
        if(done.search(minimize_params)):
          file_list.append("new_" + self.stripDotPDB(self.filename) + "-min.pdb")
	try:
	    solvation_params = solvation.models.solvationParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    solvation_params = ''
        if(done.search(solvation_params)):
           file_list.append("new_" + self.stripDotPDB(self.filename) + "-solv.pdb")
           try:
              os.stat(self.location + "new_" + self.stripDotPDB(self.filename) + "-neutralized.pdb")
              file_list.append("new_" + self.stripDotPDB(self.filename) + "-neutralized.pdb")
           except:
              pass
	try:
	    md_status = dynamics.models.mdParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    md_status = ''
        if(done.search(md_status)):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-md.pdb")
	
        try:
	    nma_status = normalmodes.models.nmodeParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    nma_status = ''
	if(done.search(nma_status)):
            nmtrjnum = getNormalModeMovieNum(self)
            for trj in range(nmtrjnum):
                trj += 1
                file_list.append('new_' + self.stripDotPDB(self.filename) + '-nma-mainmovie-' + str(trj) + '.pdb')
	try:
	    md_movie_status = dynamics.models.mdParams.objects.filter(pdb = self, selected = 'y')[0].md_movie_status
	except:
	    md_movie_status = ''
        if(md_movie_status):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-md-mainmovie.pdb")
	try:
	    ld_status = dynamics.models.ldParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    ld_status = ''
        if(done.search(ld_status)):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-ld.pdb")
	try:
	    ld_movie_status = dynamics.models.ldParams.objects.filter(pdb = self, selected = 'y')[0].ld_movie_status
	except:
	    ld_movie_status = ''
        if(ld_movie_status):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-ld-mainmovie.pdb")
	try:
	    sgld_status = dynamics.models.sgldParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    sgld_status = ''
        if(done.search(sgld_status)):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-sgld.pdb")
	try:
	    sgld_movie_status = dynamics.models.sgldParams.objects.filter(pdb = self, selected = 'y')[0].sgld_movie_status
	except:
	    sgld_movie_status = ''
        if(sgld_movie_status):
            file_list.append("new_" + self.stripDotPDB(self.filename) + "-sgld-mainmovie.pdb")
        
        return file_list 
    
    #Checks for all input files that may have been created
    #and returns it as a list
    def getInputList(self):
        input_list = []
        seg_list = self.segids.split()
        het_list = self.nongood_het.split(' ')
        tip_list = self.good_het.split(' ')
        done = re.compile('Done')
        fail = re.compile('Failed')

	for item in het_list:
	    try:
    	        os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
    	        input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
    	        input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".out")
	    except:
	        pass
        for item in seg_list:
            try:
    	        os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
    	        input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
    	        input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".out")
    	    except:
    	        pass
        for item in tip_list:
            try:
    	        os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
    	        input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".inp")
    	        input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-" + item + ".out")
    	    except:
    	        pass
        try:
    	    os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-energy" + ".inp")
    	    input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-energy" + ".inp")
    	    input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-energy" + ".out")
    	except:
    	    pass
	if(done.search(self.append_status)):
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-append.inp")
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-append.out")
	try:
	    minimize_status = minimization.models.minimizeParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    minimize_status = ''
        if done.search(minimize_status) or fail.search(minimize_status):
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-min.inp")
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-min.out")
	try:
	    solvation_status = solvation.models.solvationParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    solvation_status = ''
        if done.search(solvation_status) or fail.search(solvation_status):
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-solv.inp")
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-solv.out")
            try:
                os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-neutralize.inp")
                input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-neutralize.inp")
            except:
                pass
            try:
                os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-neutralize.out")
                input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-neutralize.out")
            except:
                pass
	try:
	    nma_status = normalmodes.models.nmodeParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    nma_status = ''
        if done.search(nma_status) or fail.search(nma_status):
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-nmodes.inp")
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-nmodes.out")
        #get parameter and check status
	try:
	    md_status = dynamics.models.mdParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    md_status = ''
        if done.search(md_status) or fail.search(md_status):
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-md.inp")
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-md.out")
            
	try:
	    ld_status = dynamics.models.ldParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    ld_status = ''
        if done.search(ld_status) or fail.search(ld_status):
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-ld.inp")
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-ld.out")
            
	try:
	    sgld_status = dynamics.models.sgldParams.objects.filter(pdb = self, selected = 'y')[0].statusHTML
	except:
	    sgld_status = ''
        if done.search(sgld_status) or fail.search(sgld_status):
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-sgld.inp")
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-sgld.out")

        try:
            os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-rmsd.inp")
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-rmsd.inp")
            os.stat(self.location + "charmm-" + self.stripDotPDB(self.filename) + "-rmsd.out")
            input_list.append("charmm-" + self.stripDotPDB(self.filename) + "-rmsd.out")
        except:
            pass

        loc = os.getcwd()
        os.chdir(self.location)
        input_list += glob.glob("redox-%s-*.inp" % self.stripDotPDB(self.filename))
        input_list += glob.glob("redox-%s-*.out" % self.stripDotPDB(self.filename))
        os.chdir(loc)

        return input_list 

    #takes the information from the Remark statement of
    #a PDB and determines the title, jrnl, and author
    def getHeader(self):
        dfp = None
        dspatchfile = self.location + "/pdb_disulfides-" + self.stripDotPDB(self.filename) + ".patch"

        remark = re.compile('REMARK')
        title = re.compile('TITLE')
        jrnl = re.compile('JRNL')
        author = re.compile('AUTHOR') 
        ref = re.compile('REF') 
        refn = re.compile('REFN')
        temp = open(self.location + "/" + self.filename, 'r')
        for line in temp:
	    #ignore remark statements
            if(remark.match(line)):            
                break
            if(title.match(line)):
                self.title = self.title + line.strip()
            elif(author.match(line)):
                self.author = self.author + line.strip()
            elif(jrnl.match(line) and (ref.search(line) or refn.search(line))):
                self.journal = self.journal + line.strip()
            elif line.startswith("SSBOND"):
                # process disulfide bridges
                if not dfp:
                    dfp = open(dspatchfile, 'w')
                    self.pdb_disul = dspatchfile
                    self.save()
                tlist = string.splitfields(line)
                try:
                    dfp.write('%s:%s-%s:%s\n' % (tlist[3].lower(), tlist[4], tlist[6].lower(), tlist[7]))
                except:
                    # ToDo, we should raise some sort of parse error here, but for now just pass
                    pass
        if dfp:
            dfp.close()
	if(self.title):
            self.title = (title.sub('',self.title))
            if len(self.title) > 249:
                self.title = self.title[0:248]
	else:
	    self.title = "No information found"
	if(self.author):
            self.author = (author.sub('',self.author))
            if len(self.author) > 249:
                self.author = self.author[0:248]
	else:
	    self.author = "No information found"
	self.journal = self.journal.strip()
	if(self.journal):
            self.journal = (jrnl.sub('',self.journal))
            self.journal = (ref.sub('',self.journal))
            self.journal = (refn.sub('',self.journal))
            if len(self.journal) > 249:
                self.journal = self.journal[0:248]
	else:
	    self.journal = "No information found"
        temp.close()
 
    #pre: 
    #post: Returns 1 if there is only  a custom 
    #topology file, 2 if there is only a custom 
    #parameter, 3 if both, -1 if none
    def ifExistsRtfPrm(self):
        if self.rtf and self.prm:
            if self.prm.startswith("genrtf") or self.prm.startswith("antechamber"):
                return -1
	    return 3
	if self.rtf:
            if self.rtf.startswith("genrtf") or self.rtf.startswith("antechamber"):
                return -1
	    return 1
	if self.prm:
            if self.prm.startswith("genrtf") or self.prm.startswith("antechamber"):
                return -1
	    return 2
	else:
	    return -1


    #pre:requires a request object
    #checks if there was a topology or parameter file
    #and if there was it will store the information
    def getRtfPrm(self,request,makeGoModel,makeBLNModel):

        if makeGoModel:
            gmp = goModel()
            gmp.pdb = self

            startdir = os.getcwd()
            os.chdir(self.location)
            nGoSegs = 0

            # check POSTDATA
            domainData = None
            try:
                uploadType = request.POST['gm_dm_uploadtype']
            except:
                uploadType = None
            if uploadType == 'auto':
                gmp.domainData = ''
            elif uploadType == 'box':
                try:
                    gmp.domainData = request.POST['gm_dm_string']
                except:
                    gmp.domainData = ''

            try:
                gmp.contactType = request.POST['gm_contact_type']
            except:
                gmp.contact_type = 'bt'

            try:
                gmp.nScale = request.POST['gm_nScale']
            except:
                gmp.nScale = '0.05'
            try:
                gmp.domainScale = request.POST['gm_domainScale']
            except:
                gmp.domainScale = '1.0'
            try:
                gmp.kBond = request.POST['gm_kBond']
            except:
                gmp.kBond = '50.0'
            try:
                gmp.kAngle = request.POST['gm_kAngle']
            except:
                gmp.kAngle = '30.0'

            # for now, there's only one Go-model object with a given PDB, so this
            # is safe ... in time we'll need to refine this scheme
            gmp.selected = 'y'
            gmp.save()

            # now run the command...
            seg_list = self.segids.split(' ')
            go_list = []
            for segment in seg_list:
                inpFilename = "new_%s-%s.pdb" % (self.stripDotPDB(self.filename),segment)

                # we only build Go models for protein segments
                if not segment.endswith('-pro'):
                    # hang on to non-protein segs in case the user wants to create a mixed model
                    #seg_list.remove(segment)
                    #os.unlink(inpFilename)
                    continue

                # For each individual segment, we don't use any domain data that we might have available since
                # not all domains will be in this segment. Instead, if there's domain data we create a new
                # topology/parameter file below with the global information.
                mycmd = "PYTHONPATH=%s %s/cg_models/alpha_carbon/alpha_carbon.py --pdb_only --input=%s --data=%s/etc/cg --output=%s --stride_path=/usr/local/bin --contacts=%s --nScale=%s" \
                        % (charmming_config.charmming_root,charmming_config.charmming_root,inpFilename,charmming_config.charmming_root,self.location,gmp.contactType,gmp.nScale)
                status, output = commands.getstatusoutput(mycmd)
                if status != 0:
                    raise "go model creation failed!"

                # re-name the go-ified PDB files
                cgInpFilename = "cg_%s" % inpFilename
                cgOutFilename = "new_%s-%s.pdb" % (self.stripDotPDB(self.filename),segment.replace('-pro','-go'))
                os.rename(cgInpFilename, cgOutFilename)
                go_list.append(segment.replace('-pro','-go'))
                nGoSegs += 1

            # re-save segids in case we deleted any of them, also get rid of the good het and nongood het segments
            # since we can't use them with the CG model
            seg_list += go_list
            self.segids = ' '.join(seg_list)

            self.save()
            os.chdir(startdir)
            return

        if makeBLNModel:
            blnp = blnModel()
            startdir = os.getcwd()
            os.chdir(self.location)
            nBLNSegs = 0

            domainData = None
            try:
                uploadType = request.POST['bln_dm_uploadtype']
            except:
                uploadType = None
            if uploadType == 'auto':
                blnp.domainData = ''
            elif uploadType == 'box':
                try:
                    blnp.domainData = request.POST['bln_dm_string']
                except:
                    blnp.domainData = ''

            try:
                blnp.nScale = request.POST['bln_nScale']
            except:
                blnp.nScale = '1.0'
            try:
                blnp.domainScale = request.POST['bln_domainScale']
            except:
                blnp.domainScale = '1.0'


            try:
                blnp.kBondHelix = request.POST['bln_kBondHelix']
            except:
                blnp.kBondHelix = '3.5'
            try:
                blnp.kBondSheet = request.POST['bln_kBondSheet']
            except:
                blnp.kBondSheet = '3.5'
            try:
                blnp.kBondCoil = request.POST['bln_kBondCoil']
            except:
                blnp.kBondCoil = '2.5'

            try:
                blnp.kAngleHelix = request.POST['bln_kAngleHelix']
            except:
                blnp.kAngleHelix = '8.37'
            try:
                blnp.kAngleSheet = request.POST['bln_kAngleSheet']
            except:
                blnp.kAngleSheet = '8.37'
            try:
                blnp.kAngleCoil = request.POST['bln_kAngleCoil']
            except:
                blnp.kAngleCoil = '5.98'

            blnp.save()

            # now run the command
            seg_list = self.segids.split(' ')
            bln_list = []
            for segment in seg_list:
                inpFilename = "new_%s-%s.pdb" % (self.stripDotPDB(self.filename),segment)

                # we only build Go models for protein segments
                if not segment.endswith('-pro'):
                    continue

                # For each individual segment, we don't use any domain data that we might have available since
                # not all domains will be in this segment. Instead, if there's domain data we create a new
                # topology/parameter file below with the global information.
                mycmd = "PYTHONPATH=%s %s/cg_models/bln/bln.py --input=%s --data=%s/etc/cg --output=%s --stride_path=/usr/local/bin --nScale=%s --domainScale=%s --kBondHelix=%s --kBondSheet=%s --kBondCoil=%s --kAngleHelix=%s --kAngleSheet=%s --kAngleCoil=%s" \
                        % (charmming_config.charmming_root,charmming_config.charmming_root,inpFilename,charmming_config.charmming_root,self.location,blnp.nScale,blnp.domainScale,blnp.kBondHelix,blnp.kBondSheet,blnp.kBondCoil, \
                           blnp.kAngleHelix,blnp.kAngleSheet,blnp.kAngleCoil)
                status, output = commands.getstatusoutput(mycmd)
                if status != 0:
                    raise "go model creation failed!"

                # re-name the bln-ified PDB files
                cgInpFilename = "bln_%s" % inpFilename
                cgOutFilename = "new_%s-%s.pdb" % (self.stripDotPDB(self.filename),segment.replace('-pro','-bln'))
                os.rename(cgInpFilename, cgOutFilename)
                bln_list.append(segment.replace('-pro','-bln'))
                nBLNSegs += 1

            seg_list += bln_list
            self.segids = ' '.join(seg_list)
            self.save()

            os.chdir(startdir)
            return

        postfiles = request.FILES
        if postfiles.has_key('rtf_file'):
            self.rtf = "rtf-" + self.stripDotPDB(self.filename) + ".rtf"
            #Check to see if the user uploaded an rtf and/or prm and if they did
            #See if they wanted it appended or replaced
            self.rtf_append_replace = request.POST['rtf_append_or_replace']
	    self.save()
            temp = open(self.location + self.rtf,'w') 
	    for fchunk in postfiles['rtf_file'].chunks():
                temp.write(fchunk)
	    temp.close()
	else:
	    self.rtf = ""

        if postfiles.has_key('prm_file'):
            self.prm = "prm-" + self.stripDotPDB(self.filename) + ".prm"
            #Check to see if the user uploaded an rtf and/or prm and if they did
            #See if they wanted it appended or replaced
            self.prm_append_replace = request.POST['prm_append_or_replace']
            temp = open(self.location + self.prm,'w') 
            for fchunk in postfiles['prm_file'].chunks():
                temp.write(fchunk)
	    temp.close()
	    self.save()
	else:
	    self.prm = ""

    #pre:requires a reference to itself
    #Returns the path of a rtf/prm file depending on
    #whether the user uploaded one of not
    def getRtfPrmPath(self):
        rtf_prm_dict = {}
	# in the case where genRTF has been made, it should be appended and the RTF/PRM
	# should not be changed either
	rtf_prm_dict["rtf"] = "/usr/local/charmming/toppar/top_all27_prot_na.rtf"
        rtf_prm_dict["prm"] = "/usr/local/charmming/toppar/par_all27_prot_na.prm"
	return rtf_prm_dict

    #pre:requires a string
    #If a line contains segid-goodhet, this function will
    #return the segid
    def segidFromGoodhet(self,line_with_goodhet):
        strip_goodhet = re.compile('-goodhet')
	list_without_goodhet = strip_goodhet.sub('',line_with_goodhet)
	line_without_goodhet = ' '.join(list_without_goodhet)
	return line_without_goodhet

    #pre: requires a string
    #This will return a string of the filename without
    #the .crd
    def stripDotCOR(self, old_filename):
        new_filename = old_filename
	pdb = re.compile('.cor')
        new_filename = 	pdb.sub('',new_filename)
        return new_filename
    
    #pre: requires a string
    #This will return a string of the filename without
    #the .crd
    def stripDotCRD(self, old_filename):
        new_filename = old_filename
	pdb = re.compile('.crd')
        new_filename = 	pdb.sub('',new_filename)
        return new_filename
    
    #pre: requires a string
    #This will return a string of the filename without
    #the .pdb
    def stripDotPDB(self, old_filename):
        new_filename = old_filename
	pdb = re.compile('.pdb')
        new_filename = 	pdb.sub('',new_filename)
	return new_filename

    #pre: requires a string
    #This will return a string of the filename without
    #the prefix new_
    def stripNew(self, old_filename):
        new_filename = old_filename
	new = re.compile('new_')
        new_filename = 	new.sub('',new_filename)
	return new_filename
	
    #pre: requires a filename in the format new_temp-seg.pdb
    #This strips out the new_temp-
    def getSegIdFromFilename(self, old_filename):
        tempregex = re.compile('new_\w*-\d*-')
        new_filename = 	tempregex.sub('',old_filename)
	new_filename = self.stripDotPDB(new_filename)
        nfl = new_filename.split("-")[1:]
        if len(nfl) > 1:
           if "het" == nfl[1] or "goodhet" == nfl[1] or "pro" == nfl[1] or "rna" == nfl[1] or "dna" == nfl[1] or "go" == nfl[1] or "bln" == nfl[1]:
              return '-'.join(nfl[:2])
           else:
              return nfl[0]
        else:
           return nfl[0]

    # Tim Miller, make hetid RTF/PRM using antechamber
    def makeAntechamber(self,hetids,doHyd):
        os.putenv("ACHOME", "/usr/local/charmming/antechamber")
        cwd = os.getcwd()
        os.chdir(self.location)
        for het in hetids:
            fname = "new_" + self.stripDotPDB(self.filename) + "-" + het + ".pdb"
            acbase = "antechamber-" + self.stripDotPDB(self.filename) + "-" + het 
            rgbase = "new_" + self.stripDotPDB(self.filename) + "-" + het 
            acname = acbase + ".ac"

            if doHyd:
                os.system("/usr/local/bin/babel -ipdb " + fname + " -opdb " + fname + " -h")

            cmd = "/usr/local/charmming/antechamber/exe/antechamber -fi pdb -i " + \
                  fname + " -fo ac -o " + acname 
            status = os.system(cmd)
            if status != 0:
                # debug purposes only
                raise("Antechamber screwed up!")

            try:
                # see if output file exists
                os.stat(acname)
            except:
                # nope, run antechamber w/ -j 5
                cmd = "/usr/local/charmming/antechamber/exe/antechamber -j 5 -fi pdb -i " + \
                      fname + " -fo ac -o " + acname
                status = os.system(cmd)
                if status != 0:
                    # debug purposes only
                   raise("Antechamber screwed up!")

            cmd = "/usr/local/charmming/antechamber/exe/charmmgen -f ac -i " + \
                  acname + " -o " + acbase + " -s " + het
            status = os.system(cmd)
            if status != 0:
                # also for debug purposes only
                raise("Charmmgen screwed up!")
            # run through CHARMM since Antechamber changes the residue IDs
            cmd = "/usr/local/charmming/gfortran-xxlg-qc.one < " + acbase + ".inp > " + \
                  acbase + ".out"
            os.system(cmd)
            os.rename(acbase + ".psf", rgbase + "-final.psf")
            os.rename(acbase + ".pdb", rgbase + "-final.pdb")
            os.rename(acbase + ".crd", rgbase + "-final.crd")
            os.system("rm -f ANTECHAMBER* ATOMTYPE.INF")

            # set the newly generated RTF/PRM to be used.
            self.rtf = acbase + ".rtf"
            self.prm = acbase + ".prm"

        os.chdir(cwd)    

    #pre:requires a list of het segids
    #This will run genRTF through the non-compliant hetatms
    #and make topology files for them
    def makeGenRTF(self,hetids,doHyd):
        cwd = os.getcwd()
	os.chdir(self.location)

        for het in hetids:
            filebase = "new_"+ self.stripDotPDB(self.filename) + "-" + het 
            filename = filebase + '.pdb'

            if doHyd:
                os.system('/usr/local/bin/babel -ipdb ' + filename + ' -oxyz ' + filebase + '.xyz -h')
                os.system('/usr/local/charmming/genrtf-v3.3 -s ' + het + ' -x ' + filebase + '.xyz')

                # Run the GENRTF generated inp script through CHARMMing b/c the residue names/ids have changed.
                # so we can't trust the original PDB any more.
                os.system('/usr/local/charmming/gfortran-xxlg-qc.one < ' + filebase + '.inp > ' + filebase + '.out')
            else:
	        os.system('/usr/local/charmming/genrtf-v3.3 ' +filename)

	    continue_to_write = ""

	    #The new rtf filename will look like genrtf-filename.rtf
	    self.rtf = "genrtf-" + self.stripNew(self.stripDotPDB(filename)) + ".rtf"
	    self.save()
	    genrtf_handle = open(self.location + self.stripDotPDB(filename) + ".inp",'r')
	    write_rtf_handle = open(self.location + "genrtf-" +\
	                       self.stripNew(self.stripDotPDB(filename)) + ".rtf",'w')

	    #Now the input file will be parsed for the rtf lines
	    #The RTF files start with read rtf card append and end with End
	    write_rtf_handle.write('* Genrtf \n*\n')
            write_rtf_handle.write('   28\n\n') # need the CHARMM version under which this was generated, c28 is OK.
	    mass = re.compile("MASS")
	    end = re.compile("End")
	    for line in genrtf_handle:
	        if(mass.search(line)):
		    write_rtf_handle.write(line)
		    continue_to_write = 1
		elif(end.search(line)):
		    write_rtf_handle.write(line)
		    continue_to_write = 0
		    break
	        elif(continue_to_write):
		    write_rtf_handle.write(line)
	    write_rtf_handle.close()

	    #now for the Paramter files
	    continue_to_write = ""
	    #The new prm filename will look like genrtf-filename.prm
	    self.prm = "genrtf-" + self.stripNew(self.stripDotPDB(filename)) + ".prm"
	    self.save()
	    write_prm_handle = open(self.location + self.prm,'w')
	    #Now the input file will be parsed for the prm lines
	    #The PRM files start with read prm card append and end with End
	    write_prm_handle.write('* Genrtf \n*\n\n')
	    bonds_dont_believe = re.compile("bonds ! Don't believe these ones")
	    for line in genrtf_handle:
	        if(bonds_dont_believe.search(line)):
		    write_prm_handle.write(line)
		    continue_to_write = 1
		elif(end.search(line)):
                    write_prm_handle.write(line)
		    continue_to_write = 0
		    break
		elif(continue_to_write):
                    write_prm_handle.write(line)
	    
            genrtf_handle.close()
	    write_prm_handle.close()

        os.chdir(cwd)

    
    #pre: requires a list of segids
    #post: reveals if there is already a basic input for 
    #the segids. This is done so a person can solvate
    #before minimization
    
    def ifBasicInputExists(self,seg):
        cpath = self.location + 'new_' + self.stripDotPDB(self.filename) + '-' + seg + '-final.crd'
        ppath = self.location + 'new_' + self.stripDotPDB(self.filename) + '-' + seg + '-final.psf'
	try:
	   os.stat(cpath)
	   os.stat(ppath)
	except:
	   return False
	return True

    #pre: Requires a valid CHARMM title
    #This will make and return the header files for a CHARMM file such as the paramter/topology
    #paths and CHARMM title
    def makeCHARMMInputHeader(self, charmm_title, postdata):
	rtf_prm_dict = self.getRtfPrmPath()
        if self.rtf_append_replace == None or self.rtf_append_replace == 'append':
            charmm_inp = """* """ + charmm_title + """
*
bomlev -2

! Read in Topology and  Parameter files

open unit 1 card read name """ + rtf_prm_dict['rtf'] + """ 
read RTF card unit 1
close unit 1"""
            if(self.rtf):
	        charmm_inp = charmm_inp + """

open unit 1 card read name """ + self.rtf + """ 
read rtf card append unit 1
close unit 1"""
        else:
            charmm_inp = """* """ + charmm_title + """
*
bomlev -2

! Read in Topology and  Parameter files

open unit 1 card read name """ + self.rtf + """ 
read RTF card unit 1
close unit 1"""


        if postdata.has_key("useqmmm"):
            charmm_inp += """
read rtf card append
* QM/MM topology (for link atom)
*
   22     1
mass    98 QQH    1.00800 ! link atom

end
"""

        if not self.prm_append_replace or self.prm_append_replace == 'append':
            charmm_inp = charmm_inp + """

open unit 1 card read name """ + rtf_prm_dict['prm'] + """ 
read PARA card unit 1
close unit 1
"""
            if(self.prm):
	        charmm_inp = charmm_inp + """

open unit 1 card read name """ + self.prm + """ 
read para card append unit 1
close unit 1

""" 
        else:
            charmm_inp = charmm_inp + """

open unit 1 card read name """ + self.prm + """ 
read PARA card unit 1
close unit 1
"""

        if postdata.has_key("useqmmm"):
            charmm_inp = charmm_inp + """
read param card append
* QM/MM parameters (for link atom)
*

BONDS
qqh  cc      0.0       1.0    ! Link atom

ANGLES
qqh  cc   ct1    0.0        0.0    ! Link atom

NONBONDED
qqh    0.000000   0.000000     0.000000 ! Link atom

end
"""

	   
        try:
            postdata['solvate_implicitly']
            charmm_inp = charmm_inp + """

!This is for implicit solvation
!Scpism.inp is really a parameter file
open read unit 84 card name /usr/local/charmming/solvation/scpism.inp
"""
        except:
            pass
 
	return charmm_inp
	
    def getTopologyList(self):
        rtf_prm_dict = self.getRtfPrmPath()
        rlist = []
        if self.rtf_append_replace == None or self.rtf_append_replace == 'append':
            rlist.append(rtf_prm_dict['rtf'])
            if self.rtf:
                rlist.append(self.rtf)
        else:
            rlist.append(self.rtf)
        return rlist

    def getParameterList(self):
        rtf_prm_dict = self.getRtfPrmPath()
        rlist = []
        if not self.prm_append_replace or self.prm_append_replace == 'append':
            rlist.append(rtf_prm_dict['prm'])
            if self.prm:
                rlist.append(self.prm)
        else:
            rlist.append(self.prm)
        return rlist

    def makeBasicInput_tpl(self,seg,postdata,scriptlist):
        goodhet = re.compile('goodhet')
        line_is_goodhet = 0
        segid = seg

        #output line allows for shorthand notation so when
        #secnding charmm output write commands, output_line will be
        #a-goodhet or a
        output_line = "-" + segid
        #charmm_inp = self.makeCHARMMInputHeader('Run Segment Through CHARMM',postdata)

        # template dictionary passes the needed variables to the template
        template_dict = {}
        template_dict['topology_list'] = self.getTopologyList()
        template_dict['parameter_list'] = self.getParameterList()
        template_dict['filebase'] = self.stripDotPDB(self.filename)
        template_dict['seg'] = seg
        template_dict['suffix'] = ''

        #Custom user sequences are treated differently
        if(seg == 'sequ-pro'):
            #The sequ filename stores in the PDBInfo filename
            #differs from the actual filename due to comaptibility
            #problems otherwise so sequ_filename is the actual filename
            sequ_filename = "new_" + self.stripDotPDB(self.filename) + "-sequ-pro.pdb"
            sequ_handle = open(self.location + sequ_filename,'r')
            sequ_line = sequ_handle.read()
            sequ_line.strip()
            sequ_list = sequ_line.split(' ')
            number_of_sequences = len(sequ_list)
            template_dict['sequ_line'] = sequ_line
            template_dict['number_of_sequences'] = `number_of_sequences`  
        else:
            #if the line resembles a-tip
            #this will make it become a and set a flag
            if(goodhet.search(seg)):
                segid = self.segidFromGoodhet(seg)
                output_line = "-" + segid + "-goodhet"
                line_is_goodhet = 1
                template_dict['suffix'] = '-goodhet'

        template_dict['line_is_goodhet'] = line_is_goodhet  

        # handles patching if it exists, but we don't want to patch goodhet segments
        # since we haven't generated them yet.
        # Tim Miller: Mar 9, 2009 -- we also want to make sure we don't patch hetatom
        # segments along with regular protein statement, hence the third clause of the
        # following if...
        template_dict['larr_list'] = []
        template_dict['segpatch_name'] = self.segpatch_name 
        if(self.segpatch_name and not line_is_goodhet and not segid.endswith("-het")):
            segre = re.compile('SEGMENT ([a-z]+)')

            patch_file = open(self.location + self.segpatch_name,'r')
            for line in patch_file:
                larr = line.split('::')
                if len(larr) != 2:
                    # ToDo, really ought to throw an error here
                    continue
                # check for the correct line for this segment.
                m = segre.match(larr[0])
                if not m:
                    # ToDo, really ought to throw an error here
                    continue
                if m.group(1) == segid:
                    template_dict['larr_list'].append(larr[1])

        template_dict['endswith_het'] = segid.endswith("-het")
        template_dict['endswith_dna'] = segid.endswith("-dna")
        template_dict['endswith_rna'] = segid.endswith("-rna")
        template_dict['pproto'] = ''
        template_dict['location'] = self.location
        if not line_is_goodhet:
            # we need to handle protonation patching (if it exists)
            try:
                os.stat(self.location + "pproto_" + self.stripDotPDB(self.filename) + "-" + segid + ".str")
                template_dict['pproto'] = 'y'
            except:
                pass


        # if we have bad hetatms, use the genrtf supplied coordinate file...
        #unless they have their own replace topology parameter file
        template_dict['rtf_append_replace'] = self.rtf_append_replace
        template_dict['prm_append_replace'] = self.prm_append_replace
        template_dict['segid'] = segid
        template_dict['output_line'] = output_line

        t = get_template('%s/mytemplates/input_scripts/makeBasicInput_template.inp' % charmming_config.charmming_root)
        charmm_inp = output.tidyInp(t.render(Context(template_dict)))

        user_id = self.owner.id
        os.chdir(self.location)
        charmm_inp_filename = "charmm-"  + self.stripDotPDB(self.filename) + output_line + ".inp"
        charmm_inp_file = open(charmm_inp_filename, 'w')
        charmm_inp_file.write(charmm_inp)
        charmm_inp_file.close()
        #send to job queue
        si = schedInterface()
        #si.submitJob(user_id,self.location,charmm_inp_filename)
        scriptlist.append(charmm_inp_filename)

    def handlePatching(self,postdata):
        seg_list = self.segids.split()
        #used in disulfide bond patching
        disustartsegid = [] 
        disuendsegid = [] 
        disustart = [] 
        disuend = [] 
	patch_list = []
	patch_line = ""
        for seg in seg_list:
	    try:
	        #this part of the patch contains the CTER/NTER type patches
	        patch_list.append(postdata['first_patch'+seg])
	        patch_list.append(postdata['last_patch'+seg])
                patch_line = patch_line + "SEGMENT " + seg + " :: generate setu " + seg + " first " +\
		             postdata['first_patch'+seg] + " last " + postdata['last_patch'+seg] + "\n"
	    except:
	        pass
        if patch_line:
            self.segpatch_name = 'new_' +  self.stripDotPDB(self.filename) + '-segpatch.txt'
            self.save()
            segpatch_fp = open(self.location + self.segpatch_name,'w')
            segpatch_fp.write(patch_line)
            segpatch_fp.close()

	#This deals with the disulfide bond patching
        patch_line = ""
	try:
	    num_patches = int(postdata['num_patches'])
            for i in range(num_patches):
	        disustartsegid.append(postdata['disustartsegid' + str(i)])
	        disustart.append(postdata['disustart' + str(i)])
	        disuend.append(postdata['disuend' + str(i)])
	        disuendsegid.append(postdata['disuendsegid' + str(i)])
	        patch_line = patch_line + "patch disu " + \
                             disustartsegid[i] + " " + disustart[i] + " " +\
                             disuendsegid[i] + " " + disuend[i] + \
                             " setup\n"
        except:
   	    pass
        if patch_line:
            self.patch_name = 'new_' + self.stripDotPDB(self.filename) + '-patch.txt'
            self.save()
            patch_file = open(self.location + self.patch_name,'w')
            patch_file.write(patch_line)
            patch_file.close()

    #returns the warnings as a list if a warnings file exists
    #false otherwise
    def ifWarningsExist(self):
        try:
	    warnings = []
	    warningsfp = open(self.location + self.stripDotPDB(self.filename) + "-warnings.txt","r")
	    for line in warningsfp:
	        warnings.append(line)
	    warningsfp.close()
	    return warnings
	except:
	    return ''

    #check request data for malicious code
    def checkRequestData(self,request):
        for parameter in request.POST:
	    self.checkForMaliciousCode(request.POST[parameter],request)
        for parameter in request.GET:
	    self.checkForMaliciousCode(request.GET[parameter],request)

    #checks data for possible malicious code
    def checkForMaliciousCode(self,text,request):
        #syst is how charmm executes system commands
        syst = re.compile('syst')
	semicolon = re.compile(';')
	osdot = re.compile('os\.')
	if(syst.search(text) or semicolon.search(text) or osdot.search(text)):
	    msg = "User :" + request.user.username + "\n tried to execute malicous code with:\n " + text + "\n"
	    mail_admins('Malcious Code Warning',msg,fail_silently=False)
	    sys.exit()
	    return "Dangerous Data! Attempt has been logged."
	return text
        
  
    #gets the restraints and returns them as a string
    def handleRestraints(self,request):
        try:
            num_of_restraints = int(request.POST['num_restraints'])
	except:
	    return ""
	restraints = ""
	for i in range(num_of_restraints):
	    #con_sele = postdata['cons_selection' + `i`]
	    con_sele = self.checkForMaliciousCode(request.POST['cons_selection'+`i`],request)
	    restraints = restraints + "cons harm bestfit mass force 100.0 select " +\
	                  con_sele + ' end\n'
	return restraints

    
        
class PDBFileForm(forms.Form):
    pdbid = forms.CharField(max_length=5)   
    sequ = forms.CharField(widget=forms.widgets.Textarea())   
    pdbupload = forms.FileField()
    psfupload = forms.FileField()
    rtf_file = forms.FileField()
    prm_file = forms.FileField()

    # Go model stuffs
    gm_dm_file = forms.FileField()
    gm_dm_string = forms.CharField(max_length=25)
    gm_nScale = forms.CharField(max_length=10,initial="0.05")
    gm_domainScale = forms.CharField(max_length=10,initial="1.0")
    gm_kBond = forms.CharField(max_length=10,initial="50.0")
    gm_kAngle = forms.CharField(max_length=10,initial="30.0")

    # BLN model stuffs
    bln_dm_file = forms.FileField()
    bln_dm_string = forms.CharField(max_length=25)
    bln_nScale = forms.CharField(max_length=10,initial="1.0")
    bln_domainScale = forms.CharField(max_length=8,initial="1.0")
    bln_kBondHelix = forms.CharField(max_length=8,initial="3.5")
    bln_kBondSheet = forms.CharField(max_length=8,initial="3.5")
    bln_kBondCoil  = forms.CharField(max_length=8,initial="2.5")
    bln_kAngleHelix = forms.CharField(max_length=8,initial="8.37")
    bln_kAngleSheet = forms.CharField(max_length=8,initial="8.37")
    bln_kAngleCoil  = forms.CharField(max_length=8,initial="5.98")


class ParseException(Exception):
    reason = "No Reason"
    def __init__(self,reason):
        self.reason = reason

class energyParams(models.Model):
    pdb = models.ForeignKey(PDBFile,null=True)
    finale = models.DecimalField(max_digits=12,decimal_places=5,null=True)
    selected = models.CharField(max_length=1)
    usepbc = models.CharField(max_length=1)
    useqmmm = models.CharField(max_length=1)
    qmmmsel = models.CharField(max_length=250)

class goModel(models.Model):
    selected = models.CharField(max_length=1)
    pdb = models.ForeignKey(PDBFile)
    contactType = models.CharField(max_length=10)
    domainData  = models.CharField(max_length=20)
    nScale      = models.DecimalField(max_digits=6,decimal_places=3)
    domainScale = models.DecimalField(max_digits=6,decimal_places=3)
    kBond       = models.DecimalField(max_digits=6,decimal_places=3)
    kAngle      = models.DecimalField(max_digits=6,decimal_places=3)

class blnModel(models.Model):
    domainData  = models.CharField(max_length=20)
    nScale      = models.DecimalField(max_digits=6,decimal_places=3)
    domainScale = models.DecimalField(max_digits=6,decimal_places=3)
    kBondHelix  = models.DecimalField(max_digits=6,decimal_places=3)
    kBondCoil   = models.DecimalField(max_digits=6,decimal_places=3)
    kBondSheet  = models.DecimalField(max_digits=6,decimal_places=3)
    kAngleHelix = models.DecimalField(max_digits=6,decimal_places=3)
    kAngleCoil  = models.DecimalField(max_digits=6,decimal_places=3)
    kAngleSheet = models.DecimalField(max_digits=6,decimal_places=3)
