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
import charmming_config, input
import commands, datetime, sys, re, os, glob, smtplib
import lesson1, normalmodes, dynamics, minimization
import solvation, lessonaux, apbs
import string, output, charmming_config
import toppar.Top, toppar.Par, lib.Etc

class Structure(models.Model):

    owner = models.ForeignKey(User)
    lesson_type = models.CharField(max_length=50,null=True)
    lesson_id = models.PositiveIntegerField(default=0,null=True)

    natom  = models.PositiveIntegerField(default=0)
    name   = models.CharField(max_length=100)
    pickle = models.CharField(max_length=100)

    pdb_disul = models.CharField(max_length=100)
    location = models.CharField(max_length=200) 
    title = models.CharField(max_length=250) 
    author = models.CharField(max_length=250) 
    journal = models.CharField(max_length=250)
    pub_date = models.DateTimeField(default=datetime.datetime.now)
    selected = models.CharField(max_length=1)  
    append_status = models.CharField(max_length=1) 
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

    #Returns a list of files not specifically associated with the structure
    def getNonStructureFiles(self):
        file_list = []
	file_list.append("/solvation/water.crd")
        file_list.append("/scripts/savegv.py")
        file_list.append("/scripts/savechrg.py")
        file_list.append("/scripts/calcewald.pl")
        return file_list

    #takes the information from the Remark statement of
    #a PDB and determines the title, jrnl, and author
    def getHeader(self,pdbHeader):
        remark = re.compile('REMARK')
        title = re.compile('TITLE')
        jrnl = re.compile('JRNL')
        author = re.compile('AUTHOR') 
        ref = re.compile('REF') 
        refn = re.compile('REFN')
        for line in pdbHeader:
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
                dspatch = structure.Patch()
                dspatch.structure = self
                dspatch.patch_name = 'disul'
                try:
                    dspatch.patch_atoms = '%s:%s-%s:%s\n' % (tlist[3].lower(), tlist[4], tlist[6].lower(), tlist[7])
                    dspatch.save()
                except:
                    # ToDo, we should raise some sort of parse error here, but for now just pass
                    pass

	if self.title:
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


    ### CHARMMing > 0.9 methods go here. Eventually, all of the methods above this point ###
    ### should be brought below it and cleaned up, or eliminated.                        ###

    def setupSeg(self,seg,postdata,scriptlist):
        # template dictionary passes the needed variables to the template
        template_dict = {}
        template_dict['topology_list'] = seg.rtf_list.split(',')
        template_dict['parameter_list'] = seg.prm_list.split(',')
        template_dict['segname'] = seg.name

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

        template_dict['line_is_goodhet'] = seg.name.endswith('good')

        # handles patching if it exists, but we don't want to patch goodhet segments
        # since we haven't generated them yet.
        # Tim Miller: Mar 9, 2009 -- we also want to make sure we don't patch hetatom
        # segments along with regular protein statement, hence the third clause of the
        # following if...
        template_dict['larr_list'] = []
        template_dict['segpatch_name'] = self.segpatch_name 
        if(self.segpatch_name and not line_is_goodhet and not seg.name.endswith("-het")):
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
        template_dict['segid'] = seg.name
        template_dict['output_line'] = output_line

        t = get_template('%s/mytemplates/input_scripts/makeBasicInput_template.inp' % charmming_config.charmming_root)
        charmm_inp = output.tidyInp(t.render(Context(template_dict)))

        user_id = self.owner.id
        os.chdir(self.location)
        charmm_inp_filename = "build-"  + seg.name + ".inp"
        charmm_inp_file = open(charmm_inp_filename, 'w')
        charmm_inp_file.write(charmm_inp)
        charmm_inp_file.close()
        #send to job queue
        si = schedInterface()
        #si.submitJob(user_id,self.location,charmm_inp_filename)
        scriptlist.append(charmm_inp_filename)



    def getMoleculeFiles(self):
        """Get all molecule files associated with this structure. Note that "inherent"
           structural objects (i.e. the individual segment files) do not have StructureFile
           objects --- they're handled through the segment objects. The file name and
           description is returned as a tuple.
        """
        rlist = []

        # get list of all .psf files associated w/ this structure ... it is assumed that if the
        # psf exists then the crd will as well.
        for x in StructureFile.objects.filter(structure=self):
            if x.path.endswith('.psf'):
                rlist.append((x.path.replace(".psf",""), x.description))

        # now append the inherent files
        seglst = Segment.objects.filter(structure=self)
        for s in seglst:
            rlist.append((s.name, "Segment %s" % s.name))

        return rlist
    

    def getCHARMMFiles(self):
        """Get all CHARMM input, output, and stram files associated with this structure.
        """
        return [x for x in structure.models.StructureFile.objects.filter(structure=self) if x.endswith(".inp") or x.endswith(".out")]

class StructureFile(models.Model):
    structure   = models.ForeignKey(Structure)
    path        = models.CharField(max_length=100)
    version     = models.PositiveIntegerField(default=1)
    type        = models.CharField(max_length=20)
    description = models.CharField(max_length=500)

    # temp variable
    fd    = None

    def open(self,mode):
        self.fd = open(path,mode)

    def close(self):
        self.fd.close()

class Segment(models.Model):
    structure   = models.ForeignKey(Structure)
    is_appended = models.CharField(max_length=1)
    name        = models.CharField(max_length=6)
    type        = models.CharField(max_length=10)
    patch_first = models.CharField(max_length=100)
    patch_last  = models.CharField(max_length=100)
    rtf_list    = models.CharField(max_length=500)
    prm_list    = models.CharField(max_length=500)

class Patch(models.Model):
    structure   = models.ForeignKey(Structure)
    patch_name  = models.CharField(max_length=10)
    patch_atoms = models.CharField(max_length=100)

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
    pdb = models.ForeignKey(Structure,null=True)
    finale = models.DecimalField(max_digits=12,decimal_places=5,null=True)
    selected = models.CharField(max_length=1)
    usepbc = models.CharField(max_length=1)
    useqmmm = models.CharField(max_length=1)
    qmmmsel = models.CharField(max_length=250)

class goModel(models.Model):
    selected = models.CharField(max_length=1)
    pdb = models.ForeignKey(Structure)
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
