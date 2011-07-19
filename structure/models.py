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
import normalmodes, dynamics, minimization
import solvation, lessons.models, apbs
import string, output, charmming_config
import toppar.Top, toppar.Par, lib.Etc
import cPickle
import pychm.io, pychm.lib

class Structure(models.Model):

    owner = models.ForeignKey(User)
    selected = models.CharField(max_length=1,default='n')
    lesson_type = models.CharField(max_length=50,null=True)
    lesson_id = models.PositiveIntegerField(default=0,null=True)

    natom  = models.PositiveIntegerField(default=0)
    name   = models.CharField(max_length=100)
    pickle = models.CharField(max_length=100)

    pdb_disul = models.CharField(max_length=250)
    location = models.CharField(max_length=200) 
    title = models.CharField(max_length=250) 
    author = models.CharField(max_length=250) 
    journal = models.CharField(max_length=250)
    pub_date = models.DateTimeField(default=datetime.datetime.now)
    domains = models.CharField(max_length=250,default='')
    fes4list = models.CharField(max_length=250,default='')
    
    #Returns a list of files not specifically associated with the structure
    def getNonStructureFiles(self):
        file_list = []
	file_list.append("/solvation/water.crd")
        file_list.append("/scripts/savegv.py")
        file_list.append("/scripts/savechrg.py")
        file_list.append("/scripts/calcewald.pl")
        return file_list

    ### CHARMMing GUTS-ified methods go here. Eventually, all of the methods above this point ###
    ### should be brought below it and cleaned up, or eliminated.                             ###

    # GUTS-ified disulfide list builder
    # Note: this returns a list of tuples
    def getDisulfideList(self):
        if not self.pdb_disul:
            return None

        logfp = open('/tmp/getDisul.log', 'w')

        dsl = self.pdb_disul.split()
        logfp.write('dsl = %s\n' % dsl)
        logfp.flush()
        n = 0
        rarr = []
        logfp.write('len dsl = %d\n' % len(dsl))
        logfp.flush()
        while n < len(dsl):
            idx = dsl[n]
            resn1 = dsl[n+1]
            segn1 = dsl[n+2] + '-pro' # protein segments have disulfides
            resi1 = dsl[n+3]
            resn2 = dsl[n+4]
            segn2 = dsl[n+5] + '-pro'
            resi2 = dsl[n+6]
            n += 7
            logfp.write('got %s: %s %s %s disul to %s %s %s\n' % (idx,segn1,resn1,resi1,segn2,resn2,resi2))
            logfp.flush()
            rarr.append((segn1,resn1,resi1,segn2,resn2,resi2))

        logfp.close()
        return rarr


    #takes the information from the Remark statement of
    #a PDB and determines the title, jrnl, and author
    def getHeader(self,pdbHeader):
        for line in pdbHeader:
            if line.startswith('title'):
                line = line.replace('title','')
                self.title += line.strip()
            elif line.startswith('author'):
                line = line.replace('author','')
                self.author += line.strip()
            elif line.startswith('jrnl') or line.startswith('ref') or line.startswith('refn'):
                line = line.replace('jrnl','')
                line = line.replace('ref','')
                line = line.replace('refn','')
                self.journal += line.strip()
            elif line.startswith('ssbond'):
                # process disulfide bridges
                line = line.replace('ssbond','')
                line = ' '.join(line.split()[:7]) # grab the first seven elements of the line

                self.pdb_disul += ' %s' % line.strip()

	if self.title:
            if len(self.title) > 249:
                self.title = self.title[0:248]
	else:
	    self.title = "No information found"
	if self.author:
            if len(self.author) > 249:
                self.author = self.author[0:248]
	else:
	    self.author = "No information found"
	self.journal = self.journal.strip()
	if self.journal:
            if len(self.journal) > 249:
                self.journal = self.journal[0:248]
	else:
	    self.journal = "No information found"
        if self.pdb_disul:
            if len(self.pdb_disul) > 249:
                raise AssertionError('Too many disulfides in PDB')

        self.save()
 

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

class Segment(models.Model):
    structure   = models.ForeignKey(Structure)
    name        = models.CharField(max_length=6)
    type        = models.CharField(max_length=10)
    default_patch_first = models.CharField(max_length=100)
    default_patch_last  = models.CharField(max_length=100)
    rtf_list    = models.CharField(max_length=500)
    prm_list    = models.CharField(max_length=500)
    is_working  = models.CharField(max_length=1,default='n')

    def set_default_patches(self,firstres):
        if self.type == 'pro':
            if firstres == 'gly':
                self.default_patch_first = 'GLYP'
            else:
                self.default_patch_first = 'NTER'
            self.default_patch_last  = 'CTER'
        elif self.type == 'dna' or self.type == 'rna':
            self.default_patch_first = '5TER'
            self.default_patch_last  = '3TER'
        else:
            self.default_patch_first = 'NONE'
            self.default_patch_last = 'NONE'
        self.save()

    def getProtonizableResidues(self,model=None):
        pfp = open(self.structure.pickle, 'r')
        pdb = cPickle.load(pfp)

        if not model: 
            mol = pdb.iter_models().next() # grab first model
        else:
            mol = pdb[model]

        # ToDo: if there is more than one model, the protonizable
        # residues could change. We don't handle this right at
        # the moment, but it may not matter too much since the
        # residue sequence tends not to change between models, just
        # the atom positions.

        found = False
        for s in mol.iter_seg():
            if s.segid == self.name:
                found = True
                break

        if not found:
            raise AssertionError('Could not find right segment')

        rarr = []

        # todo, make the actual one the default
        for m in s.iter_res():
            nm = m.resName.strip()
            if nm in ['hsd','hse','hsp','his']:
                rarr.append((self.name,m.resid, \
                            [('hsd','Neutral histadine with proton on the delta carbon'), \
                             ('hse','Neutral histadine with proton on the epsilon carbon'),
                             ('hsp','Positive histadine with protons on both the delta and epsilon carbons')]))
            if nm in ['asp','aspp']:
                rarr.append((self.name,m.resid,[('asp','-1 charged aspartic acid'),('aspp','Neutral aspartic acid')]))
            if nm in ['glu','glup']:
                rarr.append((self.name,m.resid,[('glu','-1 charged glutamic acid'),('glup','Neutral glutamic acid')]))
            if nm in ['lys','lsn']:
                rarr.append((self.name,m.resid,[('lys','+1 charged Lysine'),('lsn','Neutral Lysine')]))

        pfp.close()
        return rarr

    @property
    def possible_firstpatches(self):
        if self.type == 'pro':
            return ['NTER','GLYP','PROP','ACE','ACP','NONE']
        elif self.type == 'dna' or self.type == 'rna':
            return ['5TER','5MET','5PHO','5POM','CY35']
        else:
            return ['NONE']

    @property
    def possible_lastpatches(self):
        if self.type == 'pro' :
            return ['CTER','CT1','CT2','CT3','ACP','NONE']
        elif self.type == 'dna' or self.type == 'rna':
            return ['3TER','3PHO','3POM','3CO3','CY35']
        else:
            return ['NONE']


class WorkingSegment(Segment):
    isBuilt     = models.CharField(max_length=1)
    patch_first = models.CharField(max_length=100)
    patch_last  = models.CharField(max_length=100)
    builtPSF    = models.CharField(max_length=100)
    builtCRD    = models.CharField(max_length=100)
    tpMethod    = models.CharField(max_length=20,default="standard")

    def set_terminal_patches(self,postdata):
        if postdata.has_key('first_patch' + self.name):
            self.patch_first = postdata['first_patch' + self.name]
        if postdata.has_key('last_patch' + segobj.name):
            self.patch_last = postdata['first_patch' + self.name]
        self.save()

    # This method will handle the building of the segment, including
    # any terminal patching
    def build(self,mdlname,workstruct,scriptlist):
        template_dict = {}

        # write out the PDB file in case it doesn't exist
        fp = open(self.structure.pickle, 'r')
        mol = (cPickle.load(fp))[mdlname]
        for s in mol.iter_seg():
            if s.segid == self.name:
                fname = self.structure.location + "/" + "segment-" + self.name + ".pdb"
                s.write(fname, outformat="charmm")
                template_dict['pdb_in'] = fname
                break
        fp.close()

        # see if we need to build any topology or param files
        if self.tpMethod == 'genrtf':
            self.makeGenRTF()
        elif self.tpMethod == 'antechamber':
            self.makeAntechamber()
        elif self.tpMethod == 'cgenff':
            self.makeCGenFF()

        # template dictionary passes the needed variables to the template
        template_dict['topology_list'] = self.rtf_list.split(' ')
        template_dict['parameter_list'] = self.prm_list.split(' ')
        template_dict['patch_first'] = self.patch_first
        template_dict['patch_last'] = self.patch_last
        template_dict['segname'] = self.name
        template_dict['outname'] = self.name + '-' + str(self.id)
        if self.type == 'good':
            template_dict['doic']      = False
            template_dict['noangdihe'] = 'noangle nodihedral'
        else:
            template_dict['doic']      = True
            template_dict['noangdihe'] = ''

        # Custom user sequences are treated differently
        # NB: this has not been gutsified at all...
        if(self.name == 'sequ-pro'):
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

        # Right around here is where we would want to handle any
        # protonation patching for this segment.
        patch_lines = ''
        patches = Patch.objects.filter(structure=workstruct)
        for patch in patches:
            pseg = patch.patch_segid
            if pseg and pseg.name == self.name:
                # this is a patch that we want to apply
                if patch.patch_name.startswith('hs'):
                   patch_lines += 'rename resn %s sele resid %s end\n' % (patch.patch_name,patch.patch_segres.split()[1])
                else:
                   patch_lines += 'patch %s %s\n' % (patch.patch_name,patch.patch_segres)

        if patch_lines: template_dict['patch_lines'] = patch_lines

        # write out the job script
        t = get_template('%s/mytemplates/input_scripts/buildSeg.inp' % charmming_config.charmming_root)
        charmm_inp = output.tidyInp(t.render(Context(template_dict)))
        charmm_inp_filename = self.structure.location + "/build-"  + template_dict['outname'] + ".inp"
        charmm_inp_file = open(charmm_inp_filename, 'w')
        charmm_inp_file.write(charmm_inp)
        charmm_inp_file.close()

        user_id = self.structure.owner.id
        os.chdir(self.structure.location)

        # send to job queue
        si = schedInterface()
        scriptlist.append(charmm_inp_filename)

        self.builtPSF = template_dict['outname'] + '.psf'
        self.builtCRD = template_dict['outname'] + '.crd'
        self.isBuilt = 'y'
        self.save()

    def makeCGenFF(self):
        """
        Connects to dogmans.umaryland.edu to build topology and
        parameter files using CGenFF
        """
        pass

    # Tim Miller, make hetid RTF/PRM using antechamber
    def makeAntechamber(self):
        os.putenv("AMBERHOME", charmming_config.data_home + "/amber11")
        os.chdir(self.structure.location)

        fname = "segment-" + self.name
        acbase = "antechamber-" + fname 
        acname = acbase + ".ac"

        cmd = charmming_config.data_home + "/amber11/bin/antechamber -j 5 -fi pdb -i " + \
              fname + ".pdb -fo ac -o " + acbase + ".ac"
        status = os.system(cmd)
        if status != 0:
            # debug purposes only
            raise("Antechamber screwed up!")

        try:
            # see if output file exists
            os.stat(acbase + '.ac')
        except:
            raise("Antechamber screwed up!")

        cmd = "/usr/local/charmming/antechamber/exe/charmmgen -f ac -i " + \
              acbase + ".ac -o " + acbase + " -s " + self.name
        status = os.system(cmd)
        if status != 0:
            # also for debug purposes only
            raise("Charmmgen screwed up!")
        # run through CHARMM since Antechamber changes the residue IDs
        cmd = charmming_config.charmm_exe + " < " + acbase + ".inp > " + \
              acbase + ".out"
        os.system(cmd)

        # set the newly generated RTF/PRM to be used.
        self.rtf_list = acbase + ".rtf"
        self.prm_list = acbase + ".prm"
        self.save()

    #pre:requires a list of het segids
    #This will run genRTF through the non-compliant hetatms
    #and make topology files for them
    def makeGenRTF(self):
        os.chdir(self.structure.location)

        filebase = 'segment-' + self.name
        filename = filebase + '.pdb'

        # BABEL gets called to put Hydrogen positions in for the bonds
        os.system('/usr/bin/babel -ipdb ' + filename + ' -oxyz ' + filebase + '.xyz -h')
        os.system(charmming_config.data_home + '/genrtf-v3.3 -s ' + self.name + ' -x ' + filebase + '.xyz')

        # Run the GENRTF generated inp script through CHARMMing b/c the residue names/ids have changed.
        # so we can't trust the original PDB any more.
        os.system(charmming_config.charmm_exe + ' < ' + filebase + '.inp > ' + filebase + '.out')

        #The new rtf filename will look like genrtf-filename.rtf
        self.rtf_list = "genrtf-" + self.name + "-" + str(self.id) + ".rtf"
        self.save()
        genrtf_handle = open(filebase + ".inp",'r')
        rtf_handle = open(self.structure.location + '/' + self.rtf_list, 'w')

        #Now the input file will be parsed for the rtf lines
        #The RTF files start with read rtf card append and end with End
	rtf_handle.write('* Genrtf \n*\n')
        rtf_handle.write('   28\n\n') # need the CHARMM version under which this was generated, c28 is OK.
        mass_start = False
        logfp = open('/tmp/genrtf.txt', 'w')

        for line in genrtf_handle:
            logfp.write(line)
            if "MASS" in line:
                logfp.write('got mass!\n')
                mass_start = True
                rtf_handle.write(line)
            elif "End" in line or "END" in line:
                logfp.write('got end!\n')
		rtf_handle.write(line)
		break
	    elif mass_start:
                rtf_handle.write(line)
        rtf_handle.close()

        #now for the Paramter files
        #The new prm filename will look like genrtf-filename.prm
        self.prm_list = "genrtf-" + self.name + '-' + str(self.id) + ".prm"
        self.save()
        prm_handle = open(self.structure.location + '/' + self.prm_list, 'w')

        #Now the input file will be parsed for the prm lines
        #The PRM files start with read prm card append and end with End
        prm_handle.write('* Genrtf \n*\n\n')
        continue_to_write = False
       
        for line in genrtf_handle:
            logfp.write(line)

            if line.upper().startswith("BOND"):
                logfp.write('Got BOND Line!\n')
                prm_handle.write(line)
                continue_to_write = True
            elif line.upper().startswith("END"):
                logfp.write('Got END Line!\n')
                prm_handle.write(line)
                continue_to_write = False
                break
            elif continue_to_write:
                prm_handle.write(line)
	    
        prm_handle.close()
        genrtf_handle.close()


# The idea is that the WorkingStructure class will hold structures that
# are ready to be run through minimization, dynamics, etc.
class WorkingStructure(models.Model):
    structure = models.ForeignKey(Structure)
    identifier = models.CharField(max_length=20,default='')

    selected = models.CharField(max_length=1,default='n')
    doblncharge = models.CharField(max_length=1,default='f')
    isBuilt = models.CharField(max_length=1,default='f')
    segments = models.ManyToManyField(WorkingSegment)

    modelName = models.CharField(max_length=100,default='model0')

    # final topologies and parameters (built using Frank's TOP/PAR
    # stuff).
    finalTopology = models.CharField(max_length=50,null=True)
    finalParameter = models.CharField(max_length=50,null=True)

    # job IDs
    minimization_jobID = models.PositiveIntegerField(default=0)
    solvation_jobID = models.PositiveIntegerField(default=0)
    nma_jobID = models.PositiveIntegerField(default=0)
    md_jobID = models.PositiveIntegerField(default=0)
    ld_jobID = models.PositiveIntegerField(default=0)
    redox_jobID = models.PositiveIntegerField(default=0)
    sgld_jobID = models.PositiveIntegerField(default=0)

    # lesson
    lesson = models.ForeignKey(lessons.models.Lesson,null=True)
    

    def associate(self,structref,segids,tpdict):
        for sid in segids:
            self.structure = structref
            self.save()
            segobj = Segment.objects.filter(name=sid,structure=structref)[0]

            # right now I can't think of any better way to copy all the data
            # from the segment pbject to the working segment object, so we just
            # do it by hand.
            wseg = WorkingSegment()
            wseg.is_working = 'y'
            wseg.structure = segobj.structure
            wseg.name = segobj.name
            wseg.type = segobj.type
            wseg.default_patch_first = segobj.default_patch_first
            wseg.default_patch_last = segobj.default_patch_last
            wseg.tpMethod = tpdict[sid]

            if wseg.tpMethod == 'standard':
                wseg.rtf_list = segobj.rtf_list
                wseg.prm_list = segobj.prm_list
            elif wseg.tpMethod == 'upload':
                wseg.rtf_list = self.structure.location + '/' + self.identifier + '-' + wseg.name + '.rtf'
                wseg.prm_list = self.structure.location + '/' + self.identifier + '-' + wseg.name + '.prm'
            else:
                # custom built topology/parameter files will be handled at build time
                wseg.rtf_list = '' 
                wseg.prm_list = ''

            wseg.save()

            self.segments.add(wseg)
            self.save()

    def getTopologyList(self):
        """
        Returns a list of all topology files used by all the segments
        in this WorkingStructure.
        """
        rlist = set()
        for segobj in self.segments.all():
            for rtf in segobj.rtf_list.split(' '):
                rlist.add(rtf)
        return rlist

    def getParameterList(self):
        """
        Returns a list of all parameter files used by all the segments
        in this WorkingStructure.
        """
        rlist = set()
        for segobj in self.segments.all():
            for prm in segobj.prm_list.split(' '):
                rlist.add(prm)
        return rlist

    def getAppendPatches(self):
        """
        Get a list of all the patches that need to be applied at append
        time. For now, this means disulfide patches, but it can be used
        for other types of patching such as for FeS4 clusters in the
        future.
        """
        plist = Patch.objects.filter(structure=self)
        plines = ''
        for patch in plist.all():
            if patch.patch_segid: continue # not a structure-wide patch
            plines += 'patch %s %s' % (patch.patch_name,patch.patch_segres)

        return plines

    def build(self,scriptlist):
        """
        This method replaces minimization.append_tpl() -- it is the explicit
        step that builds a new structure and appends it to the PDB object
        in charmminglib.
        """

        tdict = {}
        # step 1: check if all segments are built
        tdict['seg_list'] = []
        tdict['output_name'] = self.identifier
        tdict['blncharge'] = False # we're not handling BLN models for now
        for segobj in self.segments.all():
            if segobj.isBuilt != 't':
                segobj.build(self.modelName,self,scriptlist)
            tdict['seg_list'].append(segobj)

        tdict['topology_list'] = self.getTopologyList()
        tdict['parameter_list'] = self.getParameterList()
        tdict['patch_lines'] = self.getAppendPatches()
        t = get_template('%s/mytemplates/input_scripts/append.inp' % charmming_config.charmming_root)
        charmm_inp = output.tidyInp(t.render(Context(tdict)))

        user_id = self.structure.owner.id
        os.chdir(self.structure.location)
        charmm_inp_filename = self.structure.location + "/build-"  + self.identifier + ".inp"
        charmm_inp_file = open(charmm_inp_filename, 'w')
        charmm_inp_file.write(charmm_inp)
        charmm_inp_file.close()
        scriptlist.append(charmm_inp_filename)

        # create workingFile object for the appended structure; this has no parent
        wf = WorkingFile()
        wf.structure = self
        wf.path = self.structure.location + '/' + self.identifier + '.crd'
        wf.canonPath = wf.path
        wf.type = 'crd'
        wf.description = 'appended structure'
        wf.parentAction = 'build'
        wf.pdbkey = 'append_' + self.identifier
        wf.save()

        self.save()

        return wf

    # Updates the status of in progress operations
    def updateActionStatus(self):

        si = schedInterface()        

        pickleFile = open(self.structure.pickle, 'r+')
        pdb = cPickle.load(pickleFile)
        pickleFile.close()
        mod = False

        if self.isBuilt != 't':
            # check if the PSF and CRD files for this structure exist
            try:
                os.stat(self.structure.location + '/' + self.identifier + '.psf')
                os.stat(self.structure.location + '/' + self.identifier + '.crd')
            except:
                pass
            else:
                self.isBuilt = 't'
                mod = True
                molobj = pychm.io.pdb.get_molFromCRD(self.structure.location + '/' + self.identifier + '.crd')
                pdb['append_' + self.identifier] = molobj
                self.save()
              

        if self.minimization_jobID != 0:

            sstring = si.checkStatus(self.minimization_jobID)
	    miniparam_obj = minimization.models.minimizeParams.objects.filter(struct=self,selected='y')[0]
            miniparam_obj.statusHTML = statsDisplay(sstring,self.minimization_jobID)
	    miniparam_obj.save()

            if 'Done' in miniparam_obj.statusHTML:
                # Create a new Mol object for this
                fname = self.structure.location + '/mini-' + self.identifier + '.crd'
                mod = True

                molobj = pychm.io.pdb.get_molFromCRD(fname)
                pdb['mini_' + self.identifier] = molobj 

                wf = WorkingFile()
                wf.structure = self
                wf.path = fname
                wf.canonPath = wf.path
                wf.type = 'crd'
                wf.description = 'minimized structure'
                wf.parent = miniparam_obj.inpStruct
                wf.parentAction = 'mini'
                wf.pdbkey = 'mini_' + self.identifier
                wf.save() 

            if 'Done' in miniparam_obj.statusHTML or 'Fail' in miniparam_obj.statusHTML:
               self.minimization_jobID = 0
               if self.lesson:
                   self.lesson.onMinimizeDone(self)

        if self.solvation_jobID != 0:
            sstring = si.checkStatus(self.solvation_jobID)
	    solvparam_obj = solvation.models.solvationParams.objects.filter(structure=self,selected='y')[0]
            solvparam_obj.statusHTML = statsDisplay(sstring,self.solvation_jobID)
	    solvparam_obj.save()

            if 'Done' in solvparam_obj.statusHTML:
                # Create a new Mol object for this
                fname = self.structure.location + '/solv-' + self.identifier + '.crd'
                mod = True

                molobj = pychm.io.pdb.get_molFromCRD(fname)
                pdb['solv_' + self.identifier] = molobj

                wf = WorkingFile()
                wf.structure = self
                wf.path = fname
                wf.canonPath = wf.path
                wf.type = 'crd'
                wf.description = 'solvated structure'
                wf.parent = solvparam_obj.inpStruct
                wf.parentAction = 'solv'
                wf.pdbkey = 'solv_' + self.identifier
                wf.save()

                if solvparam_obj.salt:
                    # create working file for neutralized struct
                    wfn = WorkingFile()
                    wfn.structure = self
                    wfn.path = fname.replace('/solv-','/neut-')
                    wfn.canonPath = wfn.path
                    wfn.type = 'crd'
                    wfn.description = 'solvated and neutralized structure'
                    wfn.parent = wf
                    wfn.parentAction = 'neut'
                    wfn.pdbkey = 'neut_' + self.identifier
                    wfn.save()

                    molobj = pychm.io.pdb.get_molFromCRD(wfn.path)
                    pdb[wfn.pdbkey] = molobj

            if 'Done' in solvparam_obj.statusHTML or 'Fail' in solvparam_obj.statusHTML:
               self.solvation_jobID = 0
               if self.lesson:
                   self.lesson.onSolvationDone(self)

        if self.nma_jobID != 0:
            sstring = si.checkStatus(self.nma_jobID)
	    nmaparam_obj = normalmodes.models.nmodeParams.objects.filter(structure=self,selected='y')[0]
            nmaparam_obj.statusHTML = statsDisplay(sstring,self.nma_jobID)
	    nmaparam_obj.save()

            nma_status = statsDisplay(sstring,self.nma_jobID)

            if 'Done' in nmaparam_obj.statusHTML:
                # Create a new Mol object for this
                fname = self.structure.location + '/nmodes-' + self.identifier + '.inp'
                mod = True


                wf = WorkingFile()
                wf.structure = self
                wf.path = fname
                wf.canonPath = wf.path
                wf.type = 'inp'
                wf.description = 'normal modes'
                wf.parent = nmaparam_obj.inpStruct
                wf.parentAction = 'nma'
                wf.pdbkey = wf.parent.pdbkey
                wf.save()

            if 'Done' in nmaparam_obj.statusHTML or 'Fail' in nmaparam_obj.statusHTML:
               self.nma_jobID = 0
               if self.lesson:
                   self.lesson.onNMADone(self)


        if self.md_jobID != 0:
            sstring = si.checkStatus(self.md_jobID)
	    mdparam_obj = dynamics.models.mdParams.objects.filter(structure=self,selected='y')[0]
            mdparam_obj.statusHTML = statsDisplay(sstring,self.md_jobID)
	    mdparam_obj.save()

            if 'Done' in mdparam_obj.statusHTML:
                fname = self.structure.location + '/md-' + self.identifier + '.crd'

                # Create a new Mol object for this
                molobj = pychm.io.pdb.get_molFromCRD(fname)
                pdb['md_' + self.identifier] = molobj
                mod = True

                wf = WorkingFile()
                wf.structure = self
                wf.path = fname
                wf.canonPath = wf.path
                wf.type = 'crd'
                wf.description = 'structure after MD'
                wf.parent = mdparam_obj.inpStruct
                wf.parentAction = 'md'
                wf.pdbkey = 'md_' + self.identifier
                wf.save()

            if 'Done' in mdparam_obj.statusHTML or 'Fail' in mdparam_obj.statusHTML:
               self.md_jobID = 0
	       self.save()
               if self.lesson:
                  self.lesson.onMDDone(self)

        if self.ld_jobID != 0:
            sstring = si.checkStatus(self.ld_jobID)
	    ldparam_obj = dynamics.models.ldParams.objects.filter(structure=self,selected='y')[0]
            ldparam_obj.statusHTML = statsDisplay(sstring,self.ld_jobID)
	    ldparam_obj.save()

            if 'Done' in ldparam_obj.statusHTML:
                # Create a new Mol object for this
                fname = self.structure.location + '/ld-' + self.identifier + '.crd'
                molobj = pychm.io.pdb.get_molFromCRD(fname)
                pdb['md_' + self.identifier] = molobj
                mod = True

                wf = WorkingFile()
                wf.structure = self
                wf.path = fname
                wf.canonPath = wf.path
                wf.type = 'crd'
                wf.description = 'structure after LD'
                wf.parent = ldparam_obj.inpStruct
                wf.parentAction = 'ld'
                wf.pdbkey = 'ld_' + self.identifier
                wf.save()

            if 'Done' in ldparam_obj.statusHTML or 'Fail' in ldparam_obj.statusHTML:
               self.ld_jobID = 0
               self.save()
               if self.lesson:
                  self.lesson.onLDDone(self)

        if self.sgld_jobID != 0:
            sstring = si.checkStatus(self.sgld_jobID)
            sgldparam_obj = dynamics.models.sgldParams.objects.filter(structure=self,selected='y')[0]
            sgldparam_obj.statusHTML = statsDisplay(sstring,self.ld_jobID)
            sgldparam_obj.save() 

            if 'Done' in ldparam_obj.statusHTML:
                # Create a new Mol object for this    
                fname = self.structure.location + '/sgld-' + self.identifier + '.crd'
                molobj = pychm.io.pdb.get_molFromCRD(fname)
                pdb['md_' + self.identifier] = molobj
                mod = True

                wf = WorkingFile()
                wf.structure = self
                wf.path = fname
                wf.canonPath = wf.path
                wf.type = 'crd'
                wf.description = 'structure after SGLD'
                wf.parent = ldparam_obj.inpStruct
                wf.parentAction = 'sgld'
                wf.pdbkey = 'sgld_' + self.identifier
                wf.save()

            if 'Done' in sgldparam_obj.statusHTML or 'Fail' in sgldparam_obj.statusHTML:
               self.sgld_jobID = 0
               self.save()    
               if self.lesson: 
                  self.lesson.onSGLDDone(self)

        if self.redox_jobID != 0:
            sstring = si.checkStatus(self.redox_jobID)
            redox_obj = apbs.models.redoxParams.objects.filter(struct=self,selected='y')[0]
            redox_obj.statusHTML = statsDisplay(sstring,self.redox_jobID)
            redox_obj.save()
            # ToDo, add lesson hooks
    
        if mod:
            # blow away the old pickle file and re-create it
            # since otherwise this does not work
            os.unlink(self.structure.pickle)

            pickleFile = open(self.structure.pickle, 'w')
            cPickle.dump(pdb,pickleFile)
            pickleFile.close()

        self.save()

class Patch(models.Model):
    structure   = models.ForeignKey(WorkingStructure)

    # patches can cross multiple segments, if this field is set, it means
    # that the patch only applies to a particular segment (i.e. it is a
    # protonation patch that should be handled at segment build time rather
    # than append time.
    patch_segid  = models.ForeignKey(Segment,null=True)

    patch_name   = models.CharField(max_length=10)
    patch_segres = models.CharField(max_length=100)

class WorkingFile(models.Model,file):
    structure   = models.ForeignKey(WorkingStructure)
    pdbkey      = models.CharField(max_length=100,null=True) # if this is a structure file, want to know where to find it
    path        = models.CharField(max_length=100)
    canonPath   = models.CharField(max_length=100) # added so we can do file versioning
    version     = models.PositiveIntegerField(default=1)
    type        = models.CharField(max_length=20)
    description = models.CharField(max_length=500,null=True)
    parent      = models.ForeignKey('self',null=True)

    # We should probably figure out a better way to do this; the problem
    # is that I can't have a foreign key that points to multiple types, so
    # I've implemented the "action" hack.
    parentAction = models.CharField(max_length=6)

    @property
    def basename(self):
        bn = self.path.split('/')[-1]
        return bn.split('.')[0]
       

    @property
    def cbasename(self):
        bn = self.canonPath.split('/')[-1]
        return bn.split('.')[0]

    def backup(self):
        # eventually, this will allow us to create new
        # versions of files
        pass

# --- below this point are the classes for the various forms ---

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
    struct = models.ForeignKey(WorkingStructure,null=True)
    inpStruct = models.ForeignKey(WorkingFile,null=True)

    finale = models.FloatField(null=True)
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
