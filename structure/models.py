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
from django.http import HttpResponseRedirect, HttpResponse
from django.template.loader import get_template
from django.db import models
from django import forms
from django.contrib.auth.models import User
from django.core.mail import mail_admins
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from normalmodes.aux import parseNormalModes, getNormalModeMovieNum
from pychm.cg.sansombln import SansomBLN
from httplib import HTTPConnection
import charmming_config, input
import commands, datetime, sys, re, os, glob, smtplib
import normalmodes, dynamics, minimization
import solvation, lessons.models, apbs, lesson1, lesson2, lesson3, lesson4, lesson, lessonaux
import string, output, charmming_config
import toppar.Top, toppar.Par, lib.Etc
import cPickle, copy, traceback, socket
import pychm.io, pychm.lib, pychm.cg
import pychm.future.lib.toppar as pychm_toppar
from pychm.future.io.charmm import open_rtf, open_prm

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

        # we need to translate between old and new numbers how
        pfp = open(self.pickle, 'r')
        pdb = cPickle.load(pfp)
        pfp.close()
        mdl = pdb[0]

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

            # make sure we have the correct residue number -- hack to just
            # use model 0
            for segment in mdl.iter_seg():
                logfp.write('Checking segment %s -- segn1 = %s segn2 = %s\n' % (segment.segid,segn1,segn2))

                if segment.segid == segn1:
                    for atom in segment:
                        if atom.resid0 == int(resi1):
                            resi1 = str(atom.resid)
                            logfp.write("Renumbering: changed resid to %s\n" % resi1)

                if segment.segid == segn2:
                    for atom in segment:
                        if atom.resid0 == int(resi2):
                            resi2 = str(atom.resid)
                            logfp.write("Renumbering: changed resid to %s\n" % resi2)

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


    def putSeqOnDisk(self, sequence, path):
        """Special method used with custom sequences; writes out a PDB from a sequence.
        """

        # ToDo: check and make sure that the residues in the sequence are valid
        sequence = sequence.strip()
        x = sequence.split()
        seqrdup  = ' '.join(x)

        if len(x) > 50:
            raise AssertionError('Custom sequence is to long!')

        td = {}
        td['nres'] = len(x)
        td['sequence'] = seqrdup
        td['topology'] = '%s/toppar/top_all27_prot_na.rtf' % charmming_config.data_home
        td['parameter'] = '%s/toppar/par_all27_prot_na.prm' % charmming_config.data_home
        td['outname'] = path
        td['name'] = self.owner.username
        
        # this needs to be run directly through CHARMM
        t = get_template('%s/mytemplates/input_scripts/seqtopdb.inp' % charmming_config.charmming_root)
        charmm_inp = output.tidyInp(t.render(Context(td)))

        os.chdir(self.location)
        fp = open('seq2pdb.inp', 'w')
        fp.write(charmm_inp)
        fp.close()

        logfp = open('/tmp/seqtopdb.txt', 'w')
        logfp.write('Running %s < seq2pdb.inp > seq2pdb.out\n' % charmming_config.charmm_exe)
        os.system("LD_LIBRARY_PATH=%s %s < seq2pdb.inp > seq2pdb.out" % (charmming_config.lib_path,charmming_config.charmm_exe))
        logfp.write('Done\n')
        logfp.close()

class Segment(models.Model):
    structure   = models.ForeignKey(Structure)
    name        = models.CharField(max_length=6)
    type        = models.CharField(max_length=10)
    default_patch_first = models.CharField(max_length=100)
    default_patch_last  = models.CharField(max_length=100)
    rtf_list    = models.CharField(max_length=500)
    prm_list    = models.CharField(max_length=500)
    is_working  = models.CharField(max_length=1,default='n')
    fes4        = models.BooleanField(default=False) # a bit of a hack, but it makes the views code easier

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
        """
        Returns tuples of information. ToDo: Document this!!!
        """

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
    redox       = models.BooleanField(default=False)

    def set_terminal_patches(self,postdata):
        if postdata.has_key('first_patch' + self.name):
            self.patch_first = postdata['first_patch' + self.name]
        if postdata.has_key('last_patch' + segobj.name):
            self.patch_last = postdata['first_patch' + self.name]
        self.save()

    def getUniqueResidues(self,mol):
        """
        This routine splits up the badhet segment into each of its unique residues, so
        they can be run through CGenFF/MATCH one at a time. It also write out a MOL2 file
        with added hydrogens for each residue.
        """
        logfp = open('/tmp/genUniqueResidues.txt', 'w')
        logfp.write('my name is %s\n' % self.name)

        badResList = []
        found = False
        for seg in mol.iter_seg():
            logfp.write('comp with %s\n' % seg.segid)
            if seg.segid == self.name:
                found = True
                break
        logfp.close()
        if not found:
            raise AssertionError('Asked to operate on a nonexistent segment!')

        for residue in seg.iter_res():
            if residue.resName not in badResList:
                badResList.append(residue.resName)
                filename_noh = self.structure.location + '/' + self.name + '-badres-' + residue.resName + ".pdb"
                filename_sdf = self.structure.location + '/' + self.name + '-badres-h-' + residue.resName + ".sdf"
                filename_h = self.structure.location + '/' + self.name + '-badres-h-' + residue.resName + ".mol2"
                residue.write(filename_noh, outformat='pdborg')

                # shiv to try to get an SDF file for the residue
                conn = HTTPConnection("www.pdb.org")
                logfp = open('/tmp/sdf.txt', 'w')
                reqstring = "/pdb/files/ligand/%s_ideal.sdf" % residue.resName.upper()
                logfp.write("reqstring = %s\n" % reqstring)
                conn.request("GET", reqstring)
                resp = conn.getresponse()
                if resp.status == 200:
                    logfp.write('OK\n')
                    sdf_file = resp.read()

                    outfp = open(filename_sdf, 'w')
                    outfp.write(sdf_file)
                    outfp.close()

                    pdb_rewrite = self.structure.location + '/segment-' + seg.segid + '.pdb'
                    os.system("babel --title %s -isdf %s -omol2 %s" % (residue.resName,filename_sdf,filename_h))
                    ##os.system("babel --title %s -isdf %s -opdb %s" % (residue.resName,filename_sdf,pdb_rewrite))
                else:
                    logfp.write('No: use PDB\n')
                    os.system("babel -h --title %s -ipdb %s -omol2 %s" % (residue.resName,filename_noh,filename_h))
                logfp.close()
          
        return badResList

    # This method will handle the building of the segment, including
    # any terminal patching
    def build(self,mdlname,workstruct):
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

        # see if we need to build any topology or param files
        if self.tpMethod == 'autogen':
            logfp = open('/tmp/autogen.txt', 'w')
            bhResList = self.getUniqueResidues(mol)

            success = False
            for tpmeth in charmming_config.toppar_generators.split(','):
                logfp.write('trying method: %s\n' % tpmeth)
                if tpmeth == 'genrtf':
                    rval = self.makeGenRTF()
                elif tpmeth == 'antechamber':
                    rval = self.makeAntechamber()
                elif tpmeth == 'cgenff':
                    rval = self.makeCGenFF(bhResList)
                elif tpmeth == 'match':
                    rval = self.makeMatch(bhResList)

                logfp.write('got rval = %d\n' % rval)
                if rval == 0:
                    success = True
                    break
     
            logfp.close()
            if not success:
               raise AssertionError('Unable to build topology/parameters')

        elif self.tpMethod == 'dogmans':
            rval = self.makeCGenFF(bhResList)
        elif self.tpMethod == 'match':
            rval = self.makeMatch(bhResList)
        elif self.tpMethod == 'antechamber':
            rval = self.makeAntechamber()
        elif self.tpMethod == 'genrtf':
            rval = self.makeGenRTF()
        

        # done top/par file building
        fp.close()

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

        self.builtPSF = template_dict['outname'] + '.psf'
        self.builtCRD = template_dict['outname'] + '.crd'
        self.isBuilt = 'y'
        self.save()

        return charmm_inp_filename

    def makeCGenFF(self, badResList):
        """
        Connects to dogmans.umaryland.edu to build topology and
        parameter files using CGenFF
        """
        logfp = open('/tmp/makecgenff.txt', 'w')
        logfp.write('Try to use CGenFF\n')
        self.rtf_list = '%s/toppar/top_all36_cgenff.rtf' % charmming_config.data_home
        self.prm_list = '%s/toppar/par_all36_cgenff.prm' % charmming_config.data_home

        header = '# USER_IP 165.112.184.52 USER_LOGIN %s\n\n' % self.structure.owner.username
        for badRes in badResList:
            fp = open('%s/%s-badres-h-%s.mol2' % (self.structure.location,self.name,badRes), 'r')
            payload = fp.read()
            content = header + payload
            fp.close()

            recvbuff = ''
            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.connect((charmming_config.cgenff_host, charmming_config.cgenff_port))
                s.sendall(content)
                s.shutdown(socket.SHUT_WR)
                while True:
                    data = s.recv(1024)
                    if not data:
                        break
                    recvbuff += data

            except Exception, e:
                logfp.write('crap something went wrong::\n')
                traceback.print_exc(file=logfp)
                logfp.close()
                return -1

            fp = open('%s/%s-%s-dogmans.txt' % (self.structure.location,self.name,badRes), 'w') 
            fp.write(recvbuff)
            fp.close()
            logfp.write('Dumped CGenFF result for analysis.\n')

            # parse file returned from dogmans, and make sure that it has both topology and
            # parameters.
            rtffound = False
            prmfound = False
            inrtf = False
            inprm = False
            rtffp = open('%s/%s-%s-dogmans.rtf' % (self.structure.location,self.name,badRes), 'w')
            prmfp = open('%s/%s-%s-dogmans.prm' % (self.structure.location,self.name,badRes), 'w')
            for line in recvbuff.split('\n'):
                if inrtf:
                    rtffp.write(line + '\n')
                if inprm:
                    prmfp.write(line + '\n')

                line = line.strip()
                if line == 'end' or line == 'END':
                    inrtf = False
                    inprm = False
                if line.startswith('read rtf'):
                    rtffound = True
                    inrtf = True
                if line.startswith('read param'):
                    prmfound = True
                    inprm = True
            
            if not (rtffound and prmfound):
                logfp.write('Aw crap ... did not find both topology and parameters.\n')
                logfp.close()
                return -1

            rtffp.close()
            prmfp.close()
       
            # we need to adjust the naming conventions in the PDB            
            with open_rtf('%s/%s-%s-dogmans.rtf' % (self.structure.location,self.name,badRes)) as rtffp:
                new_tp = pychm_toppar.Toppar()
                rtffp.export_to_toppar(new_tp)

                name_list = []
                if new_tp.residue is not None:
                    name_list += [ resi.name for resi in new_tp.residue ]
                if new_tp.patch is not None:
                    name_list += [ resi.name for resi in new_tp.patch ]

                pdbname = self.structure.location + '/segment-' + self.name + '.pdb'
                molobj = pychm.io.pdb.get_molFromPDB(pdbname)
                for res in molobj.iter_res():
                    if res.resName in name_list:
                         res._dogmans_rename()

                molobj.write(pdbname)


            self.rtf_list += ' %s-%s-dogmans.rtf' % (self.name,badRes)
            self.prm_list += ' %s-%s-dogmans.prm' % (self.name,badRes)
        
        
        logfp.write('All looks OK\n')
        logfp.close()
        return 0



    def makeMatch(self, badResList):
        """
        Uses the match program from the Charlie Brooks group to
        try to build topology and parameter files.
        """

        logfp = open('/tmp/match.txt', 'w')

        os.putenv("PerlChemistry","%s/MATCH_RELEASE/PerlChemistry" % charmming_config.data_home)
        os.putenv("MATCH","%s/MATCH_RELEASE/MATCH" % charmming_config.data_home)
        os.chdir(self.structure.location)
       
        self.rtf_list = '%s/toppar/top_all27_prot_na.rtf' % charmming_config.data_home
        self.prm_list = '%s/toppar/par_all27_prot_na.prm' % charmming_config.data_home
        for badRes in badResList:
            exe_line = '%s/MATCH_RELEASE/MATCH/scripts/MATCH.pl %s-badres-h-%s.mol2' % (charmming_config.data_home,self.name,badRes)
            logfp.write('Xcute: %s\n' % exe_line)
            status, output = commands.getstatusoutput(exe_line)
            if status != 0:
                logfp.write("sorry ... nonzero status :-(\n")
                logfp.close()
                return -1

            logfp.write("OK!\n")

            # ToDo: add the rtf/param files that have been generated to a master topology
            # and parameter files for the entire segment.
            self.rtf_list += ' %s-badres-h-%s.rtf' % (self.name,badRes)
            self.prm_list += ' %s-badres-h-%s.prm' % (self.name,badRes)

        logfp.close()
        return 0

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
            #raise("Antechamber screwed up!")
            return -1

        try:
            # see if output file exists
            os.stat(acbase + '.ac')
        except:
            return -1

        cmd = "/usr/local/newcharmming/antechamber/exe/charmmgen -f ac -i " + \
              acbase + ".ac -o " + acbase + " -s " + self.name
        status = os.system(cmd)
        if status != 0:
            # also for debug purposes only
            #raise("Charmmgen screwed up!")
            return -1
        # run through CHARMM since Antechamber changes the residue IDs
        cmd = charmming_config.charmm_exe + " < " + acbase + ".inp > " + \
              acbase + ".out"
        os.system(cmd)

        # set the newly generated RTF/PRM to be used.
        self.rtf_list = acbase + ".rtf"
        self.prm_list = acbase + ".prm"
        self.save()
        return 0

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
        return 0


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

    qmRegion = models.CharField(max_length=250,default='none')

    # final topologies and parameters (built using Frank's TOP/PAR
    # stuff).
    finalTopology = models.CharField(max_length=50,null=True)
    finalParameter = models.CharField(max_length=50,null=True)

    # lesson
    lesson = models.ForeignKey(lessons.models.Lesson,null=True)

    # extra setup stream files that are needed (useful for BLN model
    # and possibly elsewhere)
    extraStreams = models.CharField(max_length=120,null=True)

    @property
    def dimension(self):
        xmax = -9999.
        ymax = -9999.
        zmax = -9999.
        xmin = 9999.
        ymin = 9999.
        zmin = 9999.

        # get the names of all segments that are part of this working structure
        segnamelist = []
        for wseg in self.segments.all():
            segnamelist.append(wseg.name)

        pickle = open(self.structure.pickle,'r')
        pdb = cPickle.load(pickle)
        pickle.close()

        mol = pdb[0] # potentially dangerous: assume we're dealing with model 0
        for seg in mol.iter_seg():
            if seg.segid in segnamelist:
                for atom in seg:
                    x, y, z = atom.cart
                    if x > xmax: xmax = x
                    if y > xmax: ymax = y
                    if z > xmax: zmax = z
                    if x < xmin: xmin = x
                    if y < xmin: ymin = y
                    if z < xmin: zmin = z

        logfp = open('/tmp/dimensions.txt', 'w')
        logfp.write('x = %10.6f %10.6f y = %10.6f %10.6f z = %10.6f %10.6f\n' % (xmin,xmax,ymin,ymax,zmin,zmax))
        logfp.close()

        return((xmax-xmin,ymax-ymin,zmax-zmin))

    def associate(self,structref,segids,tpdict):
        """
        use @transaction.commit_manually
        and see Django docs
        """

        for sid in segids:
            self.structure = structref
            self.save()
            segobj = Segment.objects.filter(name=sid,structure=structref)[0]

            # right now I can't think of any better way to copy all the data
            # from the segment object to the working segment object, so we just
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
                logfp = open('/tmp/assoc.txt', 'w')
                logfp.write('In assoc, identifier is %s.\n' % self.identifier)
                logfp.close()
                wseg.rtf_list = self.structure.location + '/' + self.identifier + '-' + wseg.name + '.rtf'
                wseg.prm_list = self.structure.location + '/' + self.identifier + '-' + wseg.name + '.prm'
            elif wseg.tpMethod == 'redox':
                # structure will be used ONLY for oxi/reduce calculations so doesn't need
                # top/par (redox script provides these when needed).
                wseg.rtf_list = ''
                wseg.prm_list = ''
                wseg.redox = True
            else:
                # custom built topology/parameter files will be handled at build time
                wseg.rtf_list = '' 
                wseg.prm_list = ''

            wseg.save()

            self.segments.add(wseg)
            self.save()

    def getTopparList(self):
        """
        Returns a list of all parameter files used by all the segments
        in this WorkingStructure.

        BTM 20120522 -- build one parameter file to rule them all, using
        Frank's code. 
        """

        qrebuild=False

        prm_list = []
        seglist = self.segments.all()

        bhlist = [x for x in seglist if x.name.endswith('-bad')]

        # This rather obnoxious chunk of code checks if any of the bad hets are for
        # REDOX, and if so removes them from the list, since they're not really bad
        # in the sense that they require new topology and parameter files.
        oklist = []
        for i, seg in enumerate(bhlist):
            if seg.redox: 
                oklist.append(i)
        oklist.sort(reverse=True)
        for elt in oklist:
            bhlist.pop(elt)

        if not bhlist:

            # no bad hets means no custom topology or parameters, behave as normal
            tlist = set()
            plist = set()
            for segobj in self.segments.all():
                if segobj.redox: continue
                for prm in segobj.prm_list.split(' '):
                   plist.add(prm)
                for rtf in segobj.rtf_list.split(' '):
                   tlist.add(rtf)

            return tlist, plist, False
        else:
            # Shizzle my izzle we have some custom generated topologies and parameters
            # here. Rather than worry about file ordering and inter-dependencies, we
            # will just combine them with Frank's nifty code.
            qrebuild = True # we might need to rebuild the PSF

            pnlist = [x for x in seglist if x.name.endswith('-pro') or x.name.endswith('-dna') or x.name.endswith('rna')]
            ghlist = [x for x in seglist if x.name.endswith('-good')]

            # finlist is built to put the segments in the correct order for parameter
            # file merging
            finlist = pnlist + ghlist + bhlist
            final_toppar = pychm_toppar.Toppar()

            logfp = open('/tmp/debug-toppar-combo.txt', 'w')
            for segobj in finlist:
                logfp.write('Processing segment %s\n' % segobj.name)  

                # we are going to assume that RTFs and PRMs are matched pairs in their
                # respective lists.
                rtflist = segobj.rtf_list.split(' ')
                prmlist = segobj.prm_list.split(' ')
                if len(rtflist) != len(prmlist):
                    logfp.write('Problem: rtf list length %d prm list length %d\n' % (len(rtflist),len(prmlist)))
                    raise Exception("Strange list lengths")

                for i, rtf in enumerate(rtflist):
                    logfp.write('Processing rtf %s and prm %s\n' % (rtf,prmlist[i]))

                    with open_rtf(rtf) as tmp_rtf:
                        tmp_rtf.export_to_toppar(final_toppar)
                    with open_prm(prmlist[i]) as tmp_prm:
                        tmp_prm.export_to_toppar(final_toppar)
#                    rtfobj = open_rtf(rtf)
#                    prmobj = open_prm(prmlist[i])
#                    temp_toppar = pychm_toppar.Toppar()
#                    rtfobj.export_to_toppar(temp_toppar)
#                    prmobj.export_to_toppar(temp_toppar)
#
#                    final_toppar = final_toppar + temp_toppar

            logfp.write('Done\n')
            logfp.close()

            # we have now combined all of the parameter file objects, write out a new file
            newRTFFile = self.structure.location + '/' + self.identifier + '-constructed.rtf'
            newPrmFile = self.structure.location + '/' + self.identifier + '-constructed.prm'
            tfinal = open_rtf(newRTFFile, 'w')
            pfinal = open_prm(newPrmFile, 'w')
            tfinal.import_from_toppar(final_toppar)
            pfinal.import_from_toppar(final_toppar)
            tfinal.write_all()
            pfinal.write_all()
            tfinal.close()
            pfinal.close()

            tlist = set()
            plist = set()
            tlist.add(newRTFFile)
            plist.add(newPrmFile)
            return tlist, plist, True

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
            plines += 'patch %s %s\n' % (patch.patch_name,patch.patch_segres)

        return plines

    def build(self,inTask):
        """
        This method replaces minimization.append_tpl() -- it is the explicit
        step that builds a new structure and appends it to the PDB object
        in charmminglib.
        """

        tdict = {}
        # step 1: check if all segments are built
        tdict['seg_list'] = []
        tdict['output_name'] = self.identifier + '-build'
        tdict['blncharge'] = False # we're not handling BLN models for now

        logfp = open('/tmp/test.txt', 'w')
        logfp.write('My id = %d\n' % self.id)

        # check if we're a CG model.
        try:
            cgws = CGWorkingStructure.objects.get(workingstructure_ptr=self.id)
            logfp.write("Got a CG working struct\n")
        except:
            logfp.write("No CG working struct here\n")
            cgws = None
        logfp.close()

        if cgws:
            tdict['input_pdb'] = cgws.pdbname
            tdict['finalname'] = cgws.cg_type
        else:
            # We are not a CG structure... go ahead and build up the seg_list
            # as normal.
            for segobj in self.segments.all():
                # segment to be used solely for redox calculations get
                # handled differently
                if segobj.isBuilt != 't' and not segobj.redox:
                    newScript = segobj.build(self.modelName,self)
                    if inTask.scripts:
                        inTask.scripts += ',' + newScript
                    else:
                        inTask.scripts = newScript
                if not segobj.redox:
                    tdict['seg_list'].append(segobj)

        tdict['topology_list'], tdict['parameter_list'], qrebuild = self.getTopparList()
        tdict['patch_lines'] = self.getAppendPatches()

        logfp = open('/tmp/buildLog.txt', 'w')
        logfp.write('Having called getTopparList: qrebuild = %s\n' % qrebuild)

        if qrebuild:
            tdict['rebuild'] = True
            segtuplist = []
            for segobj in self.segments.all():
                if segobj.type == 'good':
                    special = 'noangle nodihedral'
                else:
                    special = ''
                segtuplist.append((segobj.name,segobj.builtCRD.replace('.crd','.pdb'),segobj.patch_first,segobj.patch_last,special))
            tdict['segbuild'] = segtuplist
        else:
            tdict['rebuild'] = False
        logfp.write('Setting template dict rebuild = %s\n' % tdict['rebuild'])
        logfp.close() 

        t = get_template('%s/mytemplates/input_scripts/append.inp' % charmming_config.charmming_root)
        charmm_inp = output.tidyInp(t.render(Context(tdict)))

        os.chdir(self.structure.location)
        charmm_inp_filename = self.structure.location + "/"  + self.identifier + "-build.inp"
        charmm_inp_file = open(charmm_inp_filename, 'w')
        charmm_inp_file.write(charmm_inp)
        charmm_inp_file.close()
        inTask.scripts += ',%s' % charmm_inp_filename
        inTask.save()

        # create a Task object for appending; this has no parent
        task = Task()
        task.setup(self)
        task.parent = None
        task.action = 'build'
        task.active = 'y'
        task.save()

        self.save()

        return task

    def addCRDToPickle(self,fname,fkey):
        pickleFile = open(self.structure.pickle, 'r+')
        pdb = cPickle.load(pickleFile)
        pickleFile.close()

        # Create a new Mol object for this
        molobj = pychm.io.pdb.get_molFromCRD(fname)
        pdb[fkey] = molobj

        # kill off the old pickle and re-create
        os.unlink(self.structure.pickle)

        pickleFile = open(self.structure.pickle, 'w')
        cPickle.dump(pdb,pickleFile)
        pickleFile.close()
 
        self.save() 

    # Updates the status of in progress operations
    def updateActionStatus(self):
        # appending is a special case, since it doesn't exist as a task unto
        # itself. So if the structure is not built, we should check and see
        # whether or not that happened.
        if self.isBuilt != 't':
            # check if the PSF and CRD files for this structure exist
            try:
                os.stat(self.structure.location + '/' + self.identifier + '-build.psf')
                os.stat(self.structure.location + '/' + self.identifier + '-build.crd')
            except:
                pass
            else:
                self.isBuilt = 't'
                self.addCRDToPickle(self.structure.location + '/' + self.identifier + '-build.crd', 'append_' + self.identifier)
                loc = self.structure.location
                bnm = self.identifier

                # add the working files to the quasi appending Task
                buildtask = Task.objects.get(workstruct=self,finished='n',action='build')

                wfinp = WorkingFile()
                wfinp.task = buildtask
                wfinp.path = loc + '/' + bnm + '-build.inp'
                wfinp.canonPath = wfinp.path
                wfinp.type = 'inp'
                wfinp.save()

                wfout = WorkingFile()
                wfout.task = buildtask
                wfout.path = loc + '/' + bnm + '-build.out'
                wfout.canonPath = wfout.path
                wfout.type = 'out'
                wfout.save()

                wfpsf = WorkingFile()
                wfpsf.task = buildtask
                wfpsf.path = loc + '/' + bnm + '-build.psf'
                wfpsf.canonPath = wfpsf.path
                wfpsf.type = 'psf'
                wfpsf.save()

                wfpdb = WorkingFile()
                wfpdb.task = buildtask
                wfpdb.path = loc + '/' + bnm + '-build.pdb'
                wfpdb.canonPath = wfpdb.path
                wfpdb.type = 'pdb'
                wfpdb.save()

                wfcrd = WorkingFile()
                wfcrd.task = buildtask
                wfcrd.path = loc + '/' + bnm + '-build.crd'
                wfcrd.canonPath = wfcrd.path
                wfcrd.type = 'crd'
                wfcrd.pdbkey = 'append_' + self.identifier
                wfcrd.save()

                buildtask.status = 'C'
                buildtask.save()

              
        tasks = Task.objects.filter(workstruct=self,active='y',finished='n')
        #YP lesson stuff 
        try:
            lnum=self.structure.lesson_type
            lesson_obj = eval(lnum+'.models.'+lnum.capitalize()+'()')
        except:
            lesson_obj=None
        #YP

        for t in tasks:
            t.query()

            if t.status == 'C' or t.status == 'F':
                if t.action == 'minimization':
                    t2 = minimization.models.minimizeTask.objects.get(id=t.id)
                    if lesson_obj: lessonaux.doLessonAct(self.structure,"onMinimizeDone",t)
                elif t.action == 'solvation' or t.action == 'neutralization':
                    t2 = solvation.models.solvationTask.objects.get(id=t.id)
                    if lesson_obj:lessonaux.doLessonAct(self.structure,"onSolvationDone",t)
                elif t.action == 'md':
                    if lesson_obj:lessonaux.doLessonAct(self.structure,"onMDDone",t)
                    t2 = dynamics.models.mdTask.objects.get(id=t.id)
                elif t.action == 'ld':
                    t2 = dynamics.models.ldTask.objects.get(id=t.id)
                elif t.action == 'sgld':
                    if lesson_obj:lessonaux.doLessonAct(self.structure,"onSGLDDone",t)
                    t2 = dynamics.models.sgldTask.objects.get(id=t.id)
                elif t.action == 'energy':
                    t2 = energyTask.objects.get(id=t.id)
                    if lesson_obj:lessonaux.doLessonAct(self.structure,"onEnergyDone",t2)
                elif t.action == 'nmode':
                    t2 = normalmodes.models.nmodeTask.objects.get(id=t.id)
                elif t.action == 'redox':
                    t2 = apbs.models.redoxTask.objects.get(id=t.id)
                    if lesson_obj:lessonaux.doLessonAct(self.structure,"onRedoxDone",t)
                else:
                    t2 = t

                t2.finish()
                t2.finished = 'y'
                t2.save()

                #YP lessons status update                
                #try:
                #    lnum=self.structure.lesson_type
                #    lesson_num_obj = eval(lnum+'.models.'+lnum.capitalize()+'()')
                #    lesson_num_class=lesson_num_obj.__class__
                #    lesson_obj=lesson_num_class.objects.filter(id=self.structure.lesson_id)[0]
                #    lesson_obj.curStep=float(float(lesson_obj.curStep)+float(0.5))
                #    lesson_obj.save()
                #except:
                #    pass
                #YP

        self.save()

class CGWorkingStructure(WorkingStructure):
    """
    This is a special WorkingStructure type that is designed for coarse
    grained models. The main difference is how the associate routine
    works (since CG models are simple, we don't worry about patching).
    """
    pdbname = models.CharField(max_length=120) # holds the CG-ified PDB
    cg_type = models.CharField(max_length=10)

    def associate(self,structref,segids,**kwargs):
        """
        This method builds the CG model from the AA coordinates. Note
        that the kwargs is designed to hold CG-model specific parameters,
        e.g. nscale etc. These get passed directly to the pychm
        constructor of the CG model.
        """

        self.structure = structref
        self.save()

        for sid in segids:
            segobj = Segment.objects.filter(name=sid,structure=structref)[0]

            # here is where there are some changes from the parent associate
            # routine ... WorkingSegments need to have their type set to
            # modelType and no patch_first or patch_last.
            wseg = WorkingSegment()
            wseg.is_working = 'y'
            wseg.structure = segobj.structure
            wseg.name = segobj.name.replace('-pro','-%s' % self.cg_type)
            wseg.type = self.cg_type
            wseg.default_patch_first = 'None'
            wseg.default_patch_last = 'None'
            wseg.tpMethod = 'standard'

            wseg.save()

            self.segments.add(wseg)
            self.save()

            # We could actually go ahead and call wseg.build() here, but I
            # don't think that it is necessary. XXX Right now the segment
            # holds the RTF and PRM lists, which aren't set right for a CG
            # model segment. Fortunately, we can just build in a per cgWorkStruct,
            # and not a per segment basis, but this will need to be re-done
            # at some point.

        # because of the way the Go model works, we need to actually build
        # the AA structure first, then build the CG model.
        if self.cg_type == 'go' or self.cg_type == 'bln':

            logfp = open('/tmp/cgassoc.txt', 'w')
            logfp.write("In CG associate: segids = %s\n" % segids)
            logfp.write("In CG associate: cg_type = %s\n" % self.cg_type)

            pdbfname = structref.location + '/cgbase-' + str(self.id) + '.pdb'

            fp = open(structref.pickle, 'r')
            mol = (cPickle.load(fp))[self.modelName]
            sct = 0
            for s in mol.iter_seg():
                logfp.write("Checking %s ... " % s.segid)
                if s.segid in segids:
                    logfp.write("OK!\n")
                    sct += 1
                    s.write(pdbfname, outformat='charmm',append=True,ter=True,end=False)
                else:
                    logfp.write("No.\n")
            fp.close()

            if sct == 0:
                return HttpResponse("No segments selected")

            fp = open(pdbfname, 'a')
            fp.write('\nEND\n')
            fp.close()

            # AA PDB file is built, now turn it into a CG model
            pdbfile = pychm.io.pdb.PDBFile(pdbfname)[0]

            # decide how this model is going to be stored on disk.
            basefname = self.identifier + '-' + self.cg_type

            logfp.write("Point OSCAR self.cg_type = '%s'\n" % self.cg_type)
            if self.cg_type == 'go':
                # need to set strideBin via kwargs
                cgm = pychm.cg.ktgo.KTGo(pdbfile, strideBin=charmming_config.stride_bin)
                if kwargs.has_key('nScale'):
                    cgm.nScale = float(kwargs['nScale'])
                if kwargs.has_key('contactSet'):
                    cgm.contactSet = kwargs['contactSet']
                if kwargs.has_key('kBond'):
                    cgm.kBond = float(kwargs['kBond'])
                if kwargs.has_key('kAngle'):
                    cgm.kAngle = float(kwargs['kAngle'])
                if kwargs.has_key('contactrad'):
                    cgm.contactrad = float(kwargs['contactrad'])

            elif self.cg_type == 'bln':
                # kwargs can be hbondstream
                logfp.write('Do BLN code.\n')
                cgm = SansomBLN(pdbfile, **kwargs)
                if kwargs.has_key('kBondHelix'):
                    cgm.kBondHelix = float(kwargs['kBondHelix'])
                if kwargs.has_key('kBondSheet'):
                    cgm.kBondHelix = float(kwargs['kBondSheet'])
                if kwargs.has_key('kBondCoil'):
                    cgm.kBondHelix = float(kwargs['kBondCoil'])
                if kwargs.has_key('kAngleHelix'):
                    cgm.kAngleHelix = float(kwargs['kAngleHelix'])
                if kwargs.has_key('kAngleSheet'):
                    cgm.kAngleSheet = float(kwargs['kAngleSheet'])
                if kwargs.has_key('kAngleCoil'):
                    cgm.kAngleCoil = float(kwargs['kAngleCoil'])

            # write out 
            logfp.write('Writing out model.\n')
            cgm.write_pdb(self.structure.location + '/' + basefname + '.pdb')
            cgm.write_rtf(self.structure.location + '/' + basefname + '.rtf')
            cgm.write_prm(self.structure.location + '/' + basefname + '.prm')

            # add the topology and parameter to each segment ... since we
            # never want to build CG model segments individually, we just
            # set these all to built
            for wseg in self.segments.all():
                wseg.isBuilt = 't'
                wseg.rtf_list = self.structure.location + '/' + basefname + '.rtf'
                wseg.prm_list = self.structure.location + '/' + basefname + '.prm'
                wseg.save()

            self.pdbname = self.structure.location + '/' + basefname + '.pdb'
            self.save()
            
            logfp.close()
        else:
            raise AssertionError('Unknown CG model')

    def getAppendPatches():
        """
        This method is overwritten because we don't do any post-appending
        patches on the Go or BLN models, so we want to make sure that we
        return an empty string
        """
        return ''

class Patch(models.Model):
    structure   = models.ForeignKey(WorkingStructure)

    # patches can cross multiple segments, if this field is set, it means
    # that the patch only applies to a particular segment (i.e. it is a
    # protonation patch that should be handled at segment build time rather
    # than append time.
    patch_segid  = models.ForeignKey(Segment,null=True)

    patch_name   = models.CharField(max_length=10)
    patch_segres = models.CharField(max_length=100)

class Task(models.Model):
    STATUS_CHOICES = (
       ('I', 'Inactive'),
       ('Q', 'Queued'),
       ('R', 'Running'),
       ('C', 'Done'),
       ('F', 'Failed'),
       ('K', 'Killed'),
    )

    workstruct  = models.ForeignKey(WorkingStructure)
    action      = models.CharField(max_length=100,null=True) # what this task is doing
    parent      = models.ForeignKey('self',null=True)
    status      = models.CharField(max_length=1, choices=STATUS_CHOICES)
    jobID       = models.PositiveIntegerField()
    scripts     = models.CharField(max_length=500,null=True)
    active      = models.CharField(max_length=1)
    finished    = models.CharField(max_length=1)

    @property
    def scriptList(self):
        if self.scripts.startswith(','):
            self.scripts = self.scripts[1:]
            self.save()
        return self.scripts.split(',')

    def start(self,**kwargs):
        st = self.workstruct.structure

        logfp = open('/tmp/startjob.txt', 'w')
        logfp.write('In job start routine.\n')

        si = schedInterface()
        if kwargs.has_key('altexe'):
            logfp.write('Got alt exe.\n')
            exedict = {}
            for inpscript in self.scriptList:
                exedict[inpscript] = kwargs['altexe']
            logfp.write('exedict = %s\n' % exedict)
            self.jobID = si.submitJob(st.owner.id,st.location,self.scriptList,exedict=exedict)
        else:
            logfp.write('No alt exe.\n')
            self.jobID = si.submitJob(st.owner.id,st.location,self.scriptList)
        logfp.close()
        if self.jobID > 0:
            self.save()
            self.query()
        else:
            raise AssertionError('Job submission fails')

    def kill(self): 
        si = schedInterface()

        if self.jobID > 0:
            si.doJobKill(self.workstruct.structure.owner.id,self.jobID)
            self.status = 'K'

    def query(self):
 
        si = schedInterface()
        if self.jobID > 0:
            sstring = si.checkStatus(self.jobID).split()[4]

            if sstring == 'submitted' or sstring == 'queued':
                self.status = 'Q'
            elif sstring == 'running':
                self.status = 'R'
            elif sstring == 'complete':
                self.status = 'C'
            elif sstring == 'failed':
                self.status = 'F'
            else:
                raise AssertionError('Unknown status ' + sstring)

            self.save()
            return sstring
        else:
            return 'unknown'

    def finish(self):
        """
        This is a pure virtual function. I expect it to be overridden in
        the subclasses of Task. In particular, it needs to decide whether
        the action succeeded or failed.
        """
        pass

    def setup(self,ws):
        models.Model.__init__(self) # call base class init

        self.status = 'I'
        self.jobID = 0
        self.finished = 'n'
        self.scripts = ''
        self.workstruct = ws

class WorkingFile(models.Model):
    path        = models.CharField(max_length=160)
    canonPath   = models.CharField(max_length=160) # added so we can do file versioning
    version     = models.PositiveIntegerField(default=1)
    type        = models.CharField(max_length=20)
    description = models.CharField(max_length=500,null=True)
    task        = models.ForeignKey(Task,null=True)
    pdbkey      = models.CharField(max_length=100,null=True) # if this is a structure file, want to know where to find it

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
    gm_kBond = forms.CharField(max_length=10,initial="50.0")
    gm_kAngle = forms.CharField(max_length=10,initial="30.0")

    # BLN model stuffs
    bln_dm_file = forms.FileField()
    bln_dm_string = forms.CharField(max_length=25)
    bln_nScale = forms.CharField(max_length=10,initial="1.0")
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

class energyTask(Task):
    finale = models.FloatField(null=True)
    usepbc = models.CharField(max_length=1)
    useqmmm = models.CharField(max_length=1)
    qmmmsel = models.CharField(max_length=250)

    def finish(self):
        """test if the job suceeded, create entries for output"""

        loc = self.workstruct.structure.location
        bnm = self.workstruct.identifier

        # There's always an input file, so create a WorkingFile
        # for it.
        wfinp = WorkingFile()
        wfinp.task = self
        wfinp.path = loc + '/' + bnm + '-ener.inp'
        wfinp.canonPath = wfinp.path
        wfinp.type = 'inp'
        wfinp.description = 'energy input script'
        wfinp.save()


        # Check if an output file was created and if so create
        # a WorkingFile for it.
        try:
            os.stat(loc + '/' + bnm + '-ener.out')
        except:
            self.status = 'F'
            return

        wfout = WorkingFile()
        wfout.task = self
        wfout.path = loc + '/' + bnm + '-ener.out'
        wfout.canonPath = wfout.path
        wfout.type = 'out'
        wfout.description = 'energy output'
        wfout.save()

        if self.status == 'F':
            return

        self.status = 'C'

class goModel(models.Model):
    selected = models.CharField(max_length=1)
    pdb = models.ForeignKey(Structure)
    contactType = models.CharField(max_length=10)
    nScale      = models.DecimalField(max_digits=6,decimal_places=3)
    kBond       = models.DecimalField(max_digits=6,decimal_places=3)
    kAngle      = models.DecimalField(max_digits=6,decimal_places=3)

class blnModel(models.Model):
    nScale      = models.DecimalField(max_digits=6,decimal_places=3)
    kBondHelix  = models.DecimalField(max_digits=6,decimal_places=3)
    kBondCoil   = models.DecimalField(max_digits=6,decimal_places=3)
    kBondSheet  = models.DecimalField(max_digits=6,decimal_places=3)
    kAngleHelix = models.DecimalField(max_digits=6,decimal_places=3)
    kAngleCoil  = models.DecimalField(max_digits=6,decimal_places=3)
    kAngleSheet = models.DecimalField(max_digits=6,decimal_places=3)
