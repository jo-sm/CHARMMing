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
from django import forms
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from pdbinfo.models import PDBFile, PDBFileForm 
from pdbinfo.qmmm import makeQChem, makeQChem_tpl, writeQMheader
from minimization.views import append_tpl
from normalmodes.aux import getNormalModeMovieNum
from normalmodes.models import nmodeParams
from account.views import isUserTrustworthy
from pdbinfo.editscripts import generateHTMLScriptEdit
from pdbinfo.aux import checkNterPatch
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
import output
import re, copy, os, shutil
import lessonaux, charmming_config

#returns atom number
def calcAtomNumber(request):
    PDBFile().checkRequestData(request)
    try:
        file =  PDBFile.objects.filter(owner=request.user,selected='y')[0]
        os.chdir(file.location)
    except:
        return HttpResponse("Please submit a structure first.")
    try:
        time_filenames = request.POST['time_filenames']
    except:
        return HttpResponse('')
    time_filenames = time_filenames.rstrip(',')
    if time_filenames == '':
        return HttpResponse('')
    time_list = time_filenames.split(',')
    total_atom_num = 0
    for pdb in time_list: 
        total_atom_num += file.countAtomsInPDB(pdb)
    x = total_atom_num
    return HttpResponse(str(int(x)))


#processes form data for ld simulations
def normalmodesformdisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    PDBFile().checkRequestData(request)
    #chooses the file based on if it is selected or not
    try:
        file =  PDBFile.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    os.chdir(file.location)
    #creates a list of filenames associated with the PDB
    #The "md" option takes away all -md,-ld, and -sgld names
    #filename_list = file.getLimitedFileList("md")
    filename_list = file.getLimitedFileList('blank')
    pdbtoobig_list = []
    count = 0
    solv_pdb = "new_" + file.stripDotPDB(file.filename) + "-solv.pdb"
    neu_pdb = "new_" + file.stripDotPDB(file.filename) + "-neutralized.pdb"
    min_pdb = "new_" + file.stripDotPDB(file.filename) + "-min.pdb"
    md_pdb = "new_" + file.stripDotPDB(file.filename) + "-md.pdb"
    ld_pdb = "new_" + file.stripDotPDB(file.filename) + "-ld.pdb"
    sgld_pdb = "new_" + file.stripDotPDB(file.filename) + "-sgld.pdb"
    het_list = file.getNonGoodHetPDBList()
    tip_list = file.getGoodHetPDBList()
    seg2_list = file.getProteinSegPDBList()
    protein_list = [x for x in seg2_list if x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb")]
    go_list = [x for x in seg2_list if x.endswith("-go.pdb") or x.endswith("-go-final.pdb")]
    bln_list = [x for x in seg2_list if x.endswith("-bln.pdb") or x.endswith("-bln-final.pdb")]
    nucleic_list = [x for x in seg2_list if x.endswith("-dna.pdb") or x.endswith("-dna-final.pdb") or x.endswith("-rna.pdb") or x.endswith("-rna-final.pdb")]
    userup_list = [x for x in seg2_list if not ( x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb") or x.endswith("-dna.pdb") or x.endswith("-dna-final.pdb") \
                     or x.endswith("-rna.pdb") or x.endswith("-rna-final.pdb") or x.endswith("het.pdb") or x.endswith("het-final.pdb") or x.endswith("-go.pdb") \
                     or x.endswith("-go-final.pdb") or x.endswith("-bln.pdb") or x.endswith("-bln-final.pdb") ) ]
    seg_list = file.segids.split(' ')
    disulfide_list = file.getPDBDisulfides()

    # Tim Miller: 03-09-2009 we need to pass a segpatch_list (containing
    # both the list of segments and their default NTER patches) to the
    # template.
    propatch_list = []
    nucpatch_list = []
    for seg in seg_list:
        if seg.endswith("-pro"):
            defpatch = checkNterPatch(file,seg)
            propatch_list.append((seg,defpatch))
        elif seg.endswith("-dna") or seg.endswith("-rna"):
            nucpatch_list.append((seg,"5TER","3TER"))

    nmode_lines = ''
    try:
        nmodefp = open('10modes-' + file.stripDotPDB(file.filename) + '.txt','r')
        for line in nmodefp:
            nmode_lines += line
        nmodefp.close()
    except:
        pass
 
    
    #Each checkbox or radio button on the ldform.html page has a name equal
    #to the filename it's associated with. This for loop tries to get the post
    #data using the name as a key
    tempid = None
    for i in range(len(filename_list)):
        try:
            tempid = request.POST[filename_list[i]]
            break
        except:
            try:
                tempid = request.POST['solv_or_min']
                break
            except:
                tempid = None

    if(tempid):
        # if the name is solv/min pdb, it does not need to be appended/patched
        # otherwise reappend them
        scriptlist = []
        if tempid != min_pdb and tempid != solv_pdb and tempid != md_pdb and tempid != ld_pdb and tempid != sgld_pdb and tempid != neu_pdb:
            seg_list = append_tpl(request.POST,filename_list,file,scriptlist)
            ld_this_file = 'new_' + file.stripDotPDB(file.filename) + "-final.pdb"
            html = applynma_tpl(request.POST,file,seg_list,ld_this_file,scriptlist)
            return HttpResponse(html)
	else:
            html = applynma_tpl(request.POST,file,seg_list,tempid,scriptlist)
            return HttpResponse(html)

    doCustomShake = 1
    if file.ifExistsRtfPrm() < 0:
        doCustomShake = 0

    if file.natom > 30000:
        solv_time_est = "at least 4 hours"
        unsolv_time_est = "at least half an hour"
    elif file.natom > 15000:
        solv_time_est = "between 2 to 4 hours"
        unsolv_time_est = "between 15 and 30 minutes"
    elif file.natom > 10000:
        solv_time_est = "from 1 to 2 hours"
        unsolv_time_est = "between 5 and 15 minutes"
    elif file.natom > 5000:
        solv_time_est = "from 30 minutes to an hour"
        unsolv_time_est = "from 1 to 5 minutes"    
    elif file.natom > 1000:
        solv_time_est = "from 15 to 30 minutes"
        unsolv_time_est = "less than a minute" 
    else:
        solv_time_est = "less than 15 minutes"
        unsolv_time_est = "less than a minute" 
    trusted = isUserTrustworthy(request.user)
    return render_to_response('html/normalmodesform.html', {'filename_list': filename_list,'min_pdb':min_pdb,'solv_pdb':solv_pdb, 'neu_pdb': neu_pdb, 'md_pdb': md_pdb, \
      'ld_pdb': ld_pdb,'sgld_pdb': sgld_pdb,'seg_list':seg_list,'protein_list':protein_list,'tip_list':tip_list,'het_list':het_list, 'file': file, 'docustshake': doCustomShake, \
      'solv_time_est': solv_time_est, 'unsolv_time_est': unsolv_time_est, 'nmode_lines':nmode_lines,'trusted':trusted, 'disulfide_list': disulfide_list, 'propatch_list': propatch_list, \
      'nucpatch_list': nucpatch_list, 'nucleic_list': nucleic_list, 'userup_list': userup_list, 'go_list': go_list, 'bln_list': bln_list})

#Generates NMA script and runs it
def applynma_tpl(postdata,file,seg_list,min_pdb,scriptlist):
#    charmm_inp = file.makeCHARMMInputHeader('Normal Mode Analysis',postdata)

    try:
        oldparam = nmodeParams.objects.filter(pdb = file, selected = 'y')[0]
        oldparam.selected = 'n'
        oldparam.save()
    except:
        pass
    nmm = nmodeParams()
    # the decimal values need to be filled in thanks to Django's stupidity -- they get reset
    # to correct values (if needed) later
    nmm.rcut = '0.0'
    nmm.kshort = '0.0'
    nmm.klong = '0.0'
    nmm.statusHTML = "<font color=yellow>Processing</font>"

    #nma will either be useenm or usevibran
    nma = postdata['nma'] 

#    charmm_inp = file.makeCHARMMInputHeader('Molecular Dynamics',postdata)
    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'] = file.getTopologyList()
    template_dict['parameter_list'] = file.getParameterList()
    template_dict['filebase'] = file.stripDotPDB(file.filename)
    template_dict['input_file'] = file.stripDotPDB(min_pdb)
    template_dict['nma'] = nma
    # save model
    nmm.pdb = file
    if nma == 'useenm':
        nmm.type = 2
        nmm.rcut = postdata['rcut']
        nmm.kshort = postdata['kshort']
        nmm.klong = postdata['klong']
    else:
        nmm.type = 1
    nmm.nmodes = int(postdata['num_normalmodes'])
    nmm.selected = 'y'

    if nma == 'useenm':
        pass
    else:
        if file.ifExistsRtfPrm() < 0:
            template_dict['rtf_prm'] = 'y'
        template_dict['restraints'] = '' 
        try:
            postdata['apply_restraints']
            template_dict['restraints'] = file.handleRestraints(postdata,request)
        except:
            pass

    numnormalmodes = postdata['num_normalmodes'] 
    typeoption = 'vibran'
    if nma == 'useenm':
        typeoption = 'enm'
        rcut = postdata['rcut']
        kshort = postdata['kshort']
        klong = postdata['klong']
        template_dict['rcut'] = rcut
        template_dict['kshort'] = kshort
        template_dict['klong'] = klong

    template_dict['useqmmm'] = postdata.has_key("useqmmm")
    if nma == 'usevibran' and postdata.has_key("useqmmm"):
        # validate input
        if postdata['qmmm_exchange'] in ['HF','B','B3']:
            exch = postdata['qmmm_exchange']
        else:
            exch = 'HF'
        if postdata['qmmm_correlation'] in ['None','LYP']:
            corr = postdata['qmmm_correlation']
        else:
            corr = 'None'
        if postdata['qmmm_basisset'] in ['STO-3G','3-21G*','6-31G*']:
            bs = postdata['qmmm_basisset']
        else:
            bs = 'sto3g'
        file.checkForMaliciousCode(postdata['qmsele'],postdata)
        qmsel = postdata['qmsele']
        if qmsel == '':
            qmsel = 'resid 1'
        if postdata['qmmm_charge'] in ['-5','-4','-3','-2','-1','0','1','2','3','4','5']:
            charge = postdata['qmmm_charge']
        else:
            charge = '0'
        if postdata['qmmm_multiplicity'] in ['0','1','2','3','4','5','6','7','8','9','10']:
            multi = postdata['qmmm_multiplicity']
        else:
            multi = '0'
        try:
            if int(postdata['num_linkatoms']) > 0:
                linkatoms = handleLinkAtoms(file,postdata)
            else:
                linkatoms = None
        except:
            linkatoms = None
        template_dict = makeQChem_tpl(template_dict, file, exch, corr, bs, qmsel, "freq", charge, multi, file.stripDotPDB(min_pdb) + ".crd", linkatoms)

    template_dict['numnormalmodes'] = numnormalmodes

    # If need be, print out trajectories for the modes the user requested.
    template_dict['gen_trj'] = postdata.has_key("gen_trj") 
    template_dict['num_trjs'] = postdata.has_key("num_trjs") 
    if postdata.has_key("gen_trj") and postdata.has_key("num_trjs"):
        try:
	    numtrj = getNormalModeMovieNum(file)
	    for i in range(numtrj+1):
	        i = i+1
	        os.remove(file.location +  "new_" + file.stripDotPDB(file.filename) + "-mtraj_" + str(i) + ".trj")
		os.remove(file.location + 'new_' + file.stripDotPDB(file.filename) + '-nma-mainmovie-' + str(i) + '.pdb')
	except:
	    pass
        file.nma_movie_req = True
        nmm.nma_movie_req = True
        try:
            ntrjs = int(postdata['num_trjs'])
        except:
            ntrjs = 5
        if ntrjs > 20:
            ntrjs = 20
        if ntrjs > int(numnormalmodes):
            ntrjs = int(numnormalmodes)
        if postdata.has_key("useqmmm"):
            headstr = writeQMheader("", "SELE " + qmsel + " END")
        else:
            headstr = "* Not using QM/MM\n"
        template_dict['ntrjs'] = str(ntrjs) 
        template_dict['headstr'] = headstr
    else:
        nmm.nma_movie_req = False
    t = get_template('%s/mytemplates/input_scripts/applynma_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))    

    nmm.save()
    user_id = file.owner.id
    os.chdir(file.location)
    nma_filename = "charmm-" + file.stripDotPDB(file.filename) + "-nmodes.inp"
    inp_out = open(file.location + nma_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()  
    #change the status of the file regarding minimization 
    scriptlist.append(file.location + nma_filename)
    if postdata.has_key("gen_trj") and postdata.has_key("num_trjs"):
        return makeNmaMovie_tpl(file,postdata,min_pdb,scriptlist,int(postdata['num_trjs']),typeoption)
#        return makeNmaMovie(file,postdata,min_pdb,scriptlist,int(postdata['num_trjs']),typeoption)
    else:
        if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
            return generateHTMLScriptEdit(charmm_inp,scriptlist,'nma')
	else:
            si = schedInterface()
            newJobID = si.submitJob(user_id,file.location,scriptlist)
            if file.lesson_type:
                lessonaux.doLessonAct(file,postdata,"onNMASubmit")
            file.nma_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            nmm.statusHTML = statsDisplay(sstring,newJobID)
	    nmm.save()

    file.save()
    return "Done."

#Generates NMA script and runs it
def applynma(postdata,file,seg_list,min_pdb,scriptlist):
    charmm_inp = file.makeCHARMMInputHeader('Normal Mode Analysis',postdata)

    try:
        oldparam = nmodeParams.objects.filter(pdb = file, selected = 'y')[0]
        oldparam.selected = 'n'
        oldparam.save()
    except:
        pass
    nmm = nmodeParams()
    nmm.statusHTML = "<font color=yellow>Processing</font>"

    #nma will either be useenm or usevibran
    nma = postdata['nma'] 
    charmm_inp = file.makeCHARMMInputHeader('Molecular Dynamics',postdata)

    # save model
    nmm.pdb = file
    if nma == 'useenm':
        nmm.type = 2
        nmm.rcut = float(postdata['rcut'])
        nmm.kshort = float(postdata['kshort'])
        nmm.klong = float(postdata['klong'])
    else:
        nmm.type = 1
    nmm.nmodes = int(postdata['num_normalmodes'])
    nmm.selected = 'y'

    if nma == 'useenm':
        charmm_inp += """
! Read special ENM atom types (masses are set to 0 here but will be assigned by
! hand in the coarsegrain stream file).
READ RTF CARD APPEND
* ENM specific topology information
*

MASS    98 CGT    0.00000 C ! Coarse grained atom type for ENM
MASS    99 DUM    0.00000 H ! dummy atom

End"""

    charmm_inp = charmm_inp + """

open read unit 2 card name """ + file.stripDotPDB(min_pdb) + """.psf
read psf card unit 2
close unit 2

open read unit 2 card name """ + file.stripDotPDB(min_pdb) + """.crd
read coor card unit 2
close unit 2
"""

    if nma == 'useenm':
        charmm_inp += """
! Waters do not work with ENM, delete them
delete atom select resname TIP3 end
"""
    else:
        if file.ifExistsRtfPrm() < 0:
	    charmm_inp += """
! We do not want to allow atoms not known natively by the CHARMM force field (i.e. those generated
! by GENRTF) to move during the minimization. Therefore, we define all residues that CHARMM's force
! field does know about and use a constraint to fix all other atoms.
DEFI notfix select RESN ALA .or. resn GLU .or. resn GLN .or. resn ASP .or. resn ASN .or. resn LEU -
 .or. resn GLY .or. resn LYS .or. resn SER .or. resn VAL .or. resn ARG .or. resn THR .or. resn PRO -
 .or. resn ILE .or. resn MET .or. resn PHE .or. resn TYR .or. resn CYS .or. resn TRP .or. resn HSD -
 .or. resn HSE .or. resn HSP .or. resn TIP3 .or. resn ZN2 .or. resn SOD .or. resn MG .or. resn FE -
 .or. resn CAL .or. resn CLA .or. resn POT .or. resn CES .or. resn GUA .or. resn ADE .or. resn CYT -
 .or. resn THY .or. resn URA .or. resn LSN .or. resn ASPP .or. resn GLUP end

! we can't use cons fix in case we're doing ewald, but we can do a harmonic constraint with an absurdly high
! force constant
cons harm bestfit force 10000. sele .not. notfix end


!print initial energy
energy

"""

        try:
            postdata['apply_restraints']
            charmm_inp = charmm_inp + file.handleRestraints(postdata,request)
        except:
            pass

    numnormalmodes = postdata['num_normalmodes'] 
    typeoption = 'vibran'
    if nma == 'useenm':
        typeoption = 'enm'
        rcut = postdata['rcut']
        kshort = postdata['kshort']
        klong = postdata['klong']
        charmm_inp += """
! Define CGREG to be all atoms. The stream file
! will delete all of these (in this case, setting the
! entire model to be coarse grained).
DEFIne CGREG SELEct ALL END

! Define CGATS to be all of the alpha carbons. These
! will be the atoms that get left behind from the
! deletion. 
DEFIne CGATS SELEct ATOM * * CA END

! I've set the VDW type for a coarse-grained atom to be
! 98 in the custom RTF. We need to tell the stream this.
SET cgtype 98

! define CG centers with appropriate mass
STREAM /usr/local/charmming/coarsegrain.str

! I think CGATS needs to be reset...
DEFIne CGATS SELEct ATOM * * CA END

! Set weighting to 1 for all selected coarse grained atoms. This
! should activate them all in the model.
SCALar WMAIN SET 1 select CGATS end

! From the ENM restraints script...
! prior to invoking, set value @RCUT to the max distance for restraints.
! prior to invoking, set value @KSHORT to the near-neighbor force constant.
! prior to invoking, set value @KLONG to the long-range force constant.
! prior to invoking, set value @ADDB to nonzero to add bonds for graphics
SET RCUT """ + rcut + """ 
SET KSHORT """ + kshort + """
SET KLONG """ + klong + """
SET ADDB 0

STREAM /usr/local/charmming/enmrestraint-multi.str

SKIP ALL EXCL RESD
"""

    if nma == 'usevibran' and postdata.has_key("useqmmm"):
        # validate input
        if postdata['qmmm_exchange'] in ['HF','B','B3']:
            exch = postdata['qmmm_exchange']
        else:
            exch = 'HF'
        if postdata['qmmm_correlation'] in ['None','LYP']:
            corr = postdata['qmmm_correlation']
        else:
            corr = 'None'
        if postdata['qmmm_basisset'] in ['STO-3G','3-21G*','6-31G*']:
            bs = postdata['qmmm_basisset']
        else:
            bs = 'sto3g'
        file.checkForMaliciousCode(postdata['qmsele'],postdata)
        qmsel = postdata['qmsele']
        if qmsel == '':
            qmsel = 'resid 1'
        if postdata['qmmm_charge'] in ['-5','-4','-3','-2','-1','0','1','2','3','4','5']:
            charge = postdata['qmmm_charge']
        else:
            charge = '0'
        if postdata['qmmm_multiplicity'] in ['0','1','2','3','4','5','6','7','8','9','10']:
            multi = postdata['qmmm_multiplicity']
        else:
            multi = '0'
        try:
            if int(postdata['num_linkatoms']) > 0:
                linkatoms = handleLinkAtoms(file,postdata)
            else:
                linkatoms = None
        except:
            linkatoms = None
        charmm_inp = makeQChem(charmm_inp, file, exch, corr, bs, qmsel, "freq", charge, multi, file.stripDotPDB(min_pdb) + ".crd", linkatoms)

    charmm_inp += "vibran nmode %s\n" % numnormalmodes
    charmm_inp += "diagonalize\n"
    charmm_inp += "print norm\n"
    charmm_inp += "open write unit 10 card name new_%s-nmodes.txt\n" % file.stripDotPDB(file.filename)
    charmm_inp += "write norm card unit 10\n"
    charmm_inp += "* Normal mode vectors\n"
    charmm_inp += "*\n\n"

    # If need be, print out trajectories for the modes the user requested.
    if postdata.has_key("gen_trj") and postdata.has_key("num_trjs"):
        try:
	    numtrj = getNormalModeMovieNum(file)
	    for i in range(numtrj+1):
	        i = i+1
	        os.remove(file.location +  "new_" + file.stripDotPDB(file.filename) + "-mtraj_" + str(i) + ".trj")
		os.remove(file.location + 'new_' + file.stripDotPDB(file.filename) + '-nma-mainmovie-' + str(i) + '.pdb')
	except:
	    pass
        file.nma_movie_req = True
        nmm.nma_movie_req = True
        try:
            ntrjs = int(postdata['num_trjs'])
        except:
            ntrjs = 5
        if ntrjs > 20:
            ntrjs = 20
        if ntrjs > int(numnormalmodes):
            ntrjs = int(numnormalmodes)
        if postdata.has_key("useqmmm"):
            headstr = writeQMheader("", "SELE " + qmsel + " END")
        else:
            headstr = "* Not using QM/MM\n"
        charmm_inp += """

! write trajectories for the lowest """ + str(ntrjs) + """ normal modes at 300K
set n = 1
label wrloop
if n .gt. """ + str(ntrjs) + """ then goto loopend
open write unit 10 unform name new_""" + file.stripDotPDB(file.filename) + """-mtraj_@n.trj
write trajectory unit 10 mode @n temp 300.
* Trajectory for normal mode @n
*

calc n = @n + 1
goto wrloop
label loopend

! save out 

end

! write out the psf
open write unit 1 card name new_""" + file.stripDotPDB(file.filename) + """-nmodes.psf
write psf card unit 1
* PSF from normal modes
""" + headstr + """*

stop
"""
    else:
        nmm.nma_movie_req = False
        charmm_inp += """
end

stop
"""
    nmm.save()
    user_id = file.owner.id
    os.chdir(file.location)
    nma_filename = "charmm-" + file.stripDotPDB(file.filename) + "-nmodes.inp"
    inp_out = open(file.location + nma_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()  
    #change the status of the file regarding minimization 
    scriptlist.append(file.location + nma_filename)
    if postdata.has_key("gen_trj") and postdata.has_key("num_trjs"):
        return makeNmaMovie(file,postdata,min_pdb,scriptlist,int(postdata['num_trjs']),typeoption)
    else:
        if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
            return generateHTMLScriptEdit(charmm_inp,scriptlist,'nma')
	else:
            si = schedInterface()
            newJobID = si.submitJob(user_id,file.location,scriptlist)
            if file.lesson_type:
                lessonaux.doLessonAct(file,postdata,"onNMASubmit")
            file.nma_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            nmm.statusHTML = statsDisplay(sstring,newJobID)
	    nmm.save()

    file.save()
    return "Done."

# --- end of applynma ---

#makes NMA Movie
def makeNmaMovie_tpl(file,postdata,filename,scriptlist,num_trjs,typeoption):
    nmm = nmodeParams.objects.filter(pdb = file, selected = 'y')[0]
    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'] = file.getTopologyList()
    template_dict['parameter_list'] = file.getParameterList()
    template_dict['filebase'] = file.stripDotPDB(file.filename)
    template_dict['typeoption'] = typeoption
    template_dict['num_trjs'] = str(num_trjs)

    #Run the script through charmm, this is not done under job queue at the moment
    #because the PDBs must be appended together after the above script has been run.
    #Once the DAG and query stuff has been implemented, this is how the following code should
    #be changed to
    #1. Create a new field in the PDBFile object called movie_status
    #2. In the status method, check to see if movie_status is done
    #3. If it is done, call the method combinePDBsForMovie(...): right below the following code.

    t = get_template('%s/mytemplates/input_scripts/makeNmaMovie_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    
    movie_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-nmamovie.inp'
    movie_handle = open(file.location + movie_filename,'w')
    movie_handle.write(charmm_inp)
    movie_handle.close()
    user_id = file.owner.id
    scriptlist.append(file.location + movie_filename)
    if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
        return generateHTMLScriptEdit(charmm_inp,scriptlist,'nma')
    else:
        si = schedInterface()
        newJobID = si.submitJob(user_id,file.location,scriptlist)
        if file.lesson_type:
            lessonaux.doLessonAct(file,"onNMASubmit",postdata,None)
        file.nma_jobID = newJobID
        sstring = si.checkStatus(newJobID)
        nmm.statusHTML = statsDisplay(sstring,newJobID)
        #nmm.nma_movie_status = "<font color=33CC00>Done</font>"
        nmm.save()
        file.save()
        return "Done."

#makes NMA Movie
def makeNmaMovie(file,postdata,filename,scriptlist,num_trjs,typeoption):
    nmm = nmodeParams.objects.filter(pdb = file, selected = 'y')[0]
    charmm_inp = """* NMA Movie making
*

bomlev -1

open read unit 2 card name new_""" + file.stripDotPDB(file.filename) + """-nmodes.psf
read psf card unit 2
close unit 2
"""
    if typeoption == 'enm':
        charmm_inp += """
! Define CGREG to be all atoms. The stream file
! will delete all of these (in this case, setting the
! entire model to be coarse grained).
DEFIne CGREG SELEct ALL END

! Define CGATS to be all of the alpha carbons. These
! will be the atoms that get left behind from the
! deletion. 
DEFIne CGATS SELEct ATOM * * CA END

! I've set the VDW type for a coarse-grained atom to be
! 98 in the custom RTF. We need to tell the stream this.
SET cgtype 98

! define CG centers with appropriate mass
STREAM /usr/local/charmming/coarsegrain.str
"""
    charmm_inp += """
set currtrj 1
calc maxtrj """ + str(num_trjs) + """ + 1
label outerloop
 open read unit 10 file name new_""" + file.stripDotPDB(file.filename) + """-mtraj_@currtrj.trj
 traj firstu 10 skip 1

 set i 1
  label loop
  traj read
  open unit 1 card write name new_""" + file.stripDotPDB(file.filename) + """-nmapremovie@i-@currtrj.pdb
  write coor pdb unit 1
  **Cords
  *
  incr i by 1
  if i lt 13 goto loop
 
 incr currtrj by 1
 if currtrj lt @maxtrj goto outerloop

stop"""
    #Run the script through charmm, this is not done under job queue at the moment
    #because the PDBs must be appended together after the above script has been run.
    #Once the DAG and query stuff has been implemented, this is how the following code should
    #be changed to
    #1. Create a new field in the PDBFile object called movie_status
    #2. In the status method, check to see if movie_status is done
    #3. If it is done, call the method combinePDBsForMovie(...): right below the following code.
    movie_filename = 'charmm-' + file.stripDotPDB(file.filename) + '-nmamovie.inp' 
    movie_handle = open(file.location + movie_filename,'w')
    movie_handle.write(charmm_inp)
    movie_handle.close()
    user_id = file.owner.id
    scriptlist.append(file.location + movie_filename)
    if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
        return generateHTMLScriptEdit(charmm_inp,scriptlist,'nma')
    else:		
        si = schedInterface()
        newJobID = si.submitJob(user_id,file.location,scriptlist)
	if file.lesson_type:
	    lessonaux.doLessonAct(file,"onNMASubmit",postdata,None)
        file.nma_jobID = newJobID
        sstring = si.checkStatus(newJobID)
        nmm.statusHTML = statsDisplay(sstring,newJobID)
        #nmm.nma_movie_status = "<font color=33CC00>Done</font>"
	nmm.save()
        file.save()
        return "Done."

#pre: Requires a PDBFile object
#Combines all the smaller PDBs make in the above method into one large PDB that
#jmol understands
#type is nma
def combineNmaPDBsForMovie(file):
    nmm = nmodeParams.objects.filter(pdb = file, selected = 'y')[0]
    ter = re.compile('TER')
    remark = re.compile('REMARK')
    #movie_handle will be the final movie created
    #mini movie is the PDBs which each time step sepearated into a new PDB
    trj_num = 0
    continueit = True
    
    while(continueit):
        try:
            pdbno = trj_num + 1
            os.stat(file.location + 'new_' + file.stripDotPDB(file.filename) + '-mtraj_' + str(pdbno) + '.trj')
            trj_num += 1 
        except:
	    continueit = False
    currtrjnum = 1
    while currtrjnum < (trj_num+1):
        movie_handle = open(file.location + 'new_' + file.stripDotPDB(file.filename) + "-nma-mainmovie-" + str(currtrjnum) + ".pdb",'a')
        for i in range(12):
            i = i+1
            minimovie_handle = open(file.location +  "new_" + file.stripDotPDB(file.filename) + "-nmapremovie" + `i` + "-" + str(currtrjnum) + ".pdb",'r')
            movie_handle.write('MODEL ' + `i` + "\n")
            for line in minimovie_handle:
                if(not remark.search(line) and not ter.search(line)):
                    movie_handle.write(line)
            movie_handle.write('ENDMDL\n')
            minimovie_handle.close()
	    os.remove(file.location +  "new_" + file.stripDotPDB(file.filename) + "-nmapremovie" + `i` + "-" + str(currtrjnum) + ".pdb")
        currtrjnum += 1

    nmm.nma_movie_status = "<font color=33CC00>Done</font>"
    nmm.make_nma_movie = False
    nmm.save()
    file.save()
