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
from account.views import isUserTrustworthy
from structure.models import Structure, goModel 
from structure.qmmm import makeQChem, makeQChem_tpl, handleLinkAtoms, writeQMheader
from structure.editscripts import generateHTMLScriptEdit
from structure.aux import checkNterPatch
from django.contrib.auth.models import User
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from minimization.models import minimizeParams
import charmming_config
import output
import lessonaux
import re
import copy
import os, shutil
import commands

#processes form data for minimization
def minimizeformdisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    Structure.checkRequestData(request)
    #chooses the file based on if it is selected or not
    try:
        file =  Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    os.chdir(file.location)
    seg_list = file.segids.split(' ')
    #descript_seg_list = file.segDescription()
    temp = file.stripDotPDB(file.filename)
    #creates a list of filenames associated with the PDB
    #filename_list = file.getLimitedFileList('minimized')
    filename_list = file.getLimitedFileList('minimization')
    #If the user wants to submit a solvated PDB for minimization
    #then we skip the append step
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
    disulfide_list = file.getPDBDisulfides()
    scriptlist = []

    # Tim Miller: 03-09-2009 we need to pass a patch_list (containing
    # both the list of segments and their default NTER patches) to the
    # template.
    propatch_list = []
    nucpatch_list = []
    for seg in seg_list:
        if seg.endswith("-pro"):
            defpatch = checkNterPatch(file,seg)
            propatch_list.append((seg,defpatch))
        elif seg.endswith("-dna") or seg.endswith("-rna"):
            # Right now, all nucleic acids get the 5TER and 3TER patches by default
            nucpatch_list.append((seg,"5TER","3TER"))

    for i in range(len(filename_list)):
        try:
            tempid = request.POST[filename_list[i]]
	    filename = request.POST[filename_list[i]]
	except:
	    try:
	        tempid = request.POST['solv']
		filename = request.POST['solv']
	    except:
	        tempid = "null"
        if tempid != "null":
            file =  Structure.objects.filter(owner=request.user,selected='y')[0]
	    try:
	        if(request.POST['usepatch']):
		    file.handlePatching(request.POST)
	    except:
	        #If there is no patch, make sure patch_name is zero
	        file.patch_name = ""
		file.save()
	    if filename != solv_pdb and filename != min_pdb and filename != md_pdb and filename != ld_pdb and filename != sgld_pdb and filename != neu_pdb:
	        append_tpl(request.POST,filename_list,file,scriptlist)
	        html = minimize_tpl(request,file,"",scriptlist)
                return HttpResponse(html)
	    #treat a solvated PDB differently as it will have different PSF/CRD files
	    else:
	        html = minimize_tpl(request,file,filename,scriptlist)
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
    return render_to_response('html/minimizeform.html', {'filename_list': filename_list, 'seg_list':seg_list,'file':file, 'min_pdb': min_pdb, 'solv_pdb':solv_pdb, \
      'neu_pdb': neu_pdb, 'md_pdb': md_pdb, 'ld_pdb': ld_pdb, 'sgld_pdb': sgld_pdb, 'het_list': het_list, 'tip_list':tip_list, 'protein_list':protein_list, \
      'docustshake': doCustomShake, 'solv_time_est': solv_time_est, 'unsolv_time_est': unsolv_time_est,'trusted':trusted, 'disulfide_list': disulfide_list, \
      'propatch_list': propatch_list, 'nucleic_list': nucleic_list, 'nucpatch_list': nucpatch_list, 'seg_list': seg_list, 'userup_list': userup_list, 'go_list': go_list, \
      'bln_list': bln_list} )

# go models need to be concatenated together and have their topology and parameter files generated
def handle_gomodel_append(file,gosegs):
    try:
        # only one go model associated with a PDB for now...
        gmp = goModel.objects.filter(pdb = file, selected = 'y')[0]
    except:
        # die noisily
	raise "no go model parameters found"

    # cat together all of the files
    appended_go_pdb = "%s/new_%s-gosegs.pdb" % (file.location, file.stripDotPDB(file.filename))
    outfp = open(appended_go_pdb, 'w')
    for seg in gosegs:
        # we want to operate on the -pro segment, not the -go segment.
        infnm = '%s/new_%s-%s.pdb' % (file.location, file.stripDotPDB(file.filename), seg.replace('-go','-pro'))
        inpfp = open(infnm, 'r')
        for line in inpfp:
            outfp.write(line)
        inpfp.close()
    outfp.close()

    # need to run the alpha_carbon command and append the Go model TOP/PAR so CHARMMing will use it
    if gmp.domainData:
        mycmd = "PYTHONPATH=%s %s/cg_models/alpha_carbon/alpha_carbon.py --input=%s --data=%s/etc/cg --output=%s --stride_path=/usr/local/bin --contacts=%s --nScale=%s --domainScale=%s --domain_string=%s --kBond=%s --kAngle=%s" \
                % (charmming_config.charmming_root,charmming_config.charmming_root,appended_go_pdb,charmming_config.charmming_root,file.location,gmp.contactType,gmp.nScale,gmp.domainScale,gmp.domainData,gmp.kBond,gmp.kAngle)
    else:
        mycmd = "PYTHONPATH=%s %s/cg_models/alpha_carbon/alpha_carbon.py --input=%s --data=%s/etc/cg --output=%s --stride_path=/usr/local/bin --contacts=%s --nScale=%s --domainScale=%f --kBond=%s --kAngle=%s" \
                % (charmming_config.charmming_root,charmming_config.charmming_root,appended_go_pdb,charmming_config.charmming_root,file.location,gmp.contactType,gmp.nScale,gmp.domainScale,gmp.kBond,gmp.kAngle)
    status, output = commands.getstatusoutput(mycmd)
    if status != 0:
        raise "alpha_carbon.py failed!"

    # rename cg_new to new_
    os.rename("%s/cg_new_%s-gosegs.pdb" % (file.location,file.stripDotPDB(file.filename)), appended_go_pdb)

    file.rtf = "cg_new_%s-gosegs.rtf" % file.stripDotPDB(file.filename)
    file.prm = "cg_new_%s-gosegs.prm" % file.stripDotPDB(file.filename)
    file.rtf_append_replace = 'append'
    file.prm_append_replace = 'append'
    file.save()
    return "gosegs"

def append_tpl(postdata,filename_list,file,scriptlist):
   min_seg_list = []
   min_tip_list = []
   min_het_list = []
   go_seg_list = []
   bln_seg_list = []
   all_segids = file.segids.split()
   rtf_prm_dict = file.getRtfPrmPath()
   charmm_inp = ""
   dohbuild = False
   prm_builder = "gennrtf"
   # template dictionary passes the needed variables to the template
   template_dict = {}

   try:
       if postdata['usepatch']:
           use_patch = 1
   except:
       use_patch = 0
   try:
       if postdata['parmgen'] == 'antechamber':
           prm_builder = "antechamber"
   except:
       pass

   for i in range(len(filename_list)):
       try:
           tempdata = postdata[filename_list[i]]
       except:
           tempdata = ''
       #The indexes of filename_list and seg_list correlate so if a filename is chosen
       #then the corresponding index of seg_list should  be minimized
       if tempdata:
           if '-goodhet' in filename_list[i]:
               min_seg_list.append(file.getSegIdFromFilename(filename_list[i]))
           elif '-het' in filename_list[i]:
               min_het_list.append(file.getSegIdFromFilename(filename_list[i]))
           else:
               if '-go' in filename_list[i]:
                   go_seg_list.append(file.getSegIdFromFilename(filename_list[i]))
               if '-bln' in filename_list[i]:
                   bln_seg_list.append(file.getSegIdFromFilename(filename_list[i]))
               min_seg_list.append(file.getSegIdFromFilename(filename_list[i]))

   # sanity check
   if len(go_seg_list) > 0 and len(bln_seg_list) > 0:
       raise "Cannot mix go and BLN models!!!"

   #runs each segment through CHARMM before appending begins
   if len(go_seg_list) > 0:
       finalgoseg = handle_gomodel_append(file,go_seg_list)
   if len(bln_seg_list) > 0:
       # copy RTF and RPM in place
       shutil.copy("%s/etc/cg/bln_topology.rtf" % charmming_config.charmming_root, "%s/bln_topology.rtf" % file.location)
       shutil.copy("%s/bln_new_%s-%s.prm" % (file.location,file.stripDotPDB(file.filename),bln_seg_list[0].replace('-bln','-pro')), "%s/bln_parameter.prm" % (file.location))
       file.rtf = "bln_topology.rtf"
       file.prm = "bln_parameter.prm"
       file.rtf_append_replace = 'replace'
       file.prm_append_replace = 'replace'
       file.doblncharge = 't'
       file.save()
   else:
       file.doblncharge = 'f'
       file.save()

   if min_het_list:
       if not file.rtf or not file.prm:
           if prm_builder == 'antechamber':
               file.makeAntechamber(min_het_list,postdata.has_key("addhyd"))
           else:
               file.makeGenRTF(min_het_list,postdata.has_key("addhyd"))

       #Now that the het atoms have been made CHARMM-compliant
       #we can add the min_het_list to min_seg_list
       min_seg_list = min_seg_list + min_het_list

   if min_seg_list:
       buildcount = 0
       for theSeg in min_seg_list:
           if not file.ifBasicInputExists(theSeg) or postdata.has_key('edit_script'):
               file.makeBasicInput_tpl(theSeg,postdata,scriptlist)
               buildcount += 1
       if buildcount > 0 or use_patch:
           dohbuild = True

   #charmm_inp = file.makeCHARMMInputHeader('Append the PDBs',postdata)
   template_dict['useqmmm'] = postdata.has_key("useqmmm")	 
   template_dict['topology_list'] = file.getTopologyList()
   template_dict['parameter_list'] = file.getParameterList()
   template_dict['filebase'] = file.stripDotPDB(file.filename)
   template_dict['min_seg_list'] = min_seg_list
   template_dict['min_seg_0'] = ''

   if min_seg_list:
       template_dict['min_seg_0'] = min_seg_list[0] 
   
   t = get_template('%s/mytemplates/input_scripts/append_template.inp' % charmming_config.charmming_root)
   # since template limits using certain python functions, here a list of dictionaries is used to pass mult-variables 
   template_dict['het_list'] = []

   for i in range(1,len(min_seg_list)):
       #if(i != 0):
       # het as a dictionary will pass four variables with four keys ('seg', 'good', 'bad' and 'other') for each 'i' to the tamplate 
       het = {}
       het['seg'] = min_seg_list[i]
       het['good'] = ''
       het['bad'] = ''
       het['other'] = ''  

       goodhet = re.compile('goodhet')
       badhet = re.compile('-het')

       #if there are goodhets then the segid is parsed from the list of sections to minimize/append

       if(goodhet.search(min_seg_list[i])):
           segid = file.segidFromGoodhet(min_seg_list[i])
           het['good'] = segid  
       elif(badhet.search(min_seg_list[i])):
           het['bad'] = min_seg_list[i]
       else:
           het['other'] = min_seg_list[i]
       # all het dictionaries are stored in the het_list
       template_dict['het_list'].append(het)

   template_dict['patch_name'] = ''
   if file.patch_name:
       template_dict['patch_name'] = 'true'

       pfp = open(file.location + file.patch_name, 'r')
       #patch_list will store the line starting with 'patch disu' in pfp
       template_dict['patch_list'] = []
       for pline in pfp:
           pline = pline.strip()
           if pline.startswith('patch disu'):
               template_dict['patch_list'].append(pline)
       pfp.close()

   template_dict['dohbuild'] = ''

   if dohbuild:
       template_dict['dohbuild'] = 'true'
       template_dict['min_het_list'] = min_het_list
       if min_het_list:
           template_dict['min_het_0'] = min_het_list[0]

   template_dict['blncharge'] = file.doblncharge == 't'

   user_id = file.owner.id
   os.chdir(file.location)
   append_filename = "charmm-" + file.stripDotPDB(file.filename) + "-append.inp"
   template_dict['headqmatom'] = 'blankme'
   template_dict['output_name'] = 'new_' + file.stripDotPDB(file.filename) + '-final'
   charmm_inp = output.tidyInp(t.render(Context(template_dict)))

   inp_out = open(append_filename,'w')
   inp_out.write(charmm_inp)
   inp_out.close()
   file.append_status = "Done"
   file.save()
   #send to job queue
   #print "Submitting append job."
   #si = schedInterface()
   #si.submitJob(user_id,file.location,append_filename)
   scriptlist.append(file.location + append_filename)

   #This will return a list of segments used
   return min_seg_list


def minimize_tpl(request,file,final_pdb_name,scriptlist):
    postdata = request.POST
    #deals with changing the selected minimize_params
    try:
        oldparam = minimizeParams.objects.filter(pdb = file, selected = 'y')[0]
	oldparam.selected = 'n'
	oldparam.save()
    except:
        pass
    
    #change the status of the file regarding minimization 
    sdsteps = postdata['sdsteps']
    abnr = postdata['abnr']
    tolg = postdata['tolg']
    rtf_prm_dict = file.getRtfPrmPath()
    os.chdir(file.location)
    
    try:
        selectedparam = minimizeParams.objects.filter(pdbfile = file,selected='y')[0]
        selectedparam.selected = 'n'
	selectedparam.save()
    except:
        #No selected minimizeparam
	pass

    # create a model for the minimization
    mp = minimizeParams(selected='y')
    try:
        mp.sdsteps = int(sdsteps)
        mp.abnrsteps = int(abnr)
        mp.tolg = tolg
    except:
        # FixMe: alert user to errors
        return "Error"        

    try:
        if postdata['usepbc']:
            mp.usepbc = 'y'
        else:
            mp.usepbc = 'n'        
    except:
        mp.usepbc = 'n'

    if postdata.has_key('useqmmm'):
        mp.useqmmm = 'y'
        file.checkForMaliciousCode(postdata['qmsele'],postdata)
        try:
            mp.qmmmsel = postdata['qmsele']
        except:
            pass
    else:
        mp.useqmmm = 'n'

    mp.save()

    # final_pdb_name is just shorthand for new_filename-final.pdb/etc
    # min_pdb_name is also shorthand for the name of the minimized pdb
    # If the user submits a PDB that isn't a solvated one
    # then change final_pdb_name to the -final format since appending has
    # yet to occur. Otherwise final_pdb_name will remain the name
    # of the solvated PDB
    apply_direct_patch = 0
    if(final_pdb_name == ""):
        final_pdb_name = "new_" + file.stripDotPDB(file.filename) + "-final"
        min_pdb_name = "new_" + file.stripDotPDB(file.filename) + "-min"
    else:
        #If the user wants a patch applied to a minimized/solvated/md/etc then directly apply
        #Tim Miller: I don't think we should actually be doing this, since all patching should be done
        #by append!
        #if request.POST.has_key('usepatch'):
        #    apply_direct_patch = 1
        final_pdb_name = file.stripDotPDB(final_pdb_name)
        min_pdb_name = "new_" + file.stripDotPDB(file.filename) + "-min"
    # template dictionary passes the needed variables to the template 
    template_dict = {}
    template_dict['topology_list'] = file.getTopologyList()
    template_dict['parameter_list'] = file.getParameterList()
    template_dict['filebase'] = file.stripDotPDB(file.filename)    
    template_dict['input_file'] = final_pdb_name 
    template_dict['restraints'] = '' 
    try:
        postdata['apply_restraints']
        template_dict['restraints'] = file.handleRestraints(request)
    except:
        pass
    
    template_dict['direct_patch'] = apply_direct_patch
    # patch_list passes all lines in patch_file below to the template
    template_dict['patch_list'] = [] 
    if apply_direct_patch:
        patch_file = open(file.location + file.patch_name,'r')
        for line in patch_file:
             template_dict['patch_list'].append(line)
    solvate_implicitly = 0
    try:
        if(postdata['solvate_implicitly']):
            solvate_implicitly = 1
    except:
        pass

    template_dict['solvate_implicitly'] = solvate_implicitly
    template_dict['fixnonh'] = 0 
    try:
        postdata['fixnonh']
        template_dict['fixnonh'] = 1
    except:
        #If there is a topology or paramter file then don't constrain anything
        if file.ifExistsRtfPrm() < 0:
            template_dict['fixnonh'] = 2 

    #handles shake
    template_dict['shake'] = request.POST.has_key('apply_shake')
    if request.POST.has_key('apply_shake'):
        template_dict['which_shake'] = postdata['which_shake']
        template_dict['qmmmsel'] = mp.qmmmsel

        if postdata['which_shake'] == 'define_shake':
            template_dict['shake_line'] = postdata['shake_line']
            if postdata['shake_line'] != '':
                file.checkForMaliciousCode(postdata['shake_line'],postdata)
#            else:
#                if mp.useqmmm == 'y':
#                    template_dict['shake_line'] = "notfix .and. .not. ( %s )" % mp.qmmmsel
#                else:
#                    template_dict['shake_line'] = "notfix"
#        else:
#            if mp.useqmmm == 'y':
#                template_dict['shake_line'] = "notfix .and. .not. ( %s )" % mp.qmmmsel
#            else:
#                template_dict['shake_line'] = "notfix"
           

    template_dict['restraints'] = ''
    try:
        postdata['apply_restraints']
        template_dict['restraints'] = file.handleRestraints(request)
#        charmm_inp = charmm_inp + file.handleRestraints(request)
    except:
        pass

    # check to see if PBC needs to be used -- if so we have to set up Ewald
    template_dict['usepbc'] = ''
    template_dict['dopbc'] = 0
    try:
        if postdata['usepbc']:
            template_dict['usepbc'] =  postdata['usepbc']
            if file.solvation_structure != 'sphere':
                dopbc = 1
                template_dict['dopbc'] = 1
            else:
                dopbc = 0
        else:
            dopbc = 0
    except:
        dopbc = 0

#    template_dict['solvate_implicitly'] = solvate_implicitly 
    template_dict['solvation_structure'] = file.solvation_structure
    template_dict['relative_boundary'] = 0
    if dopbc:
        relative_boundary = 0
        if file.solvation_structure != '' and file.crystal_x < 0:
            relative_boundary = 1

        template_dict['relative_boundary'] = relative_boundary  
        template_dict['dim_x'] = str(file.crystal_x)
        template_dict['dim_z'] = str(file.crystal_z)
        # Tim Miller test
        greaterval = max(file.crystal_x,file.crystal_z)
        template_dict['greaterval'] = str(greaterval)

        # set up images
        if file.solvation_structure == '' or solvate_implicitly:
            pass  
        else:
            # we should have a solvation file to read from
            try:
                os.stat(file.location + "new_" + file.stripDotPDB(file.filename) + ".xtl")
            except:
                # need to throw some sort of error ... for now just toss cookies
                return HttpResponse("Oops ... transfer file not found.")

        # we need to get the ewald parameter
        template_dict['file_location'] = file.location
    template_dict['useqmmm'] = postdata.has_key("useqmmm")
    template_dict['sdsteps'] = sdsteps
    template_dict['abnr'] = abnr
    template_dict['tolg'] = tolg
    template_dict['output_name'] = min_pdb_name
    
    if postdata.has_key("useqmmm"):
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
	if int(postdata['num_linkatoms']) > 0:
	    linkatoms = handleLinkAtoms(file,postdata)
	else:
	    linkatoms = None
        template_dict = makeQChem_tpl(template_dict, file, exch, corr, bs, qmsel, "Force", charge, multi, file.stripDotPDB(final_pdb_name) + ".crd", linkatoms)

    template_dict['headqmatom'] = 'blankme'
    if mp.useqmmm == 'y':
        headstr = writeQMheader("", "SELE " + qmsel + " END")
        template_dict['headqmatom'] = headstr.strip() 
    mp.statusHTML = "<font color=yellow>Processing</font>"
    mp.pdb = file
    mp.save()
    file.save()
    t = get_template('%s/mytemplates/input_scripts/minimization_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))
    
    user_id = file.owner.id
    minimize_filename = file.location + "charmm-" + file.stripDotPDB(file.filename) + "-min.inp"
    inp_out = open(minimize_filename ,'w')
    inp_out.write(charmm_inp)
    inp_out.close()	
    scriptlist.append(minimize_filename)
    file.save()
    if postdata.has_key('edit_script') and isUserTrustworthy(request.user):
        return generateHTMLScriptEdit(charmm_inp,scriptlist,'minimization')
    else:
        si = schedInterface()
        newJobID = si.submitJob(user_id,file.location,scriptlist)
	if file.lesson_type:
            lessonaux.doLessonAct(file,"onMinimizeSubmit",postdata,final_pdb_name)
        file.save()

        if newJobID < 0:
           mp.statusHTML = "<font color=red>Failed</font>"
           mp.save()
        else:
           file.minimization_jobID = newJobID
           sstring = si.checkStatus(newJobID)
           mp.statusHTML = statsDisplay(sstring,newJobID)
           mp.save()
           file.save()   
	return "Done."


