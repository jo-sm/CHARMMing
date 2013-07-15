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
import structure
import output
import django.shortcuts, django.http, django.template, django.template.loader
import minimization.views
from django.http import HttpResponseRedirect
from django.contrib import messages
from django.contrib.auth.models import User
from trajanal.models import trajAnalParams, trajanalFileForm
from structure.aux import checkNterPatch
from account.views import isUserTrustworthy
import charmming_config

def parseRMSDMatrix(fname,nframe):
    """Reads elements from the nframe dimensional square matrix
    contained in fname and returns a series of restraints.
    """
    nRead = 0
    retval = ''

    mfp = open(fname, 'r')
    for line in mfp:
        line = line.strip()
        if not line:
            continue
        nRead += 1
        if nRead > nframe:
            raise "barf size mismatch."

        # because the RMSD matrix is symmetric, we toss all but the
        # upper triangular part when making restraint
        larr = line.split()
        for i in range(nRead-1,nframe+1):
            retval += "resd blah blah dist %f\n" % (larr[i])
    mfp.close()
    return retval

def trajanalformdisplay(request):
    logfp = open('/tmp/anal_form.txt', 'a')
    logfp.write('ana form...\n')

    if not request.user.is_authenticated():
        return django.shortcuts.render_to_response('html/loggedout.html')
    #chooses the file based on if it is selected or not
    try:
        file =  structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        messages.error(request, "Please submit a structure first.")
        return HttpResponseRedirect("/charmming/analysis/trajanal/")
#	return django.http.HttpResponse("Please submit a structure first.")

    # do the usual lame copy-paste job...
    solv_pdb = "new_" + file.stripDotPDB(file.filename) + "-solv.pdb"
    neu_pdb = "new_" + file.stripDotPDB(file.filename) + "-neutralized.pdb"
    min_pdb = "new_" + file.stripDotPDB(file.filename) + "-min.pdb"
    md_pdb = "new_" + file.stripDotPDB(file.filename) + "-md.pdb"
    ld_pdb = "new_" + file.stripDotPDB(file.filename) + "-ld.pdb"
    sgld_pdb = "new_" + file.stripDotPDB(file.filename) + "-sgld.pdb"
    filename_list = file.getLimitedFileList('blank')
    seg_list = file.segids.split(' ')
    het_list = file.getNonGoodHetPDBList()
    tip_list = file.getGoodHetPDBList()
    seg2_list = file.getProteinSegPDBList()
    protein_list = [x for x in seg2_list if x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb")]
    go_list = [x for x in seg2_list if x.endswith("-go.pdb") or x.endswith("-go-final.pdb")]
    nucleic_list = [x for x in seg2_list if x.endswith("-dna.pdb") or x.endswith("-dna-final.pdb") or x.endswith("-rna.pdb") or x.endswith("-rna-final.pdb")]
    userup_list = [x for x in seg2_list if not ( x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb") or x.endswith("-dna.pdb") or x.endswith("-dna-final.pdb") \
                     or x.endswith("-rna.pdb") or x.endswith("-rna-final.pdb") or x.endswith("het.pdb") or x.endswith("het-final.pdb") or x.endswith("-go.pdb") \
                     or x.endswith("-go-final.pdb")) ]
    disulfide_list = file.getPDBDisulfides()
    scriptlist = []

    # This handles the case where the user submits the form
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
            logfp.write('processing form...\n')
            logfp.close()
            try:
                if(request.POST['usepatch']):
                    file.handlePatching(request.POST)
            except:
                #If there is no patch, make sure patch_name is zero
                file.patch_name = ""
                file.save()
            if filename != solv_pdb and filename != min_pdb and filename != md_pdb and filename != ld_pdb and filename != sgld_pdb and filename != neu_pdb:
               	minimization.views.append_tpl(request.POST,filename_list,file,scriptlist)
               	html = trajanal_tpl(request,file,scriptlist)
            # treat an appended PDB differently as it will have different PSF/CRD files
            else:
                html = trjanal_tpl(request,file,scriptlist)

            return django.http.HttpResponseRedirect("/charmming/analysis/updatetrajanal/")


    # this handles the part where the user has not uploaded anything (just display the form)
    logfp.write('not processing form...\n')
    logfp.close()
    propatch_list = []
    nucpatch_list = []
    for seg in seg_list:
        if seg.endswith("-pro"):
            defpatch = checkNterPatch(file,seg)
            propatch_list.append((seg,defpatch))
        elif seg.endswith("-dna") or seg.endswith("-rna"):
            # Right now, all nucleic acids get the 5TER and 3TER patches by default
            nucpatch_list.append((seg,"5TER","3TER"))
    trusted = isUserTrustworthy(request.user)
    form = trajanalFileForm(initial={'atomselection': 'all'})
    return django.shortcuts.render_to_response('html/trajanalform.html', {'form': form, 'filename_list': filename_list, 'seg_list':seg_list,'file':file, 'min_pdb': min_pdb, 'solv_pdb':solv_pdb, \
      'neu_pdb': neu_pdb, 'md_pdb': md_pdb, 'ld_pdb': ld_pdb, 'sgld_pdb': sgld_pdb, 'het_list': het_list, 'tip_list':tip_list, 'protein_list':protein_list, \
      'trusted': trusted, 'disulfide_list': disulfide_list, 'propatch_list': propatch_list, 'nucleic_list': nucleic_list, 'nucpatch_list': nucpatch_list, 'seg_list': seg_list, \
      'userup_list': userup_list, 'go_list': go_list})

def trajanal_tpl(request,file,scriptlist):
    logfp = open('/tmp/trajanal.log', 'w')
    logfp.write('in trajanal_tpl\n')

    try:
        oldparam = trajAnalParams.objects.filter(pdb = file, selected = 'y')[0]
        oldparam.selected = 'n'
        oldparam.save()
    except:
        pass
    
    params = trajAnalParams(selected='y')
    post = request.POST
    trjfile = request.FILES['trjfile']['content']
    trjfilename = file.location + "upload_" + file.stripDotPDB(file.filename) + "-trajanal.dcd"
    matrixfilename = file.location + "trajanal_" + file.stripDotPDB(file.filename) + "-matrix.txt"
    inp_out = open(trjfilename, 'w')
    inp_out.write(trjfile)
    inp_out.close()
    logfp.write('wrote out trjfile and matrixfile\n')

    template_dict = {}    
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/trajanal.inp' % charmming_config.charmming_root)
    try:
        template_dict['trajfilename'] = trjfilename
        template_dict['matrixfilename'] = matrixfilename
        template_dict['rmsdyn_start'] = int(post['startstep'])
        template_dict['rmsdyn_stop'] = int(post['stopstep'])
        template_dict['rmsdyn_skip'] = int(post['skip'])
        template_dict['rmsdyn_sele'] = post['atomselection']
    except:
        messages.error(request, "Error in values specified.")
        return HttpResponseRedirect("/charmming/analysis/trajanal/")
#        return django.http.HttpResponse("Error in values specified")
    logfp.write('set-up template dict...\n')

    params.trajFileName = template_dict['trajfilename']
    params.trajStartStep = template_dict['rmsdyn_start']
    params.trajStopStep = template_dict['rmsdyn_stop']
    params.trajSkip = template_dict['rmsdyn_skip']
    params.atomSelect = template_dict['rmsdyn_sele']
    params.save()

    charmm_inp = output.tidyInp(t.render(django.template.Context(template_dict)))
    trajanal_filename = file.location + "charmm-" + file.stripDotPDB(file.filename) + "-trajanal.inp"
    inp_out = open(trajanal_filename, 'w')
    inp_out.write(charmm_inp)
    inp_out.close()
    logfp.write('wrote CHARMM input\n')

    logfp.close()
    return "hello"

def updatetrajanal(request):
    logfp = open('/tmp/udta.txt', 'a')
    logfp.write('updatetrajanal called...\n')
    logfp.close()
    return django.http.HttpResponse("running trajanal")
