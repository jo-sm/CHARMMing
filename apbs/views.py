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

import re, os, cPickle, copy
import django.shortcuts, django.http, django.template.loader, django.template
import minimization.views, input, output, apbs
import charmming_config, output, scheduler
from  apbs import redox_mod
from django.http import HttpResponse
from structure.models import Structure, WorkingStructure, Task
from apbs.models import redoxTask
from apbs import redox_mod
from pychm.lib.mol import Mol
from pychm.io.rtf import RTFFile

def redoxformdisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    input.checkRequestData(request)

    #chooses the file based on if it is selected or not
    try:
         struct = Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
         return output.returnSunmission("Oxidationb/reduction", error="Please submit a structure first.")
    try:
         ws = WorkingStructure.objects.filter(structure=struct,selected='y')[0]
    except:
        return output.returnSubmission("Oxidation/reduction", error="Please visit the &quot;Build Structure&quot; page to build your structure before minimizing")

    tmpfp = open(struct.pickle, 'r')
    pdb = cPickle.load(tmpfp)
    pdb_metadata = pdb.get_metaData()
    tmpfp.close()

    if request.POST.has_key('picksite'):
        # This code path is taken if the form IS filled in
        try:
            oldtsk = redoxTask.objects.filter(workstruct=ws,active='y')[0]
            oldtsk.active = 'n'
            oldtsk.save()
        except:
            pass
        rdxtsk = redoxTask()
        rdxtsk.setup(ws)
        rdxtsk.action = 'redox'
        rdxtsk.active = 'y'
        rdxtsk.save()

        if ws.isBuilt != 't':
            return output.returnSubmission("Oxidation/reduction", error='Your working structure must be built before you perform a redox calculation')

        isBuilt = True
        pTaskID = Task.objects.get(workstruct=ws,action='build').id
        pdb = createFinalPDB(request,ws)
        return HttpResponse(redox_tpl(request,rdxtsk,ws,pdb,pdb_metadata))

    else:
        # This code path is taken if the form IS NOT filled in. We need
        # to decide which REDOX sites are present in our working structure.
        
        # use charmminglib to decide if there is a valid FeS compound in
        # this structure.
        try:
            myMol = pdb['append_%s' % ws.identifier]
        except:
            return output.returnSubmission("Oxidation/reduction", error='Your working structure must be built before you perform a redox calculation')


        pfp = open(struct.pickle, 'r')
        pdb = cPickle.load(pfp)
        pfp.close()
        # ToDo: fix so that the user does not have to build their WorkingStructure
        # before they can use redox.
        redox_segs = []
        redox_nums = {}

        het = list(myMol.iter_res(segtypes = ['bad'],resName = ['fs4','sf4']))
        for res in het:
            if redox_nums.has_key(res.chainid.upper()):
                redox_nums[res.chainid.upper()] += 1
            else:
                redox_segs.append(res.chainid.upper())
                redox_nums[res.chainid.upper()] = 1


        # if we have any segments that are marked as redox only, we need to take
        # them into consideration ... this is a bit of a hack.
        parsedMdl = pdb[0]
        for seg in ws.segments.all():
            if seg.redox:
                for res in parsedMdl.iter_res():
                    if res.segid == seg.name:
                        if redox_nums.has_key(res.chainid.upper()):
                            redox_nums[res.chainid.upper()] += 1
                        else:
                            redox_segs.append(res.chainid.upper())
                            redox_nums[res.chainid.upper()] = 1


        if len(redox_segs) < 1:
            noredox = True
        else:
            noredox = False

        # Check and see if we have any results to print out
        print_result = False
        calc_final = True
        oxipot = 'N/A'
        oxipotref = 'N/A'
        modpot = 'N/A'
        modpotref = 'N/A'
        delg = "N/A"
        delgnf = "N/A"
        finres = "N/A"

        logfp = open('/tmp/happytim.txt', 'w')

        try:
            logfp.write("looking for " + struct.location + '/redox-' + ws.identifier + '-mopot.txt\n')
            os.stat(struct.location + '/redox-' + ws.identifier + '-modpot.txt')
        except OSError, oe:
            logfp.write("OSError\n")
            calc_final = False
        else:
            logfp.write("OK\n")
            print_result = True
            fp = open(struct.location + '/redox-' + ws.identifier + '-modpot.txt', 'r')
            try:
                modpot = float(fp.readline())
            except:
                print_result = False
            fp.close()

        logfp.write('modpot = %s\n' % modpot)
        logfp.close()

        try:
            os.stat(struct.location + '/redox-' + ws.identifier + '-modpotref.txt')
        except OSError, oe: 
            calc_final = False
        else:
            print_result = True
            fp = open(struct.location + '/redox-' + ws.identifier + '-modpotref.txt', 'r')
            try:
                modpotref = float(fp.readline())
            except:
                print_result = False
            fp.close()
        try:
            os.stat(struct.location + '/redox-' + ws.identifier + '-oxipot.txt')
        except OSError, oe: 
            calc_final = False
        else:
            print_result = True
            fp = open(struct.location + '/redox-' + ws.identifier + '-oxipot.txt', 'r')
            oxipot = float(fp.readline())
            fp.close()
        try:
            os.stat(struct.location + '/redox-' + ws.identifier + '-oxipotref.txt')
        except:
            calc_final = False
        else:
            print_result = True
            fp = open(struct.location + '/redox-' + ws.identifier + '-oxipotref.txt', 'r')
            try:
                oxipotref = float(fp.readline())
            except:
                print_result = False
            fp.close()

        if calc_final:
            # figure out the correct ade value
            rp = apbs.models.redoxTask.objects.filter(workstruct=ws, action='redox', active='y', finished='y')[0]
            ade = 0.0
            delg = 0.0
            if rp.redoxsite == 'couple_oxi':
                delg = modpot - modpotref - oxipot + oxipotref
                ade = 0.273
            elif rp.redoxsite == 'couple_red':
                delg = oxipot - oxipotref - modpot + modpotref
                ade = -3.543

            delgnf = delg * (-4.184/96.485)
            finres = delgnf - 4.43 + ade 
    
        # Django's template language doesn't handle ranges well, so this is a hack to display
        # the correct number of redox sites for each segment.
        rn = {}
        for k in redox_nums.keys():
             rn[k] = range(1,redox_nums[k]+1)
        return django.shortcuts.render_to_response('html/redox.html', {'redox_segs': redox_segs, 'noredox': noredox, 'print_result': print_result, \
                                                                       'oxipot': oxipot, 'oxipotref': oxipotref, 'modpot': modpot, 'modpotref': modpotref, \
                                                                       'delg': delg, 'delgnf': delgnf, 'finres': finres, \
                                                                       'redox_nums': rn, 'ws_identifier': ws.identifier})

def genstruct_old_tpl(request,file,scriptlist):
    selectedlist = []
    selectedres  = {}
    filenamelist = file.getLimitedFileList('foo')

    for fname in filenamelist:
        if request.POST.has_key(fname): selectedlist.append(file.getSegIdFromFilename(fname))
    for elt in selectedlist:
        res = elt[0]
        if res in selectedres.keys():
            selectedres[res] += 1
        else:
            selectedres[res] = 1

    finalseglist = []
    for res in selectedres.keys():
        if selectedres[res] == 2: finalseglist.append(res)

    # This loop builds the input scripts that combine the -pro and -het segments for each of the
    # chains. At the end, they all get appended, oriented, and the protein is dropped to give the
    # the reference struct.
    for segid in finalseglist:
        td = {}
        td['segid'] = segid
        td['topology_list'] = ['%s/toppar/top_all27_prot_na.rtf' % charmming_config.data_home, \
                               '%s/toppar/top_4fsr.rtf' % charmming_config.data_home]
        td['parameter_list'] = ['%s/toppar/par_all27_prot_na.prm' % charmming_config.data_home,
                                '%s/toppar/par_4fsr.prm' % charmming_config.data_home]
        td['filebase'] = file.stripDotPDB(file.filename)

        gen_filename = "redox-" + file.stripDotPDB(file.filename) + "-" + segid + "-setup.inp"
        t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_makeseg_template.inp' % charmming_config.charmming_root)
        charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
        rseg_filename = "redox-" + file.stripDotPDB(file.filename) + "-" + segid + "-join.inp"
        fp = open(file.location + rseg_filename, 'w')
        fp.write(charmm_inp)
        fp.close() 
        scriptlist.append(file.location + rseg_filename)

    # now make the script that will append all of them together
    td = {}
    td['topology_list'] = ['%s/toppar/top_all27_prot_na.rtf' % charmming_config.data_home, \
                           '%s/toppar/top_4fsr.rtf' % charmming_config.data_home]
    td['parameter_list'] = ['%s/toppar/par_all27_prot_na.prm' % charmming_config.data_home,
                            '%s/toppar/par_4fsr.prm' % charmming_config.data_home]
    td['filebase'] = file.stripDotPDB(file.filename)
    td['seglist'] = finalseglist

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_append.inp' % charmming_config.charmming_root)        
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    rapp_filename = "redox-" + file.stripDotPDB(file.filename) + "-append.inp"
    fp = open(file.location + rapp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    scriptlist.append(file.location + rapp_filename)
        
    return finalseglist

def getdelg_tpl(request,workstruct,redoxTask,knockout):
    td = {}
    try:
        td['srad'] = float(request.POST['srad'])
    except:
        td['srad'] = 1.4
   
    td['rtf'] = '%s/toppar/top_all22_4fe4s_esp_090209.inp' % charmming_config.data_home
    td['prm'] = '%s/toppar/par_all22_4fe4s_esp_090209.inp' % charmming_config.data_home
    td['id'] = workstruct.identifier
    td['data_home'] = charmming_config.data_home
    td['filebase'] = 'redox-%s' % workstruct.identifier
    td['knockout'] = knockout

    # step 1: generate final grids for full system and redox site
    # a. full sys
    td['psf'] = 'redox-%s-oxiall.psf' % workstruct.identifier
    td['crd'] = 'redox-%s-oxiall.crd' % workstruct.identifier
    td['sizegrid'] = True
    td['grid_name'] = 'pro'
    try:
        td['pdie'] = float(request.POST['prot_diel'])
        td['sdie'] = float(request.POST['solv_diel'])
    except: 
        td['pdie'] = 4
        td['sdie'] = 78
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_mkgrid.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    inp_name = 'redox-%s-mkgrid-full.inp' % workstruct.identifier
    fp = open('%s/%s' % (workstruct.structure.location,inp_name), 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_name

    # b. redox site
    td['psf'] = 'redox-%s-oxisite.psf' % workstruct.identifier
    td['crd'] = 'redox-%s-oxisite.crd' % workstruct.identifier
    td['grid_name'] = 'ref'
    td['sizegrid'] = False
    try:
    # THE FOLLOWING TWO LINES WERE IN CORRECT - BSP
    #   td['pdie'] = float(request.POST['prot_diel'])
    #   td['sdie'] = float(request.POST['solv_diel'])
        td['pdie'] = float(request.POST['redx_diel'])
        td['sdie'] = float(request.POST['prot_diel'])
    except:
        td['pdie'] = 1
        td['sdie'] = 4
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_mkgrid.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    inp_name = 'redox-%s-mkgrid-site.inp' % workstruct.identifier
    fp = open('%s/%s' % (workstruct.structure.location,inp_name), 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_name

    # step 2: use dxmath to generate the final grids
    td['dxmath_inp'] = []
    td['dxmath'] = '%s/apbs-1.1.0/tools/mesh/dxmath' % charmming_config.data_home
    for axis in ["x", "y", "z"]:
        dxmath_scr = "combo.%s" % axis
        fp = open(workstruct.structure.location + '/' + dxmath_scr, 'w')
        fp.write("# combine mesh for %s\n" % axis)
        fp.write("%s/ref.%s 4.0 - %s/pro.%s + %s/iapbs-diel%s.dx =\n" % (workstruct.structure.location,axis,workstruct.structure.location,axis,workstruct.structure.location,axis))
        fp.close()
        td['dxmath_inp'].append(dxmath_scr)
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_combinegrid.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    inp_name = 'redox-%s-combinegrid-site.inp' % workstruct.identifier
    fp = open('%s/%s' % (workstruct.structure.location,inp_name), 'w')
    fp.write(charmm_inp)  
    fp.close()
    redoxTask.scripts += ',%s' % inp_name

    td['sizegrid'] = True

    # step 3: for each of redpot, redpotref, modpot, and modpotref generate the 
    # delta G via APBS

    td['data_home'] = charmming_config.data_home

    # oxipot
    td['input_file'] = 'redox-%s-oxiall' % workstruct.identifier
    td['rdiel'] = 'rdiel'
    td['pdie'] = 4
    td['sdie'] = 78
    td['enefile'] = 'redox-%s-oxipot.txt' % workstruct.identifier
    inp_filename = 'redox-%s-oxipot.inp' % workstruct.identifier

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_delg.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(workstruct.structure.location + '/' + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_filename

    # oxipotref
    td['input_file'] = 'redox-%s-oxisite' % workstruct.identifier
    td['rdiel'] = ''
    td['pdie'] = 1
    td['sdie'] = 1
    td['enefile'] = 'redox-%s-oxipotref.txt' % workstruct.identifier
    inp_filename = 'redox-%s-oxipotref.inp' % workstruct.identifier

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_delg.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(workstruct.structure.location + '/' + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_filename

    # modpot
    td['input_file'] = 'redox-%s-redall' % workstruct.identifier
    td['rdiel'] = 'rdiel'
    td['pdie'] = 4
    td['sdie'] = 78
    td['enefile'] = 'redox-%s-modpot.txt' % workstruct.identifier
    inp_filename = 'redox-%s-modpot.inp' % workstruct.identifier

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_delg.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(workstruct.structure.location + '/' + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_filename

    # modpotref
    td['input_file'] = 'redox-%s-redsite' % workstruct.identifier
    td['rdiel'] = ''
    td['pdie'] = 1
    td['sdie'] = 1
    td['enefile'] = 'redox-%s-modpotref.txt' % workstruct.identifier
    inp_filename = 'redox-%s-modpotref.inp' % workstruct.identifier

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_delg.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(workstruct.structure.location + '/' + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_filename

def genstruct_tpl(workstruct,redoxTask,rsite_chain,cysResList):
    td = {}
    td['rtf'] = '%s/toppar/top_all22_4fe4s_esp_090209.inp' % charmming_config.data_home
    td['prm'] = '%s/toppar/par_all22_4fe4s_esp_090209.inp' % charmming_config.data_home
    td['id'] = workstruct.identifier

    logfp = open('/tmp/cysres,txt', 'w')
    logfp.write('cysResList = %s\n' % cysResList)

    td['segresid'] = ''
    resnums = cysResList.split(',')
    if len(resnums) != 4:
        raise AssertionError('wrong number of elements is cysResList')
    for rnum in resnums:
        td['segresid'] += '%s-pro %s ' % (rsite_chain,rnum)
    td['segresid'] += '%s-bad 1' % rsite_chain # FixMe: we hope for now that the REDOX site is the only residue in the bad chain

    logfp.write('segresid = %s\n' % td['segresid'])
    logfp.close()

    if redoxTask.redoxsite == 'couple_oxi':
        oxipatch = 'PFSO'
        redpatch = 'PFSR'
    else:
        oxipatch = 'PFSR'
        redpatch = 'PFSS' 

    # step 1: oxidized all
    td['segs'] = [ seg.name for seg in workstruct.segments.all()]
    td['suffix'] = '_o'
    td['outname'] = 'redox-%s-oxiall' % workstruct.identifier
    td['needpatch'] = True
    td['presname'] = oxipatch

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_makestruct.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    inp_filename = 'redox-%s-build-oxiall.inp' % workstruct.identifier
    fp = open(workstruct.structure.location + '/' + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_filename

    # step 2: oxidized redox site only
    td['segs'] = [ '%s-bad' % rsite_chain.lower() ]   
    td['outname'] = 'redox-%s-oxisite' % workstruct.identifier
    td['needpatch'] = False
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_makestruct.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    inp_filename = 'redox-%s-build-oxisite.inp' % workstruct.identifier
    fp = open(workstruct.structure.location + '/' + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_filename

    # step 3: reduced all
    td['segs'] = [seg.name for seg in workstruct.segments.all()]
    td['suffix'] = '_r'
    td['outname'] = 'redox-%s-redall' % workstruct.identifier
    td['needpatch'] = True
    td['presname'] = redpatch
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_makestruct.inp' % charmming_config.charmming_root)
    inp_filename = 'redox-%s-build-redall.inp' % workstruct.identifier
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(workstruct.structure.location + '/' + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_filename

    # step 4: reduced redox site only
    td['segs'] = [ '%s-bad' % rsite_chain.lower() ]
    td['outname'] = 'redox-%s-redsite' % workstruct.identifier
    td['needpatch'] = False
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_makestruct.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    inp_filename = 'redox-%s-build-redsite.inp' % workstruct.identifier
    fp = open(workstruct.structure.location + '/' + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    redoxTask.scripts += ',%s' % inp_filename

def redox_tpl(request,redoxTask,workstruct,pdb,pdb_metadata):

    if not request.POST.has_key('couple'):
        raise AssertionError('No coupling found')

    if not request.POST.has_key('picksite'):
        raise AssertionError('Got to redox_tpl without picksite. How did that happen??!')

    doKnockout = ''
    if request.POST.has_key('chrgknockout'):
        #input.checkForMaliciousCode(request.POST['chrgknockout'])
        if len(request.POST['chrgknockout']) > 0:
            doKnockout = request.POST['chrgknockout']
    

    m = re.search('site_([A-Z])_([1-9])', request.POST['picksite'])
    if not m:
        raise AssertionError('picksite is not in the correct format')

    if request.POST['couple'] == 'oxi':
        clusnameo = '4fso'
        redoxTask.redoxsite = 'couple_oxi'
    else:
        # assume we want a reduction
        clusnameo = '4fsr'
        redoxTask.redoxsite = 'couple_red' 

    rtf = RTFFile('%s/toppar/top_all22_4fe4s_esp_090209.inp' % charmming_config.data_home)

    # Step 1: generate the patched structures of the redox site ... a call to Scott's
    # code replaces the old genstruct_tpl call.
    cysDict = {}
    redox_mod.fesSetup(pdb,clusnameo,rtf,int(m.group(2)),0,workstruct.structure.location,workstruct.identifier,pdb_metadata,cysDict)

    dictkey = '%s-%s' % (m.group(1).lower(),m.group(2))
    try:
        cysResList = cysDict[dictkey]
    except:
        raise AssertionError('key %s not found in cysDict -- bad Scott!' % dictkey)

    # step 2: make final PSF and CRD of the oxidized and reduced structures
    genstruct_tpl(workstruct,redoxTask,m.group(1),cysResList)

    # step 3: make dielectric grids (1. protein + SF4 2. just SF4 + hanging -CH2)
    # This script will include a system call to dxmath to make the grids and combine
    # them.
    ##NB: This functionality has been moved into gendelg_tpl -- schedule for removal
    ##gengrid_tpl(request,workstruct,redoxTask)

    # step 4: Get the four free energy values for redpot, redpotref, modpot, modpotref where mod = oxi or sr
    getdelg_tpl(request,workstruct,redoxTask,doKnockout)

    # all scripts generated, submit to the scheduler
    redoxTask.start(altexe=charmming_config.charmm_apbs_exe)
    redoxTask.save()

    return output.returnSubmission('Oxidation/reduction')

def createFinalPDB(request,workstruct):
    if not request.POST.has_key('picksite'):
        raise AssertionError('Variable picksite must be specified')

    pfp = open(workstruct.structure.pickle, 'r')
    pdb = cPickle.load(pfp)
    pfp.close()
    ormdl = pdb[0] # first model has all of the segments
    mymdl = pdb['append_%s' % workstruct.identifier]
    numdl = copy.deepcopy(mymdl)

    for seg in workstruct.segments.all():
        if seg.redox:
            for mdlseg in ormdl.iter_seg():
                if mdlseg.segid == seg.name:
                    # add this segment to the model
                    numdl = numdl + Mol(mdlseg)

    numdl = Mol(numdl)
    pdb['redox_%s' % workstruct.identifier] = numdl

    os.unlink(workstruct.structure.pickle) # to be safe
    pfp = open(workstruct.structure.pickle, 'w')
    cPickle.dump(pdb,pfp)
    pfp.close()

    return numdl
