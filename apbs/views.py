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

import re, os
import django.shortcuts, django.http, django.template.loader, django.template
import structure.models
import minimization.views
import apbs.models
import charmming_config, output, scheduler

def redoxformdisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    structure.models.Structure.checkRequestData(request)

    #chooses the file based on if it is selected or not
    try:
        file = structure.models.Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return django.http.HttpResponse("Please submit a structure first.")

    seg2_list = file.getProteinSegPDBList()
    protein_list = [x for x in seg2_list if x.endswith("-pro.pdb") or x.endswith("-pro-final.pdb")]
    het_list = [x for x in seg2_list if x.endswith("-het.pdb") or x.endswith("-het-final.pdb")]
    filename_list = file.getLimitedFileList('foo')
    scriptlist = []

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

        # This code path is taken if the form IS filled in
        if tempid != "null":
            minimization.views.append_tpl(request.POST,filename_list,file,scriptlist)
            return django.http.HttpResponse(redox_tpl(request,file,scriptlist))

    # Figure out which segments can make an FeS4 complex...
    redox_list = []
    redox_segs = []
    redox_nums = {}
    redox_nums['A'] = []
    redox_nums['B'] = []
    redox_nums['C'] = []
    redox_nums['D'] = []
    redox_nums['E'] = []
    redox_nums['F'] = []
    redox_nums['G'] = []
    if not file.fes4list:
        noredox = True
    else:
        noredox = False
        for elt in file.fes4list.split(','):
            segment,numfes = elt.split(':')
            numfes = int(numfes)
            seglet = segment.split('-')[0]
            redox_segs.append(seglet.upper())
            redox_nums[seglet.upper()] = range(1,numfes+1)
            pattern = re.compile('%s(-pro|-het)?(-final)?\.pdb' % seglet)
            for fname in filename_list:
                fname = fname.strip()
                if pattern.search(fname): redox_list.append(fname) 

    # This code path is taken if the form is NOT filled in

    # Check and see if we have any results to print out
    print_result = False
    calc_final = True
    redpot = 'N/A'
    redpotref = 'N/A'
    modpot = 'N/A'
    modpotref = 'N/A'
    delg = "N/A"
    delgnf = "N/A"
    finres = "N/A"

    try:
        os.stat(file.location + 'redox-' + file.stripDotPDB(file.filename) + '-redpot.txt')
    except:
        calc_final = False
    else:
        print_result = True
        fp = open(file.location + 'redox-' + file.stripDotPDB(file.filename) + '-redpot.txt', 'r')
        try:
            redpot = float(fp.readline())
        except:
            pass
        fp.close()
    try:
        os.stat(file.location + 'redox-' + file.stripDotPDB(file.filename) + '-redpotref.txt')
    except: 
        calc_final = False
    else:
        print_result = True
        fp = open(file.location + 'redox-' + file.stripDotPDB(file.filename) + '-redpotref.txt', 'r')
        try:
            redpotref = float(fp.readline())
        except:
            pass
        fp.close()
    try:
        os.stat(file.location + 'redox-' + file.stripDotPDB(file.filename) + '-modpot.txt')
    except: 
        calc_final = False
    else:
        print_result = True
        fp = open(file.location + 'redox-' + file.stripDotPDB(file.filename) + '-modpot.txt', 'r')
        try:
            modpot = float(fp.readline())
        except:
            pass
        fp.close()
    try:
        os.stat(file.location + 'redox-' + file.stripDotPDB(file.filename) + '-modpotref.txt')
    except:
        calc_final = False
    else:
        print_result = True
        fp = open(file.location + 'redox-' + file.stripDotPDB(file.filename) + '-modpotref.txt', 'r')
        try:
            modpotref = float(fp.readline())
        except:
            pass
        fp.close()

    if calc_final:
        # figure out the correct ade value
        rp = apbs.models.redoxParams.objects.filter(pdb = file, selected = 'y')[0]
        ade = 0.0
        delg = 0.0
        if rp.redoxsite == 'couple_oxi':
            delg = redpot - redpotref - modpot + modpotref
            ade = 0.273
        elif rp.redoxsite == 'couple_red':
            delg = modpot - modpotref - redpot + redpotref
            ade = -3.543

        delgnf = delg * (-4.184/96.485)
        finres = delgnf - 4.43 + ade 
    
    return django.shortcuts.render_to_response('html/redox.html', {'filename_list': filename_list, 'redox_list': redox_list, \
                                                                   'noredox': noredox, 'redox_segs': redox_segs, 'print_result': print_result, \
                                                                   'redpot': redpot, 'redpotref': redpotref, 'modpot': modpot, \
                                                                   'modpotref': modpotref, 'delg': delg, 'delgnf': delgnf, 'finres': finres, \
                                                                   'redox_nums': redox_nums })

def genstruct_tpl(request,file,scriptlist):
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

def gengrid_tpl(request,file,scriptlist):
    td = {}
    try:
        td['srad'] = request.POST['srad']
    except:
        td['srad'] = 1.4

    # step 1: get the grid for the full structure
    td['filebase'] = file.stripDotPDB(file.filename)
    td['input_file'] = "reduced_" + file.stripDotPDB(file.filename) + "-final"
    td['grid_name'] = "pro"
    td['data_home'] = charmming_config.data_home
    try:
        td['pdie'] = float(request.POST['prot_diel'])
        td['sdie'] = float(request.POST['solv_diel'])
    except:
        td['pdie'] = 4
        td['sdie'] = 78
    td['sizegrid'] = True

    inp_filename = "redox-" + file.stripDotPDB(file.filename) + "-mkgrid1.inp"
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_mkgrid.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(file.location + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    scriptlist.append(file.location + inp_filename)

    # step 2: get the grid for only the cluster site
    td['sizegrid'] = False
    td['input_file'] = "reduced_" + file.stripDotPDB(file.filename) + "-reference"
    td['grid_name'] = "ref"
    try:
        td['pdie'] = float(request.POST['redx_diel'])
        td['sdie'] = float(request.POST['prot_diel'])
    except:
        td['pdie'] = 1
        td['sdie'] = 4

    inp_filename = "redox-" + file.stripDotPDB(file.filename) + "-mkgrid2.inp"
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_mkgrid.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(file.location + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    scriptlist.append(file.location + inp_filename)


def modstruct_tpl(request,file,scriptlist,redox_segs):
    td = {}
    changes = {}
    patches = {}
    newtype = {}
    redox_nums = {}
    retarr = [] # returns the redox sites as a flat list

    for elt in file.fes4list.split(','):
        segment,numfes = elt.split(':')
        segid  = segment.split('-')[0]
        numfes = int(numfes)
        redox_nums[segid] = numfes


    td['redox_segs'] = redox_segs
    for segid in redox_segs:
        changes[segid] = [] 
        for sitenum in range(1,redox_nums[segid]+1):
            kval = "site_%s_%s" % (segid.upper(),sitenum)
            if request.POST['picksite'] == kval:
                if request.POST['couple'] == 'oxi':
                    changes[segid].append(1)
                    retarr.append(1)
                elif request.POST['couple'] == 'red':
                    changes[segid].append(2)
                    retarr.append(2)
                else:
                    raise "need couple_oxi or couple_red!!!"
            else:
                changes[segid].append(0)
                retarr.append(0)

    sitestochange = []
    for segid in redox_segs:
        for sitenum in range(redox_nums[segid]):
            tdict = {}
            tdict['segid'] = segid
            tdict['sitenum'] = sitenum + 1
            if changes[segid][sitenum] == 0:
                tdict['newres'] = '4FSR'
            elif changes[segid][sitenum] == 1:
                tdict['newres'] = '4FSO'
                tdict['FEAR'] = 226
                tdict['FEBR'] = 227
                tdict['SAR']  = 232
                tdict['SBR']  = 233
                tdict['SR']   = 236
                # update charges
                tdict['CFE1'] = 0.3040
                tdict['CFE2'] = 0.3040
                tdict['CFE3'] = 0.3210
                tdict['CFE4'] = 0.3210
                tdict['CS1'] = -0.2330 
                tdict['CS2'] = -0.2330 
                tdict['CS3'] = -0.2380
                tdict['CS4'] = -0.2380
                tdict['CSG1'] = -0.3980
                tdict['CCB1'] = -0.0935
                tdict['CSG2'] = -0.3980
                tdict['CCB2'] = -0.0935
                tdict['CSG3'] = -0.4290
                tdict['CCB3'] = -0.0935
                tdict['CSG4'] = -0.4290
                tdict['CCB4'] = -0.0935
            elif changes[segid][sitenum] == 2:
                tdict['newres'] = '4FSS'
                tdict['FEAR'] = 222
                tdict['FEBR'] = 223
                tdict['SAR']  = 228
                tdict['SBR']  = 229
                tdict['SR']   = 234
                # update charges
                tdict['CFE1'] = 0.4810
                tdict['CFE2'] = 0.4810
                tdict['CFE3'] = 0.4840
                tdict['CFE4'] = 0.4840
                tdict['CS1'] = -0.5230
                tdict['CS2'] = -0.5230
                tdict['CS3'] = -0.5420
                tdict['CS4'] = -0.5420
                tdict['CSG1'] = -0.7240
                tdict['CCB1'] = -0.1660
                tdict['CSG2'] = -0.7240
                tdict['CCB2'] = -0.1660
                tdict['CSG3'] = -0.7100
                tdict['CCB3'] = -0.1600
                tdict['CSG4'] = -0.7100
                tdict['CCB4'] = -0.1600
            sitestochange.append(tdict)
    td['sitestochange'] = sitestochange
    td['input_file'] = "reduced_" + file.stripDotPDB(file.filename) + "-final"
    td['filebase'] = file.stripDotPDB(file.filename)
    td['dxmath'] = "/usr/local/charmming/apbs-1.1.0/share/tools/mesh/dxmath"
    td['redox_segs'] = redox_segs
    td['topology_list'] = ['%s/toppar/top_all27_prot_na.rtf' % charmming_config.data_home, \
                           '%s/toppar/top_4fsr.rtf' % charmming_config.data_home]
    td['parameter_list'] = ['%s/toppar/par_all27_prot_na.prm' % charmming_config.data_home,
                            '%s/toppar/par_4fsr.prm' % charmming_config.data_home]

    # create inputs to dxmath
    try:
        dielectric = float(request.POST['prot_diel'])
    except:
        dielectric = 4.0
    td['dxmath_inp'] = []
    for axis in ["x", "y", "z"]:
        dxmath_scr = "combo.%s" % axis
        fp = open(file.location + dxmath_scr, 'w')
        fp.write("# combine mesh for %s\n" % axis)
        fp.write("%s/ref.%s %4.3f - %s/pro.%s + %s/iapbs-diel%s.dx =\n" % (file.location,axis,dielectric,file.location,axis,file.location,axis))
        fp.close()
        td['dxmath_inp'].append(dxmath_scr)

    inp_filename = "redox-" + file.stripDotPDB(file.filename) + "-modstruct.inp"
    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_modstruct.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(file.location + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()

    scriptlist.append(file.location + inp_filename)
    return retarr

def getdelg_tpl(request,file,scriptlist,redox_segs):
    td = {}
    try:
        td['srad'] = float(request.POST['srad'])
    except:
        td['srad'] = 1.4
   
    td['topology_list'] = ['%s/toppar/top_all27_prot_na.rtf' % charmming_config.data_home, \
                           '%s/toppar/top_4fsr.rtf' % charmming_config.data_home]
    td['parameter_list'] = ['%s/toppar/par_all27_prot_na.prm' % charmming_config.data_home,
                            '%s/toppar/par_4fsr.prm' % charmming_config.data_home]
    td['filebase'] = file.stripDotPDB(file.filename)
    td['data_home'] = charmming_config.data_home
    td['redox_segs'] = redox_segs
    td['do_resdel'] = False
    td['resdel_list'] = []

    redox_nums = {}

    for elt in file.fes4list.split(','):
        segment,numfes = elt.split(':')
        segid  = segment.split('-')[0]
        numfes = int(numfes) 
        redox_nums[segid] = numfes

    changes = {}
    for segid in redox_segs:
        changes[segid] = []
        for sitenum in range(1,redox_nums[segid]+1):
            kval = "site_%s_%s" % (segid.upper(),sitenum)
            if request.POST['picksite'] == kval:
                if request.POST['couple'] == 'oxi':
                    changes[segid].append(1)
                elif request.POST['couple'] == 'red':
                    changes[segid].append(2)
                else:
                    raise "need couple_oxi or couple_red!!!"
            else:
                td['resdel_list'].append(sitenum)
                changes[segid].append(0)

    # for each of redpot, redpotref, modpot, and modpotref generate the delta G via APBS

    # redpot
    redox_seglist = []
    for segid in redox_segs:
        tdict = {}
        tdict['segid'] = segid
        tdict['resname'] = 'resname 4FSR' # reduced structure
        redox_seglist.append(tdict)
    td['redox_seglist'] = redox_seglist

    td['redox_selectall'] = False
    td['resid'] = '4FSR'
    td['input_file'] = 'reduced_' + file.stripDotPDB(file.filename) + '-final'
    td['rdiel'] = 'rdiel'
    td['pdie'] = 4
    td['sdie'] = 78
    td['enefile'] = 'redox-' + file.stripDotPDB(file.filename) + '-redpot.txt'
    inp_filename = 'redox-' + file.stripDotPDB(file.filename) + '-redpot.inp'

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_delg.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(file.location + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    scriptlist.append(file.location + inp_filename)

    # redpotref
    redox_seglist = []
    for segid in redox_segs:
        for i in range(len(changes[segid])):
            if changes[segid][i] != 0:
                tdict = {}
                tdict['segid'] = segid
                tdict['resname'] = 'resname 4FSR'
                tdict['resname'] += ' .and. resid %d' % (i+1)
                redox_seglist.append(tdict)
    td['redox_seglist'] = redox_seglist
    if len(td['resdel_list']) > 0:
        td['do_resdel'] = True

    td['redox_selectall'] = False
    td['input_file'] = 'reduced_' + file.stripDotPDB(file.filename) + '-reference'
    td['rdiel'] = ''
    td['pdie'] = 1
    td['sdie'] = 1
    td['enefile'] = 'redox-' + file.stripDotPDB(file.filename) + '-redpotref.txt'
    inp_filename = 'redox-' + file.stripDotPDB(file.filename) + '-redpotref.inp'

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_delg.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(file.location + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    scriptlist.append(file.location + inp_filename)

    # figure out new resids
    redox_seglist = []

    for segid in redox_segs:
        tdict = {}
        tdict['segid'] = segid
        tdict['resname'] = ''

        myiter = 0
        for codenum in changes[segid]:
            myiter += 1
            if codenum == 1:
                if len(tdict['resname']) > 0:
                    tdict['resname'] += ' .or. ( resname 4FSO .and. resid %d )' % myiter
                else:
                    tdict['resname'] += '( resname 4FSO .and. resid %d )' % myiter
            elif codenum == 2:
                if len(tdict['resname']) > 0:
                    tdict['resname'] += ' .or. ( resname 4FSS .and. resid %d )' % myiter
                else:
                    tdict['resname'] += '( resname 4FSS .and. resid %d )' % myiter
        redox_seglist.append(tdict)

    td['redox_seglist'] = redox_seglist

    # modpot
    td['do_resdel'] = False
    td['redox_selectall'] = True
    td['input_file'] = 'modified_' + file.stripDotPDB(file.filename) + '-final'
    td['rdiel'] = 'rdiel'
    td['pdie'] = 4
    td['sdie'] = 78
    td['enefile'] = 'redox-' + file.stripDotPDB(file.filename) + '-modpot.txt'
    inp_filename = 'redox-' + file.stripDotPDB(file.filename) + '-modpot.inp'

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_delg.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(file.location + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    scriptlist.append(file.location + inp_filename)

    # modpotref
    td['redox_selectall'] = False
    td['input_file'] = 'modified_' + file.stripDotPDB(file.filename) + '-reference'
    td['rdiel'] = ''
    td['pdie'] = 1
    td['sdie'] = 1
    if len(td['resdel_list']) > 0:
        td['do_resdel'] = True

    td['enefile'] = 'redox-' + file.stripDotPDB(file.filename) + '-modpotref.txt'
    inp_filename = 'redox-' + file.stripDotPDB(file.filename) + '-modpotref.inp'

    t = django.template.loader.get_template('%s/mytemplates/input_scripts/redox_delg.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(django.template.Context(td)))
    fp = open(file.location + inp_filename, 'w')
    fp.write(charmm_inp)
    fp.close()
    scriptlist.append(file.location + inp_filename)

def redox_tpl(request,file,scriptlist):

    try:
        oldparam = apbs.models.redoxParams.objects.filter(pdb = file, selected = 'y')[0]
        oldparam.selected = 'n'
        oldparam.save()
    except:
        pass

    # preliminaries, set up the model
    rp = apbs.models.redoxParams(selected='y')
    rp.pdb = file
    rp.statusHTML = "<font color=yellow>Processing</font>"
    rp.save()

    # step 1: make reduced CRD & PSF and oxidized/superreduced CRD & PSF for protein and reference
    redox_segs = genstruct_tpl(request,file,scriptlist)

    # step 2: make dielectric grids (1. protein + SF4 2. just SF4 + hanging -CH2)
    # This script will also include a system call to dxmath to make the grids
    gengrid_tpl(request,file,scriptlist)

    # step 3: Get combined grids, then change reduced to oxidized or superreduced structure
    redox_sites = modstruct_tpl(request,file,scriptlist,redox_segs)
    rp.segments = "%s" % redox_segs
    rp.redoxsite = ""
    if request.POST['couple'] == 'oxi':
        rp.redoxsite = 'couple_oxi'
    elif request.POST['couple'] == 'red':
        rp.redoxsite = 'couple_red' 
    rp.save()

    # step 4: Get the four free energy values for redpot, redpotref, modpot, modpotref where mod = oxi or sr
    getdelg_tpl(request,file,scriptlist,redox_segs)

    # all scripts generated, submit to the scheduler
    user_id = file.owner.id
    si = scheduler.schedInterface.schedInterface()
    newJobID = si.submitJob(user_id,file.location,scriptlist)
    # NB: we'll need to add lesson hooks here
    #if file.lesson_type:
    #    lessonaux.doLessonAct(file,"onMinimizeSubmit",postdata,final_pdb_name)

    if newJobID < 0:
       rp.statusHTML = "<font color=red>Failed</font>"
       rp.save()
    else:
       file.redox_jobID = newJobID
       sstring = si.checkStatus(newJobID)
       rp.statusHTML = scheduler.statsDisplay.statsDisplay(sstring,newJobID)
       rp.save()
       file.save()


    return "foo"
