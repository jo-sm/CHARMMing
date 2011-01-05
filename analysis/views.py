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

from django.template.loader import get_template
from django.http import HttpResponse
from django.shortcuts import render_to_response
from structure.models import Structure
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from django.template import *
import output, charmming_config, lessonaux
import os, re, time

def parseRMSD(file,mlist):
    time.sleep(2)

    rmsds = []
    for i in range(len(mlist)):
        rmsds.append({})

    rmsdsig = re.compile("! RMS ([a-z]+) and ([a-z]+)")
    rmsdval = re.compile("THUS RMS DIFF IS[ ]+([0-9]+\.[0-9]+)")

    rfd = open(file.location + 'charmm-' + file.stripDotPDB(file.filename) + '-rmsd.out', 'r')
    cidx = -1
    comp = None
    for line in rfd:
        m = rmsdsig.search(line)
        if m:
            base = "-" + m.group(1)
            comp = m.group(2)
            cidx = -1
            for i in range(len(mlist)):
                if base in mlist[i]:
                    cidx = i
                    break
            if cidx == -1:
                # To-Do: handle this error more gracefully
                raise "parseRMSD couldn't find %s in %s." % (base,mlist)
        m = rmsdval.search(line)
        if m:
            if comp and cidx >= 0:
                rmsds[cidx][comp] = m.group(1)
            else:
                raise "parseRMSD had invalid RMS info."
            comp = None
            cidx = -1
    rfd.close()

    fbase = re.compile("-([0-9]+)-([0-9]+)-(.*?)\.pdb")
    rmsdhtml = '<table cellspacing="3" cellpadding="2">\n'
    rmsdhtml += '<tr><td></td>'
    fnbases = []
    for i in range(len(mlist)):
        m = fbase.search(mlist[i])
        fname = m.group(3)
        if fname == 'final':
            fname = 'append'
        rmsdhtml += '<td>%s</td>'% fname
        fnbases.append(fname)
    rmsdhtml += '</tr>\n'
    for i in range(len(mlist)):
        m = fbase.search(mlist[i])
        fname = m.group(3)
        if fname == 'final':
            fname = 'append'
        rmsdhtml += '<tr><td>%s</td>' % fname
        for cbase in fnbases:
            if cbase == fname:
                rmsdhtml += '<td>0.000000</td>'
            elif rmsds[i].has_key(cbase):
                rmsdhtml += '<td>%s</td>' % rmsds[i][cbase]
            else:
                rmsdhtml += '<td></td>'
        rmsdhtml += '</tr>\n'
    rmsdhtml += '</table>\n'

    ofp = open(file.location + "/rmsd-" + file.stripDotPDB(file.filename) + ".html", "w")
    ofp.write(rmsdhtml)
    ofp.close()

    return render_to_response('html/displayrmsd.html', {'rmsdhtml': rmsdhtml})

def calcRMSD_tpl(file,postdata,mlist):
    fbase = re.compile("-([0-9]+)-([0-9]+)-(.*?)\.pdb")
    done = re.compile('Done')
    failed = re.compile('Failed')
    # template dictionary passes the needed variables to the template
    template_dict = {}

    if postdata.has_key('comp_allatom'):
        try:
            cmplev = int(postdata['comp_allatom'])
        except:
            cmplev = 1
    else:
        cmplev = 1

    if cmplev == 0:
        selpart = ""
        if postdata.has_key('bb_n'):
            selpart += 'type N'
        if postdata.has_key('bb_ca'):
            if selpart:
                selpart += ' .or. '
            selpart += 'type CA'
        if postdata.has_key('bb_c'):
            if selpart:
                selpart += ' .or. '
            selpart += 'type C'
        if postdata.has_key('bb_o'):
            if selpart:
                selpart += ' .or. '
            selpart += 'type O'
        if postdata.has_key('bb_hn'):
            if selpart:
                selpart += ' .or. '
            selpart += 'type HN'
        if postdata.has_key('bb_ha'):
            if selpart:
                selpart += ' .or. '
            selpart += 'type HA'
        template_dict['selpart'] = selpart
  
    os.chdir(file.location)
    template_dict['topology_list'] = file.getTopologyList()
    template_dict['parameter_list'] = file.getParameterList()
    template_dict['filebase'] = file.stripDotPDB(file.filename)
    template_dict['cmplev'] = cmplev
    # since template limits using certain python functions, here a list of dictionaries is used to pass mult-variables
    template_dict['dict_lis'] = []

    for i in range(len(mlist)-1):
        m = fbase.search(mlist[i])
        fname = m.group(3)
        dict = {}
        # another list of dictionaries is used to pass mult-variables 
        dict['lis'] = []
        dict['base'] = file.stripDotPDB(mlist[i]) 
        dict['fname'] = fname

        for j in range(i+1,len(mlist)):
            n = fbase.search(mlist[j])
            cname = n.group(3)

            # if we're using a sovlated structure, ignore the bulk solvent.
            selectstr = ''
            if "solv" in cname or "neut" in cname:
                selectstr += "select .not. segid bwat "
            if cmplev == 0:
                if selectstr:
                    selectstr += ".and. ( %s ) " % selpart
                elif selpart:
                    selectstr = "select " + selpart
            if selectstr:
                selectstr += " end"
            sub_dict = {}
            sub_dict['base'] = file.stripDotPDB(mlist[j])
            sub_dict['selectstr'] = selectstr
            sub_dict['cname'] = cname
            dict['lis'].append(sub_dict)   
        template_dict['dict_lis'].append(dict)
    
    t = get_template('%s/mytemplates/input_scripts/calcRMSD_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))

    rmsd_inp = file.location + 'charmm-' + file.stripDotPDB(file.filename) + '-rmsd.inp'
    inp_out = open(rmsd_inp,'w')
    inp_out.write(charmm_inp)
    inp_out.close()

    si = schedInterface()
    rmsdJobID = si.submitJob(file.owner.id,file.location,[rmsd_inp])
    sstring = si.checkStatus(rmsdJobID)
    while(not done.search(statsDisplay(sstring,rmsdJobID)) and not failed.search(statsDisplay(sstring,rmsdJobID))):
        sstring = si.checkStatus(rmsdJobID)
    file.save()
    if failed.search(statsDisplay(sstring,rmsdJobID)):
        return render_to_response('html/displayrmsd.html', {'rmsdhtml': '<b>RMSD calculation failed.</b>'})
    return parseRMSD(file,mlist)

def rmsformdisplay(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    Structure.checkRequestData(request)
    #chooses the file based on if it is selected or not
    try:
        file =  Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")
    os.chdir(file.location)

    feedback = ''

    #creates a list of filenames associated with the PDB
    #The "md" option takes away all -md,-ld, and -sgld names
    #filename_list = file.getLimitedFileList("md")
    #filename_list = file.getLimitedFileList('blank')
    append_pdb = "new_" + file.stripDotPDB(file.filename) + "-final.pdb"
    solv_pdb = "new_" + file.stripDotPDB(file.filename) + "-solv.pdb"
    neu_pdb = "new_" + file.stripDotPDB(file.filename) + "-neutralized.pdb"
    min_pdb = "new_" + file.stripDotPDB(file.filename) + "-min.pdb"
    md_pdb = "new_" + file.stripDotPDB(file.filename) + "-md.pdb"
    md_avg_pdb = "new_" + file.stripDotPDB(file.filename) + "-mdavg.pdb"
    ld_pdb = "new_" + file.stripDotPDB(file.filename) + "-ld.pdb"
    sgld_pdb = "new_" + file.stripDotPDB(file.filename) + "-sgld.pdb"

    filename_list = []

    # make sure the appended PDB exists.
    try:
        os.stat(file.location + append_pdb)
        filename_list.append(append_pdb) # This is not in the limited file list
    except:
        return render_to_response('html/rmsdnoappend.html')
    try:
        os.stat(file.location + solv_pdb)
        filename_list.append(solv_pdb) # This is not in the limited file list
    except:
        pass
    try:
        os.stat(file.location + neu_pdb)
        filename_list.append(neu_pdb) # This is not in the limited file list
    except:
        pass
    try:
        os.stat(file.location + min_pdb)
        filename_list.append(min_pdb) # This is not in the limited file list
    except:
        pass
    try:
        os.stat(file.location + md_pdb)
        filename_list.append(md_pdb) # This is not in the limited file list
    except:
        pass
    try:
        os.stat(file.location + md_avg_pdb)
        filename_list.append(md_avg_pdb) # This is not in the limited file list
    except:
        pass
    try:
        os.stat(file.location + ld_pdb)
        filename_list.append(ld_pdb) # This is not in the limited file list
    except:
        pass
    try:
        os.stat(file.location + sgld_pdb)
        filename_list.append(sgld_pdb) # This is not in the limited file list
    except:
        pass

    # build a list of the shizzle fo' da matrix
    matrix_list = []
    if request.POST.has_key(append_pdb):
        matrix_list.append(append_pdb)
    if request.POST.has_key(solv_pdb):
        matrix_list.append(solv_pdb)
    if request.POST.has_key(neu_pdb):
        matrix_list.append(neu_pdb)
    if request.POST.has_key(min_pdb):
        matrix_list.append(min_pdb)
    if request.POST.has_key(md_pdb):
        matrix_list.append(md_pdb)
    if request.POST.has_key(md_avg_pdb):
        matrix_list.append(md_avg_pdb)
    if request.POST.has_key(ld_pdb):
        matrix_list.append(ld_pdb)
    if request.POST.has_key(sgld_pdb):
        matrix_list.append(sgld_pdb)

    if len(matrix_list) == 1:
        feedback = 'You must select at least two structures.'
    elif len(matrix_list) > 1:
        if file.lesson_type:
            lessonaux.doLessonAct(file,"onRMSDSubmit",request.POST)
        file.save()
        return calcRMSD_tpl(file,request.POST,matrix_list)
#        return calcRMSD(file,request.POST,matrix_list)

    prior_matrix = ''
    try:
        os.stat(file.location + '/rmsd-' + file.stripDotPDB(file.filename) + '.html')
    except:
        pass
    else:
        fp = open(file.location + '/rmsd-' + file.stripDotPDB(file.filename) + '.html')
        for line in fp:
            prior_matrix += line
        fp.close() 

    return render_to_response('html/rmsdform.html', {'filename_list': filename_list, 'feedback': feedback, 'prior_matrix': prior_matrix, 'append_pdb': append_pdb, \
                              'solv_pdb': solv_pdb, 'neu_pdb': neu_pdb, 'min_pdb': min_pdb, 'md_pdb': md_pdb, 'ld_pdb': ld_pdb, 'sgld_pdb': sgld_pdb, \
                              'md_avg_pdb': md_avg_pdb })

import dynamics.models

def domdprop(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    Structure.checkRequestData(request)
    #chooses the file based on if it is selected or not
    try:
        file =  Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")

    os.chdir(file.location)
    fname = "charmm-%s-md.out" % file.stripDotPDB(file.filename)
    cmdline = "%s/prop.py %s %s-mdproperties.dat" % (charmming_config.data_home,fname,file.stripDotPDB(file.filename))
    
    nprop = 0
    if request.POST.has_key('getener'):
        cmdline += " averener"
    if request.POST.has_key('getvolu'):
        cmdline += " avervolu"
    if request.POST.has_key('getpressi'):
        cmdline += " averpressi"
    if request.POST.has_key('gettemp'):
        cmdline += " avertemp"
    if request.POST.has_key('gettote'):
        cmdline += " avertote"
    if request.POST.has_key('gettotk'):
        cmdline += " avertotk"

    logfp = open('/tmp/mdprop.txt', 'w')
    logfp.write('cmdline = %s\n' % cmdline)
    logfp.close()
    os.system(cmdline)
    return HttpResponse("Properties generated, please download %s-mdproperties.dat from the download files page." % (file.stripDotPDB(file.filename)))

def getmdprop(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
    Structure.checkRequestData(request)
    #chooses the file based on if it is selected or not
    try:
        file =  Structure.objects.filter(owner=request.user,selected='y')[0]
    except:
        return HttpResponse("Please submit a structure first.")

    try:
        mdp = dynamics.models.mdParams.objects.filter(pdb=file, selected='y')[0]
    except:
        return HttpResponse("Please perform Molecular Dynamics before trying to get properties")
    if not "Done" in mdp.statusHTML:
        return HttpResponse("Please wait until Molecular Dynamics finished before trying to get properties")

    return render_to_response('html/mdanalysis.html')

