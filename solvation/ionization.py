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
import os
from pdbinfo.models import PDBFile, PDBFileForm
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from solvation.models import solvationParams
from django.template.loader import get_template
from django.template import *
import output, lessonaux, charmming_config

# This function constructs an input script with template
def neutralize_tpl(file,postdata,scriptlist):
    sp = solvationParams.objects.filter(pdb=file,selected='y')[0]
    cation = "POT"
    try:
        if postdata['salt'] == "nacl":
            cation = "SOD"
        elif postdata['salt'] == "cacl2":
            cation = "CAL"
        elif postdata['salt'] == "mgcl2":
            cation = "MG"
        elif postdata['salt'] == "cscl":
            cation = "CES" 
    except:
        pass
    sp.salt = cation
    sp.save()
    try:
        concentration = postdata['concentration']
    except:
        concentration = "0.15"
    sp.concentration = str(concentration)
    try:
        ntrials = postdata['ntrials']
    except:
        ntrials = 10
    sp.ntrials = ntrials
    sp.save()
    os.chdir(file.location)
    n_filename = "new_" + file.stripDotPDB(file.filename) + "-solv" 

    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'] = file.getTopologyList()
    template_dict['parameter_list'] = file.getParameterList()
    template_dict['filebase'] = file.stripDotPDB(file.filename)
    template_dict['n_filename'] = n_filename
    template_dict['file_location'] = file.location
    template_dict['cation'] = cation
    template_dict['concentration'] = concentration
    template_dict['ntrials'] = ntrials
    template_dict['solvation_structure'] = file.solvation_structure
	
    charmm_inp = file.makeCHARMMInputHeader('Neutralize system',postdata)
    t = get_template('%s/mytemplates/input_scripts/neutralize_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))

    neut_filename = file.location + "charmm-" + file.stripDotPDB(file.filename) + "-neutralize.inp"
    inp_out = open(neut_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()
 
    if file.solvation_structure != 'sphere':
        # set up the crystal file
        #cryst_inp = "* Crystal stream file\n*\n\n"
        relative_boundary = 0
        if float(file.crystal_x) < 0:
            relative_boundary = 1
        template_dict['relative_boundary'] = relative_boundary
        dim_x = file.crystal_x
        template_dict['dim_x'] = float(file.crystal_x)
        dim_z = file.crystal_z
        template_dict['dim_z'] = float(file.crystal_z)
        greaterval = max(dim_x,dim_z)
        template_dict['greaterval'] = str(greaterval)

        t = get_template('%s/mytemplates/input_scripts/cryst_template.inp' % charmming_config.charmming_root)
        cryst_inp = output.tidyInp(t.render(Context(template_dict)))

        cry_filename = "charmm-" + file.stripDotPDB(file.filename) + "-crystl.str"
        str_out = open(file.location + cry_filename,'w+')
        str_out.write(cryst_inp)
        str_out.close()


    # copy the addions stream file
    os.system("cp /usr/local/charmming/addions.str %s/new_%s-addions.str" % (file.location,file.stripDotPDB(file.filename)))

    # OK, let's let 'er rip...
    user_id = file.owner.id
    file.solvation_status = "<font color=yellow>Processing</font>"
    file.save()
    si = schedInterface()
    scriptlist.append(neut_filename)
    newJobID = si.submitJob(user_id,file.location,scriptlist)    
    if file.lesson_type:
        lessonaux.doLessonAct(file,"onSolvationSubmit",postdata,None)
    file.solvation_jobID = newJobID
    sstring = si.checkStatus(newJobID)
    file.solvation_status = statsDisplay(sstring,newJobID)
    file.save()
