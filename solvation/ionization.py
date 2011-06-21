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
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from solvation.models import solvationParams
from django.template.loader import get_template
from django.template import *
from shutil import copy
import output, lessonaux, charmming_config

# This function constructs an input script with template
def neutralize_tpl(workingstruct,sp,postdata,scriptlist):
    cation = "POT"
    if postdata['salt'] == "nacl":
        cation = "SOD"
    elif postdata['salt'] == "cacl2":
        cation = "CAL"
    elif postdata['salt'] == "mgcl2":
        cation = "MG"
    elif postdata['salt'] == "cscl":
        cation = "CES" 
    sp.salt = cation

    try:
        concentration = float(postdata['concentration'])
    except:
        concentration = 0.15
    sp.concentration = str(concentration)
    try:
        ntrials = int(postdata['ntrials'])
    except:
        ntrials = 3
    sp.ntrials = ntrials
    sp.save()

    os.chdir(workingstruct.structure.location)

    # template dictionary passes the needed variables to the template
    template_dict = {}
    template_dict['topology_list'] = workingstruct.getTopologyList()
    template_dict['parameter_list'] = workingstruct.getParameterList()
    template_dict['infile'] = 'solv-' + workingstruct.identifier
    template_dict['cation'] = cation
    template_dict['concentration'] = concentration
    template_dict['ntrials'] = ntrials
    template_dict['solvation_structure'] = sp.solvation_structure
    template_dict['fname'] = workingstruct.identifier
	
    t = get_template('%s/mytemplates/input_scripts/neutralize_template.inp' % charmming_config.charmming_root)
    charmm_inp = output.tidyInp(t.render(Context(template_dict)))

    neut_filename = workingstruct.structure.location + "/neutralize-" + workingstruct.identifier + ".inp"
    inp_out = open(neut_filename,'w')
    inp_out.write(charmm_inp)
    inp_out.close()
 
    if sp.solvation_structure != 'sphere':
        # set up the crystal file
        #cryst_inp = "* Crystal stream file\n*\n\n"
        t2_dict = {}
        t2_dict['shape'] = sp.solvation_structure
        t2_dict['dim_x'] = sp.xtl_x
        t2_dict['dim_y'] = sp.xtl_y
        t2_dict['dim_z'] = sp.xtl_z
        t2_dict['angles'] = "%10.6f %10.6f %10.6f" % (sp.angles[0],sp.angles[1],sp.angles[2])

        t = get_template('%s/mytemplates/input_scripts/cryst_template.inp' % charmming_config.charmming_root)
        cryst_inp = output.tidyInp(t.render(Context(t2_dict)))

        cry_filename = workingstruct.structure.location + '/' + workingstruct.identifier + "-crystl.str"
        str_out = open(cry_filename,'w+')
        str_out.write(cryst_inp)
        str_out.close()


    # copy the addions stream file
    copy("%s/addions.str" % charmming_config.data_home, "%s/%s-addions.str" % (workingstruct.structure.location,workingstruct.identifier))

    # OK, let's let 'er rip...
    user_id = workingstruct.structure.owner.id

    scriptlist.append(neut_filename)
    si = schedInterface()
    newJobID = si.submitJob(user_id,workingstruct.structure.location,scriptlist)    

    # lessons are still borked
    #if file.lesson_type:
    #    lessonaux.doLessonAct(file,"onSolvationSubmit",postdata,None)

    workingstruct.solvation_jobID = newJobID
    sstring = si.checkStatus(newJobID)
    sp.statusHTML = statsDisplay(sstring,newJobID)
    sp.save()
    workingstruct.save()
