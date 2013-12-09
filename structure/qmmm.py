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

import copy
import re
import selection.models
import charmming_config

# Yee-haw, a recursive function
def writeQMheader(charmm_inp, selection):
    tlen = len(selection)

    if tlen > 1500:
      # we're not even gonna check something this long
      return charmm_inp
    if tlen <= 65:
      charmm_inp += "* QM REGION: %s\n" % selection.upper()
      return charmm_inp
    else:
      stoppt = 64
      while stoppt >= 0:
         tchr = selection[stoppt]
         if tchr == ' ':
            break
         stoppt -= 1
      if stoppt == 0:
         charmm_inp += "* Error, select token too long!\n"
         return charmm_inp
      else:
         charmm_inp += "* QM REGION: %s\n" % selection[:stoppt].upper()
         return writeQMheader(charmm_inp, selection[stoppt+1:])


# This function returns 3 possible values:
# 0: QM region has changed
# 1: QM region is the same
# 2: No previous QM region to check
def checkQMregion(ws, selection):
    logfp = open('/tmp/checkQMregion.txt', 'w')

    qmsele = ws.qmRegion.lower()

    if qmsele == 'none':
        logfp.write('return 2\n')
        logfp.close()
        return 2
    if qmsele != selection.lower():
        logfp.write('return 0\n')
        logfp.close()
        return 0

    logfp.write('return 1\n')
    logfp.close()
    return 1

def wrtQchemInp(fname, exchange, correlation, basis, jobtype, charge, multi, fullhess = False):
    logfp = open('/tmp/wrtqchem.txt', 'w')
    logfp.write('writing qchem input file to %s\n' % fname)
    logfp.close()

    q=open(fname, 'w')

    print >>q, '$COMMENT'
    print >>q, 'Q-Chem control file written by CHARMMing.org'
    print >>q, '$END'
    print >>q, ' '

    print >>q, '$REM'
    if exchange == 'B3' and correlation == 'LYP':
       print >>q, 'EXCHANGE              ',exchange+correlation
    elif exchange == 'HF':
       print >>q, 'EXCHANGE              ',exchange
    else:
       print >>q, 'EXCHANGE              ',exchange
       print >>q, 'CORRELATION           ',correlation

    print >>q, 'BASIS                 ',basis
    print >>q, 'JOBTYPE               ',jobtype
    print >>q, 'QM_MM                  TRUE'
    print >>q, 'SYMMETRY               OFF'
    print >>q, 'SYM_IGNORE             TRUE'
    print >>q, 'PRINT_INPUT            FALSE'
    print >>q, 'QMMM_PRINT             TRUE'
    if exchange == 'B' or exchange == 'B3':
       print >>q, 'XC_GRID                1'

    if jobtype == 'freq' and fullhess:
       print >>q, 'QMMM_FULL_HESSIAN      TRUE'
    print >>q, '$END'
    print >>q, ' '

    print >>q, '$MOLECULE'
    print >>q, charge, multi
    print >>q, '$END'
    print >>q, ' '
    q.close()

#Note: This is made to take two args so as to make it easy to process.
#Since we may implement more models than ONIOM and QM/MM, it's better to leave modelType as a field
#Instead of using properties of the database object (especially since the exceptions can get weird).
#For ONIOM, just loop over this thing multiple times.
def makeQchem_val(modelType,atomselection,**kwargs):
    prms = {}

    logfp = open('/tmp/makeqchem_val.txt', 'w')
    logfp.write('In makeqchem_val\n')
#   Since we already save the POST data to the database in objects it makes no sense
#   To re-process it here.
#   This code is left here for historical purposes.

    #TODO: Change all instances of this method to use modelType!
    if modelType == "qmmm":
        #TODO: pre-validate this data in atomselection_aux.validateInputs!
        prms['exch'] = atomselection.exchange
        prms['corr'] = atomselection.correlation
        prms['bs'] = atomselection.basis_set
        ##prms['qmsel'] = atomselection.selectionstring
        prms['charge'] = atomselection.charge
        prms['multi'] = atomselection.multiplicity
        if atomselection.linkatom_num > 0:
            prms['linkat_list'] = handleLinkAtoms(atomselection)
        else:
            prms['linkat_list'] = None
    elif modelType == "oniom":
        if atomselection.isQM:
            #If not, just blank this out, but we should check elsewhere. Make the check anyway.
            prms['exch'] = atomselection.exchange
            prms['corr'] = atomselection.correlation
            prms['bs'] = atomselection.basis_set
            ##prms['qmsel'] = atomselection.selectionstring
            prms['charge'] = atomselection.charge
            prms['multi'] = atomselection.multiplicity
            if atomselection.linkatom_num > 0:
                prms['linkat_list'] = handleLinkAtoms(atomselection)
            else:
                prms['linkat_list'] = None

    prms['bynum_list'] = []
    prms['region_list'] = []
    selar = atomselection.selectionstring.replace('.or.',"").split('bynum')
    logfp.write('selar = %s\n' % selar)
    for bnp in selar:
        bnp = bnp.strip()
        if bnp:
            bynumdict = {}
            # ToDo -- need some error checking in here
            logfp.write('bnp = %s\n' % bnp)
            bynumdict['numa'] = bnp.split(':')[0]
            try:
                bynumdict['numb'] = bnp.split(':')[1]
            except IndexError, e:
                bynumdict['numb'] = -9999
            if modelType == "qmmm":
                prms['region_list'].append(1)
            elif modelType == 'oniom':
                pass
            prms['bynum_list'].append(bynumdict)

    logfp.close()
    return prms

##def makeQChem_tpl(template_dict, file, exch, corr, bs, qmsel, jobtype, charge, multi, last_struct, linkatoms = None):
#Fname is included because we want to be able to call wrtQchemInp with it.
#Task is included because we need a decent naming convention. 
#TODO: Adapt this to multi-task file name conventions
def makeQChem_tpl(template_dict,qmparms,fname,task):
    template_dict['bad_exchange'] = 0
    template_dict['bad_correlation'] = 0
    template_dict['bad_basis'] = 0
    template_dict['bad_jobtype'] = 0
    if qmparms['exch'] != "HF" and qmparms['exch'] != "B" and qmparms['exch'] != "B3LYP":
      template_dict['bad_exchange'] = 1
      return template_dict
    if qmparms['corr'] != "None" and qmparms['corr'].upper() != "LYP":
      template_dict['bad_correlation'] = 1
      return template_dict
    if qmparms['bs'] != "STO-3G" and qmparms['bs'] != "3-21G*" and qmparms['bs'] != "6-31G*":
      template_dict['bad_basis'] = 1
      return template_dict
    if qmparms['jobtype'] != "Force" and qmparms['jobtype'] != "SP" and qmparms['jobtype'] != "freq":
      template_dict['bad_jobtype'] = 1
      return template_dict
    if not fname: #In the interests of safety, make fname start off as False. TODO: Ask Tim about kwargs.
        qchemin  = "qchem-" + task.workstruct.identifier + "-" + task.action + ".in"
        qcheminp = "qchem-" + task.workstruct.identifier + "-" + task.action + ".inp"
        qchemout = "qchem-" + task.workstruct.identifier + "-" + task.action + ".out"
    else:
        qchemin = fname + ".in"
        qcheminp = fname + ".inp"
        qchemout = fname + ".out"
    ## BTM -- I vaguely remember why we put this in here, but I am not sure how it fits in
    ## with the new GUTS-style CHARMMing, so I am going to leave it commented out for now.
    ##qmcheck = checkQMregion(task.workstruct,qmparms['qmsel'])

    template_dict['qmcheck'] = False
    ##if qmcheck == 0:
    ##    return template_dict

    fullhess = False
    template_dict['jobtype'] = qmparms['jobtype']
    if template_dict['jobtype'] == 'freq':
       tqm = qmparms['qmsel'].strip()
       tqm = tqm.lower()
       template_dict['tqm'] = tqm
       if tqm != 'all':
           fullhess = True
    if not fname:
        wrtQchemInp(task.workstruct.structure.location + '/' + qchemin, qmparms['exch'], qmparms['corr'], \
               qmparms['bs'], qmparms['jobtype'], qmparms['charge'], qmparms['multi'], fullhess) #This will write to a specific filename, which helps with ONIOM.
    else:
        wrtQchemInp(task.workstruct.structure.location + '/' + qchemin, qmparms['exch'], qmparms['corr'], \
               qmparms['bs'], qmparms['jobtype'], qmparms['charge'], qmparms['multi'], fullhess) #This will write to a specific filename, which helps with ONIOM.
    ##template_dict['qmsel'] = qmparms['qmsel']

    template_dict['region_list'] = copy.deepcopy(qmparms['region_list'])
    template_dict['bynum_list'] = copy.deepcopy(qmparms['bynum_list'])
    template_dict['linkat_list'] = copy.deepcopy(qmparms['linkat_list'])

    template_dict['qchemin'] = qchemin
    template_dict['qcheminp'] = qcheminp
    template_dict['qchemout'] = qchemout
    template_dict['data_home'] = charmming_config.data_home

    return template_dict


# pre:requires postdata that has some patch information
#file is not used. removing from the method...
#def handleLinkAtoms(file,postdata):
#Doesn't need postdata either. Just an AtomSelection. For OniomSelections, just loop.
def handleLinkAtoms(atomselection):

    #used in disulfide bond patching
    linkqmsegid = ''
    linkmmsegid = ''
    linkqm = ''
    qmatomtype = ''
    qmatomtype = ''
    linkmm = ''
    link_list = []
    link_tuple = ()

    lonepairs = selection.models.LonePair.objects.filter(selection=atomselection) #All atomselections have that in them, and it's an int field so there's no reason to worry.

    for lonepair in lonepairs:
        ldict = {} 
        ldict['segq'] = lonepair.qmsegid
        ldict['segm'] = lonepair.mmsegid
        ldict['resq'] = lonepair.qmresid
        ldict['resm'] = lonepair.mmresid
        ldict['atypq'] = lonepair.qmatomname
        ldict['atypm'] = lonepair.mmatomname

        link_list.append(ldict)

    return link_list

