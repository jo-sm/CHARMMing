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

import re
import selection.models

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
#Instead o fusing properties of the database object (especially since the exceptions can get weird).
#For ONIOM, just loop over this thing multiple times.
def makeQchem_val(modelType,atomselection):
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
        prms['qmsel'] = atomselection.selectionstring
        prms['charge'] = atomselection.charge
        prms['multi'] = atomselection.multiplicity
        if atomselection.linkatom_num > 0:
            prms['linkatoms'] = handleLinkAtoms(atomselection) #TODO: Modify this so it uses the database!
        else:
            prms['linkatoms'] = None
    elif modelType == "oniom":
        if atomselection.isQM:
            #If not, just blank this out, but we should check elsewhere. Make the check anyway.
            prms['exch'] = atomselection.exchange
            prms['corr'] = atomselection.correlation
            prms['bs'] = atomselection.basis_set
            prms['qmsel'] = atomselection.selectionstring
            prms['charge'] = atomselection.charge
            prms['multi'] = atomselection.multiplicity
            if atomselection.linkatom_num > 0:
                prms['linkatoms'] = handleLinkAtoms(atomselection) #TODO: Modify this so it uses the database!
            else:
                prms['linkatoms'] = None
#    if postdata['qmmm_exchange'] in ['HF','B','B3']:
#        prms['exch'] = postdata['qmmm_exchange']
#    else:   
#        prms['exch'] = 'HF'
#    if postdata['qmmm_correlation'] in ['None','LYP']:
#        prms['corr'] = postdata['qmmm_correlation']
#    else:   
#        prms['corr'] = 'None'
#    if postdata['qmmm_basisset'] in ['STO-3G','3-21G*','6-31G*']:
#        prms['bs'] = postdata['qmmm_basisset']
#    else:   
#        prms['bs'] = 'sto3g'
#    prms['qmsel'] = qmsel
#    if qmsel == '':
#        prms['qmsel'] = 'resid 1'
#    if postdata['qmmm_charge'] in ['-5','-4','-3','-2','-1','0','1','2','3','4','5']:
#        prms['charge'] = postdata['qmmm_charge']
#    else:   
#        prms['charge'] = '0'
#    if postdata['qmmm_multiplicity'] in ['0','1','2','3','4','5','6','7','8','9','10']:
#        prms['multi'] = postdata['qmmm_multiplicity']
#    else:   
#        prms['multi'] = '0'
#
#    logfp.write('num_linkatoms = %s\n' % postdata['num_linkatoms'])
#
#    if int(postdata['num_linkatoms']) > 0:
#        logfp.write('Call handleLinkAtoms\n')
#        logfp.flush()
#        prms['linkatoms'] = handleLinkAtoms(file,postdata) #file is never defined...how does this code even work?
#        logfp.write('Done\n')
#    else:   
#        prms['linkatoms'] = None

    logfp.close()
    return prms

##def makeQChem_tpl(template_dict, file, exch, corr, bs, qmsel, jobtype, charge, multi, last_struct, linkatoms = None):
#Fname is included because we want to be able to call wrtQchemInp with it.
def makeQChem_tpl(template_dict,qmparms,workstruct,fname):
    template_dict['bad_exchange'] = 0
    template_dict['bad_correlation'] = 0
    template_dict['bad_basis'] = 0
    template_dict['bad_jobtype'] = 0
    if qmparms['exch'] != "HF" and qmparms['exch'] != "B" and qmparms['exch'] != "B3":
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
        qchemin  = "qchem-" + workstruct.identifier + ".in"
        qcheminp = "qchem-" + workstruct.identifier + ".inp"
        qchemout = "qchem-" + workstruct.identifier + ".out"
    else:
        qchemin = fname + ".in"
        qcheminp = fname + ".inp"
        qchemout = fname + ".out"
    ## BTM -- I vaguely remember why we put this in here, but I am not sure how it fits in
    ## with the new GUTS-style CHARMMing, so I am going to leave it commented out for now.
    qmcheck = checkQMregion(workstruct,qmparms['qmsel'])

    template_dict['qmcheck'] = qmcheck 
    if qmcheck == 0:
      return template_dict

    fullhess = False
    template_dict['jobtype'] = qmparms['jobtype']
    if template_dict['jobtype'] == 'freq':
       tqm = qmparms['qmsel'].strip()
       tqm = tqm.lower()
       template_dict['tqm'] = tqm
       if tqm != 'all':
           fullhess = True
    if not fname:
        wrtQchemInp(workstruct.structure.location + '/' + qchemin, qmparms['exch'], qmparms['corr'], \
               qmparms['bs'], qmparms['jobtype'], qmparms['charge'], qmparms['multi'], fullhess) #This will write to a specific filename, which helps with ONIOM.
    else:
        wrtQchemInp(workstruct.structure.location + '/' + qchemin, qmparms['exch'], qmparms['corr'], \
               qmparms['bs'], qmparms['jobtype'], qmparms['charge'], qmparms['multi'], fullhess) #This will write to a specific filename, which helps with ONIOM.
    template_dict['qmsel'] = qmparms['qmsel']

    # user-specified link atoms are passed in as a list of 6-ples (qm-segid, qm-resid, qm-type,
    # mm-segid, mm-resid, mm-type)
    # If qmcheck is 1, that means that the correct link atoms are already in the structure and
    # we should not re-add them.
    template_dict['linkatom_list'] = []
    template_dict['linkatoms'] = ''
    if qmcheck != 1:
      template_dict['linkatoms'] = qmparms['linkatoms']
      if qmparms['linkatoms']:
         for las in qmparms['linkatoms']:
            try:
               qmseg, qmres, qmtyp, mmseg, mmres, mmtyp = las[0:6]
            except:
               return -1
            # TODO, validate that these are actually OK...
            atom = "%s %s %s %s %s %s" % (qmseg,qmres,qmtyp,mmseg,mmres,mmtyp)
            template_dict['linkatom_list'].append(atom)
    template_dict['qchemin'] = qchemin
    template_dict['qcheminp'] = qcheminp
    template_dict['qchemout'] = qchemout

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

    logfp = open('/tmp/handlelinkatoms.txt', 'w')
    logfp.write('processing link atoms.\n')

    lonepairs = selection.models.LonePair.objects.filter(selection=atomselection) #All atomselections have that in them, and it's an int field so there's no reason to worry.

    for lonepair in lonepairs:
        linkqmsegid = lonepair.qmsegid
        linkqm = lonepair.qmresid
        qmatomtype = lonepair.qmatomname

        linkmmsegid = lonepair.mmsegid
        linkmm = lonepair.mmresid
        mmatomtype = lonepair.mmatomname

#    num_linkatoms = int(postdata['num_linkatoms'])
#    for i in range(num_linkatoms):
#        logfp.write('Doing link atom %d\n' % i)
#        #Watch for the _! The _ make all the difference. This is modified from the original format in lobos!
#        linkqmsegid = postdata['linkqmsegid_' + str(i)]
#        linkqm = postdata['linkqm_' + str(i)]
#        qmatomtype = postdata['qmatomtype_' + str(i)]
#
#        linkmmsegid = postdata['linkmmsegid_' + str(i)]
#        linkmm = postdata['linkmm_' + str(i)]
#        mmatomtype = postdata['mmatomtype_' + str(i)]
#
        link_tuple = (linkqmsegid,linkqm,qmatomtype,linkmmsegid,linkmm,mmatomtype)
        link_list.append(link_tuple)

    logfp.write('Done\n')
    logfp.close()
    return link_list

