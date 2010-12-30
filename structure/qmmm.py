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
def checkQMregion(file, fname, selection):
   qmre = re.compile("QM REGION: (.*?)$")
   fpath = file.location + fname

   try:
      fp = open(fname, 'r')
   except:
      return 2
   if not fp:
      return 2

   qmsele = ""
   for line in fp:
      if line.startswith('*'):
         line = line.strip()
         if "QM REGION" in line:
            m = qmre.search(line)
            if m:
               qmsele = qmsele + " " + m.group(1)
      else:
         break 
   fp.close()

   # remove multiple/unneeded whitespaces
   " ".join(qmsele.split(" "))
   qmsele = qmsele.strip()

   if not qmsele:
      return 2
   if qmsele != selection.upper():
      return 0
   return 1

def wrtQchemInp(fname, exchange, correlation, basis, jobtype, charge, multi, fullhess = False):
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

def makeQChem_tpl(template_dict, file, exch, corr, bs, qmsel, jobtype, charge, multi, last_struct, linkatoms = None):
   template_dict['bad_exchange'] = 0 
   template_dict['bad_correlation'] = 0
   template_dict['bad_basis'] = 0
   template_dict['bad_jobtype'] = 0  
   if exch != "HF" and exch != "B" and exch != "B3":
      template_dict['bad_exchange'] = 1
      return template_dict
   if corr != "None" and corr != "Lyp":
      template_dict['bad_correlation'] = 1  
      return template_dict
   if bs != "STO-3G" and bs != "3-21G*" and bs != "6-31G*":
      template_dict['bad_basis'] = 1
      return template_dict
   if jobtype != "Force" and jobtype != "SP" and jobtype != "freq":
      template_dict['bad_jobtype'] = 1
      return template_dict

   qchemin  = "qchem-" + file.stripDotPDB(file.filename) + ".in"
   qcheminp = "qchem-" + file.stripDotPDB(file.filename) + ".inp"
   qchemout = "qchem-" + file.stripDotPDB(file.filename) + ".out"
   qmcheck = checkQMregion(file, last_struct, "SELE " + qmsel + " END")

   template_dict['qmcheck'] = qmcheck 
   if qmcheck == 0:
      return template_dict

   fullhess = False
   template_dict['jobtype'] = jobtype 
   if jobtype == 'freq':
       tqm = qmsel.strip()
       tqm = tqm.lower()
       template_dict['tqm'] = tqm
       if tqm != 'all':
           fullhess = True

   wrtQchemInp(file.location + qchemin, exch, corr, bs, jobtype, charge, multi, fullhess)
   template_dict['qmsel'] = qmsel

   # user-specified link atoms are passed in as a list of 6-ples (qm-segid, qm-resid, qm-type,
   # mm-segid, mm-resid, mm-type)
   # If qmcheck is 1, that means that the correct link atoms are already in the structure and
   # we should not re-add them.
   template_dict['linkatom_list'] = []
   template_dict['linkatoms'] = ''
   if qmcheck != 1:
      template_dict['linkatoms'] = linkatoms
      if linkatoms:
         for las in linkatoms:
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

# this function sets up Q-Chem for QM-MMin an input script...
# args:
# charmm_inp
# file: PDBFile struct
# exch: exchange (HF, B, B3)
# corr: correlation ("None" or Lyp)
# bs: basis set (STO-3G, 3-21G*, 6-31G*)
# qmsel: QM region, written as an atom selection without the sele ... end
# jobtype: "force" or "sp"
# charge: charge on QM region
# multi: multiplicity
# linkatoms (optional parameter): user specified link atoms -- if not passed,
#                                 we attempt to figure these out for ourselves.
def makeQChem(charmm_inp, file, exch, corr, bs, qmsel, jobtype, charge, multi, last_struct, linkatoms = None):

   if exch != "HF" and exch != "B" and exch != "B3":
      charmm_inp += "bomlev 3\n"
      charmm_inp += "achtung, bad QM exchange!\n"
      return charmm_inp
   if corr != "None" and corr != "Lyp":
      charmm_inp += "bomlev 3\n"
      charmm_inp += "achtung, bad QM correlation!\n"
      return charmm_inp
   if bs != "STO-3G" and bs != "3-21G*" and bs != "6-31G*":
      charmm_inp += "bomlev 3\n"
      charmm_inp += "achtung, bad QM basis set!\n"
      return charmm_inp
   if jobtype != "Force" and jobtype != "SP" and jobtype != "freq":
      charmm_inp += "bomlev 3\n"
      charmm_inp += "achtung, bad QM job type!\n"
      return charmm_inp
   qchemin  = "qchem-" + file.stripDotPDB(file.filename) + ".in"
   qcheminp = "qchem-" + file.stripDotPDB(file.filename) + ".inp"
   qchemout = "qchem-" + file.stripDotPDB(file.filename) + ".out"

   qmcheck = checkQMregion(file, last_struct, "SELE " + qmsel + " END")
   if qmcheck == 0:
      charmm_inp += "bomlev 3\n"
      charmm_inp += "achtung, you changed the QM region from a previous calculation. This isn't allowed!\n"
      return charmm_inp

   fullhess = False
   if jobtype == 'freq':
       tqm = qmsel.strip()
       tqm = tqm.lower()
       if tqm != 'all':
           fullhess = True
           charmm_inp += "! check to make sure the QM region is not all atoms (if it's not sele all end)\n"
           charmm_inp += "! unless you select all the Q-Chem script assumes that there are some MM atoms.\n"
           charmm_inp += "define junk sele %s end\n" % qmsel
           charmm_inp += "if ?nsel .ne. ?natom then goto qmok\n"
           charmm_inp += "achtung, QM region is all atoms, please use ALL as your selection!\n"
           charmm_inp += "label qmok\n\n"

   wrtQchemInp(file.location + qchemin, exch, corr, bs, jobtype, charge, multi, fullhess)

   charmm_inp += "define qmregion sele " + qmsel + " end\n"
   charmm_inp += "define mmregion sele .not. qmregion end\n"

   # user-specified link atoms are passed in as a list of 6-ples (qm-segid, qm-resid, qm-type,
   # mm-segid, mm-resid, mm-type)
   # If qmcheck is 1, that means that the correct link atoms are already in the structure and
   # we should not re-add them.
   if qmcheck != 1:
      if linkatoms:
         for las in linkatoms:
            try:
               qmseg, qmres, qmtyp, mmseg, mmres, mmtyp = las[0:6]
            except:
               return -1
            # TODO, validate that these are actually OK...
            charmm_inp += "coor copy comp\n"
            charmm_inp += "coor shake\n"
            charmm_inp += "addl qqh1 %s %s %s %s %s %s\n" % (qmseg,qmres,qmtyp,mmseg,mmres,mmtyp)
      else:
         # we need to figure out link atoms for ourselves
         charmm_inp += """
! Define MM Host for Link Atom
define mmh select .bonded. qmregion .and. mmregion show end

! loop over bonded atoms
set nmmh ?nsel
if @nmmh .eq. 0 then goto linkdone

set i 0
label mmhloop
incr i by 1
definte junk select mmh .subset. @i show end
set mmhost@i ?SELATOM
if @i .lt. @nmmh goto mmhloop

define qmh select .bonded. mmh .and. qmregion show end

set nqmh ?nsel
set i 0
label qmhloop
incr i by 1
define junk select qmh .subset. @i show end
set qmhost@i ?SELATOM
if @i .lt. @nqmh goto qmhloop

set i 1
set j 1
label qmhlink
   label mmhlink

      addl qqh@i bynum @qmhost@@j bynum @mmhost@@i
      lonepair colinear scaled -.7261 bynum qqh@i bynum @qmhost@@j bynum @mmhost@@i

      incr i by 1
      decr nmmh by 1
   if @nmmh .gt. 0 then goto mmhlink
   incr j by 1
   decr nqmh by 1
if @nqmh .gt. 0 then goto qmhlink

coor copy comp
coor shake

label linkdone
"""
   else:
      # since link atoms are already present, we need to do the coor shake bit
      charmm_inp += "coor copy comp\n"
      charmm_inp += "coor shake\n"

   charmm_inp += """

!---------- Needed to define Q-Chem env. vars. ----------
 envi qchemcnt  \"""" + qchemin + """\"
 envi qcheminp  \"""" + qcheminp + """\"
 envi qchemexe  "qchem"
 envi qchemout  \"""" + qchemout + """\"
!--------------------------------------------------------

qchem remove sele qmregion end

"""

   return charmm_inp


# pre:requires postdata that has some patch information
# Writes patch to an external file and processes it
def handleLinkAtoms(file,postdata):
    seg_list = file.segids.split()
    #used in disulfide bond patching
    linkqmsegid = ''
    linkmmsegid = ''
    linkqm = ''
    qmatomtype = ''
    qmatomtype = ''
    linkmm = ''
    link_list = []
    link_tuple = ()
    for seg in seg_list:
        num_linkatoms = int(postdata['num_linkatoms'])
        for i in range(num_linkatoms):
            linkqmsegid = postdata['linkqmsegid' + str(i)]
            linkqm = postdata['linkqm' + str(i)]
            qmatomtype = postdata['qmatomtype' + str(i)]
            qmatomtype = postdata['qmatomtype' + str(i)]
            linkmm = postdata['linkmm' + str(i)]
            linkmmsegid = postdata['linkmmsegid' + str(i)]
            linkmmsegid = postdata['linkmmsegid' + str(i)]
            link_tuple = (linkqmsegid,linkqm,qmatomtype,linkmmsegid,linkmm,qmatomtype)
            link_list.append(link_tuple)
    return link_list

