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

# tideInp takes an input script from a template a neatens it up (removes
# leading spaces, excessive blank lines, and any line with just "blankme").

from django.shortcuts import render_to_response

def returnSubmission(jobtitle, error=None):
    td = {}
    if error:
        td['joberror'] = True
        td['errortext'] = error
    else:
        td['joberror'] = False

    return render_to_response('html/jobsubmit.html', td)

def tidyInp(intxt):
    blcnt = 0
    larr = intxt.split('\n')
    outp = ''

    for line in larr:
        line = line.strip()
        if not line:
            blcnt += 1
        else:
            blcnt = 0
        if blcnt > 1 or line == 'blankme':
            continue
        if line.startswith("CHARMM-scripts"):
            line = line[14:]  
        outp += "%s\n" % line

    return outp

def tidyTxt(intxt):
    blcnt = 0
    larr = intxt.split('\n')
    outp = ''

    for line in larr:
        if line.strip() == 'blankme':
            continue
        if line:
            outp += "%s\n" % line
        else:
            outp += "\n"  
    return outp
