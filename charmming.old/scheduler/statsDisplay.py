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
def statsDisplay(statline,jobid):
    try:
       statarr = statline.split()
    except:
       return "Unknown"
    status = statarr[4] 
    if status == 'submitted':
        return "<span style='color:FFCC33'>Submitted</span>"
    elif status == 'queued':
        #return "<span style='color:FFCC33'>Queued</span> <a href=\"/charmming/killjob/%d\">X</a>" % jobid
        return "<span style='color:FFCC33'>Queued</span> <a href=\"javascript:killJob(%d);\">X</a>" % jobid
    elif status == 'running':
        #return "<span style='color:3333FF'>Running</span> <a href=\"/charmming/killjob/%d\">X</a>" % jobid
        return "<span style='color:3333FF'>Running</span> <a href=\"javascript:killJob(%d);\">X</a>" % jobid
    elif status == 'complete':
        return "<span style='color:33CC00'>Done</span>"
    elif status == 'failed':
        return "<font color=red>Failed</span>"

    print "Achtung! Couldn't figure out this status (%s)!!!!" % status
    return "<font color=\"red\">Unknown</span>"

