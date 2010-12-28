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
import os

#returns the number of normal mode movies
def getNormalModeMovieNum(file):
    trjcount = 1
    continueit = True
    while(continueit):
        try:
	    os.stat(file.location + 'new_' + file.stripDotPDB(file.filename) + '-nma-mainmovie-' + str(trjcount) + '.pdb')
	    trjcount += 1
	except:
	    continueit = False
    return trjcount-1

def parseNormalModes(file):
    #read in the normal modes output file

    infp = open(file.location + 'charmm-' + file.stripDotPDB(file.filename) + '-nmodes.out','r')
    outfp = open(file.location + '10modes-' + file.stripDotPDB(file.filename) + '.txt','w')
    vibrationmode = re.compile('VIBRATION MODE')
    linecount = 0
    vibrationcount = 0
    for line in infp:
        if linecount == 3:
	    linecount = 0
        #second and third lines after NORMAL MODES has the first 10 modes
        elif linecount == 1 or linecount == 2:
	    outfp.write(line.lstrip())
	    linecount += 1
	elif vibrationmode.search(line):
	    if vibrationcount == 10:
	        break
	    outfp.write(line.lstrip())
	    vibrationcount += 1
	    linecount +=1
	    continue
        else:
            continue
    infp.close()
    outfp.close()

