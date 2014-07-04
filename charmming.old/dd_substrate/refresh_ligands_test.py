#!/usr/bin/python
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

import sys, os, string, getopt, commands, glob, common
from dd_substrate import common
#dd stuff
from dd_substrate.models import ligands 
from dd_infrastructure.models import files,files_objects,file_types, sources
###


# main routine
if __name__ == '__main__':
    
   


def getMoleculeTitle(file):
    while 1:
        line = file.readline()
        if not line:
            break
        if line.strip() == "@<TRIPOS>MOLECULE":
            titleline.readline()
            title=titleline.strip()
            break
    return title
