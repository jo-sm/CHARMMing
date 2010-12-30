import os
import re
from django import forms
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response

def checkNterPatch(file,segment):

    # custom sequences must be handled differently
    if segment == "sequ-pro":
        sequ_filename = "new_" + file.stripDotPDB(file.filename) + "-sequ-pro.pdb"
        sequ_handle = open(file.location + sequ_filename,'r')
        sequ_line = sequ_handle.read()
        sequ_line.strip()
        resid = sequ_line.split(' ')[0]
    else:
        # we can't count on the -final files to have been created yet, but the parsed PDB
        # files will do OK
        try:
            segpdb = open(file.location + "new_" + file.stripDotPDB(file.filename) + "-" + segment + ".pdb", "r")
        except:
            # NTER is gonna be the default if something goes wrong -- lame, huh?
            return "nter"
        if not segpdb:
            return "nter"
        # read until the first line beginning with atom (this should be the very first line of the
        # the file, but this is defensive programming, kiddos.
        for line in segpdb.readlines():
            if line.startswith('ATOM'):
                break
        segpdb.close()
        try:
            resid = line[18:21]
        except:
            # should probably throw an error instead
            return "nter"

    # figured out the first resid in this segment, now let's decide what patch it needs
    resid = resid.upper()
    if resid == "GLY":
        return "glyp"
    elif resid == "PRO":
        return "prop"
    else:
        return "nter"

def parseEnergy(file,output_filename):
    outfp = open(file.location + output_filename,'r')
    ener = re.compile('ENER ENR')
    extern = re.compile('EXTERN')
    charmm = re.compile('CHARMM>')
    initiate = 0
    #the second extern marks the end of the energy output
    extern_occurance = 0
    energy_lines = ''
    for line in outfp:
        if(extern_occurance == 2 or (initiate >0 and line.strip().startswith('CHARMM>'))):
            if(line.strip().startswith('CHARMM>')):
                break
            energy_lines += line.lstrip()
            break
        if(ener.search(line) or initiate > 0):
            if(initiate == 0):
                initiate = 1
            energy_lines += line.lstrip()
        if(extern.search(line)):
            extern_occurance += 1
    outfp.close()
    writefp = open(file.location + 'energy-' + file.stripDotPDB(file.filename) + '.txt','w')
    writefp.write(energy_lines)
    writefp.close()
    return render_to_response('html/displayenergy.html',{'linelist':energy_lines})

