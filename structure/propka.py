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

import structure.models
import charmming_config, input
import commands, datetime, sys, re, os, glob, shutil
import pychm.io

#This file exists because of PROPKA for the Build Structure page ONLY!
#All functions in here are to do with using PROPKA for that page
#for the PROPKA page, see charmming/propka/views.py

def calculate_propka_residues(struct):
    #Outline goes as follows:
    #run PROPKA on the PDB from the protein (if any...if there isn't one, i.e. just PSF/CRD, we can't do anything since
    #babel used to have CHARMM format in 2002, but 12 years later they still haven't re-added support for it for openbabel 2.0)
    #Anyway!
    #we run PROPKA, and use data garnished/stolen from organic chem textbooks on the pKa of these things
    # (we can probably cite him if need be)
    # and use that to predict the proto state. If < the pKa of a certain proto state, it's protonated, else, it's deprotonated
    # if pKa = the turning point, tell the user that there's a problem
    #basically we're going to append one last thing to the tuple containing the resname to mark what species it should be
    #now we check for whether we have a PDB file
    glup_list = [] #allocate these just in case python decides to keep them inside the if statement
    hsp_list = []
    lsn_list = []
    asp_list = []
    glup_pka = 4.07 #side_chain pKa for glutamates
    hsp_pka = 6.00 #side_chain pKa for histidines
    lsn_pka = 10.80
    asp_pka = 4.10 #should be close enough to glup_pka
    os.chdir(struct.location) #have to change to here because PROPKA only outputs to the dir it's run in
    propka_output = open(struct.original_name + "-propka.pka")
    #we make a list holding data on each of the four residues we care about
    #and then we make a regex for parsing them
    space_regex = re.compile(r' +')
    found_summary = False #We use the "SUMMARY OF THIS PREDICTION" lines because they're formatted better
    user_decision = False
    for line in propka_output:
        if line.startswith("       Group"): #We use this so that we're on the right line...using SUMMARY puts us one line above the useless table headers
            found_summary = True
            continue
        if found_summary:
            if line.startswith("------------"): #close enough to the terminator line
                break
            else:
                line_pieces = space_regex.split(line)
                #Now we dump the resno pKa value and chain and make lists
                info_list = []
                info_list.append(int(line_pieces[2]))
                info_list.append(line_pieces[3].lower())
                info_list.append(line_pieces[4])
                species_pka = float(line_pieces[4]) #make our lives easier here
                #Now we do checks, and then add statuses depending on the case
                #if it should stay at the default value, we write 0
                # if it should be protonated, we write 1
                # if the user needs to decide (equal pKa), we write 2
                #these numbers will be used in the web template to modify the HTML depending.
                if line_pieces[1] == "ASP":
                    if species_pka < asp_pka:
                        info_list.append(1) #neutral
                    elif species_pka > asp_pka:
                        info_list.append(0) #-1 charge
                    else:
                        info_list.append(2) #User needs to decide
                        user_decision = True
                    asp_list.append(info_list)
                    continue #non-op, but we may as well be safe
                elif line_pieces[1] == "GLU":
                    if species_pka < glup_pka:
                        info_list.append(1) #neutral
                    elif species_pka > glup_pka:
                        info_list.append(0) #-1 charge
                    else:
                        info_list.append(2) #User needs to decide
                        user_decision = True
                    glup_list.append(info_list)
                    continue
                elif line_pieces[1] == "HIS":
                    if species_pka < hsp_pka:
                        info_list.append(1) #protonated
                    else:
                        info_list.append(2) #we don't know how to decide between HSD and HSE
                        user_decision = True
                    hsp_list.append(info_list)
                    continue
                elif line_pieces[1] == "LYS":
                    if species_pka < lsn_pka:
                        info_list.append(1) #protonated
                    elif species_pka > lsn_pka:
                        info_list.append(0) #neutral
                    else:
                        info_list.append(2) #User needs to decide
                        user_decision = True
                    lsn_list.append(info_list)
                    continue
    propka_output.close()
    return [asp_list,glup_list,hsp_list,lsn_list],user_decision
    #These are all sequentially ordered, so there should be no problems later
