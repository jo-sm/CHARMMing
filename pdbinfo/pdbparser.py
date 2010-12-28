import sys
import re
import os 
import copy

#variables:
#The variables listed below remain global
#pdb_filename stores the name of the temp-pdb being read in
#segmented_pdb_filenames stores the pdb filenames of the different segments made
pdb_filename = sys.argv[1];
#check to see if they include a file location, if not it is assumed to be the current
#working directory
try:
    fileLocation = sys.argv[2]
    if(fileLocation == "c35" or fileLocation == "c34"):
        fileLocation = "./"
        #if the user doesnt specify a file location but specifies a charmm format
        #then charmm format is actually argv[2]
        argvTwoIsCharmmFormat = 1
    else:
        fileLocation += '/'
        argvTwoIsCharmmFormat = 0
except:
    fileLocation = "./"

#will make directory if it does not exist
try:
    os.mkdir(fileLocation)
except:
    pass

#stores the user specified CHARMM format
#charmm format values can be c35, c34
#For c35, charmm_format is blank, this is done for naming conventions
#as c34 will get appended to the PDB name
try:
    #check to see if argv[2] actually has the charmm format
    if(argvTwoIsCharmmFormat):
        charmm_format = sys.argv[2]
    else:
        charmm_format = sys.argv[3];
    if(charmm_format == "c35"):
        charmm_format = ""
    else:
        #charmm_format = "-c34"
        charmm_format = ""
except:
    #blank charmm format signifies c35
    charmm_format = ""

segmented_pdb_filenames = []
pdb_line_dict = {}
listOfFileResnames = []

#gets basic information of the PDB line
def getPDBLineInfo(pdb_line):
    #pdb_line_list will convert the pdb_line into a list of characters
    #atom_number stores the atom number of the current line
    #residue_number stores the residue number of the current line
    #resname stores the residue to the current line
    #tag will be either ATOM or HETATM if the pdb file is conventional
    pdb_line_list = list(pdb_line)
    pdb_line_dict['tag'] = ''.join(pdb_line_list[0:6]).strip()
    if(pdb_line_dict['tag'] != 'HETATM' and pdb_line_dict['tag'] != 'ATOM'):
        return pdb_line_dict
    try:
        pdb_line_dict['atom_number'] = int(''.join(pdb_line_list[6:11]).strip())
    except:
        pass
    pdb_line_dict['atom_type'] = ''.join(pdb_line_list[12:16]).strip()
    try:
        pdb_line_dict['residue_number'] = int(''.join(pdb_line_list[22:27]).strip())
    except:
        pass
    pdb_line_dict['resname'] = ''.join(pdb_line_list[17:21]).strip()
    pdb_line_dict['remainder'] = ''.join(pdb_line_list[27:])
    pdb_line_dict['segid'] = pdb_line_list[21].lower().strip()
    return pdb_line_dict 

def pad_left(s, n):
   npad = n - len(s)
   if npad < 0:
      raise "Bad padding specified for %s." % s
   elif npad == 0:
      return s
   else:
      return ' ' * npad + s

def pad_right(s, n):
   npad = n - len(s)
   if npad < 0:
      raise "Bad padding specified for %s." % s
   elif npad == 0:
      return s
   else:
      return s + ' ' * npad

def make_pdb_line(pdb_line_dict):
   #if(pdb_line_dict['tag'] == 'TER' or pdb_line_dict['tag'] == 'END'):
   #    return pdb_line_dict['tag']
   new_line = pad_right(pdb_line_dict['tag'], 6)
   new_line += (pad_left(str(pdb_line_dict['atom_number']), 5) + ' ')
   # extra space blanks out the altLoc field, which we don't need
   if len(pdb_line_dict['atom_type']) == 4:
      new_line += (pdb_line_dict['atom_type'] + " ")
   else:
      new_line += (" " + pad_right(pdb_line_dict['atom_type'], 3) + " ")
   # TIP3 has 4 characters, so we may need to knock out the extra space
   new_line += (pad_right(pdb_line_dict['resname'], 4) + pdb_line_dict['segid'].upper())
   if(charmm_format == "-c34"):
       new_line += (pad_left(str(pdb_line_dict['residue_number']), 5)) 
       new_line += pdb_line_dict['remainder']
   else:
       new_line += (pad_left(str(pdb_line_dict['residue_number']), 4)) 
       new_line += " " + pdb_line_dict['remainder']
   return new_line

#removes .pdb form a file name
def stripDotPDB(old_filename):
    new_filename = old_filename
    pdb = re.compile('.pdb')
    new_filename =  pdb.sub('',new_filename)
    return new_filename

#removes .pdb form a file name
def stripDotTMP(old_filename):
    new_filename = old_filename
    pdb = re.compile('.tmp')
    new_filename =  pdb.sub('',new_filename)
    return new_filename

#renames goodhet atoms and residues to make themm CHARMM-compliant
def renameGoodHetAtomsAndResidues(pdb_line_dict):
    if(pdb_line_dict['tag'] == 'ATOM'):
        if(pdb_line_dict['resname'] == 'HOH'):
            pdb_line_dict['resname'] = 'TIP3'
        if(pdb_line_dict['resname'] == 'TIP3'):
            pdb_line_dict['atom_type'] = 'OH2'
        elif(pdb_line_dict['atom_type'] == 'ZN'):
            pdb_line_dict['resname'] = 'ZN2'
        elif(pdb_line_dict['atom_type'] == 'NA'):
            pdb_line_dict['resname'] = 'SOD'
            pdb_line_dict['atom_type'] = 'SOD'
        elif(pdb_line_dict['atom_type'] == 'CS'):
            pdb_line_dict['resname'] = 'CES'
            pdb_line_dict['atom_type'] = 'CES'
        elif(pdb_line_dict['atom_type'] == 'CL'):
            pdb_line_dict['resname'] = 'CLA'
            pdb_line_dict['atom_type'] = 'CLA'
        elif(pdb_line_dict['atom_type'] == 'CA'):
            pdb_line_dict['resname'] = 'CAL'
            pdb_line_dict['atom_type'] = 'CAL'
        elif(pdb_line_dict['atom_type'] == 'K'):
            pdb_line_dict['resname'] = 'POT'
            pdb_line_dict['atom_type'] = 'POT'
    return pdb_line_dict

#renames atom lines 
def renameAtomsAndResidues(pdb_line_dict):
    if(pdb_line_dict['resname'] == 'HIS'):
        pdb_line_dict['resname'] = 'HSD'
    if(pdb_line_dict['resname'] == 'ILE' and pdb_line_dict['atom_type'] == 'CD1'):
        pdb_line_dict['atom_type'] = 'CD'
    if(pdb_line_dict['atom_type'] == 'OXT'):
        pdb_line_dict['atom_type'] = 'OT2'

#splits the main pdb file into segmented PDB files
def splitIntoSectionFiles():
    #variables:
    #tempfile opens up the pdb file to read in
    #previous segid stores the previous line's segid
    #currrent_segid is the current segid of the line
    #list_of_segids stores the the segids tha thave passsed through
    #resname stores the residue to the current line
    #tag will be either ATOM or HETATM if the pdb file is conventional
    #fileHandleDict has segids as keys and values as opened file handles
    tempfile = open(fileLocation + pdb_filename,'r')
    previous_segid = 'null';
    list_of_segids = [];
    fileHandleDict = {}
    for pdb_line in tempfile:
        getPDBLineInfo(pdb_line)
        #removes extraneous lines, like Header lines by just ignoring them
        if(pdb_line_dict['tag'] != 'HETATM' and pdb_line_dict['tag'] != 'ATOM'):
            continue
        #Current resname extracted just for code's cleanliness sake
        resname = pdb_line_dict['resname']
        #if the current segid doesn't equal the segid from the previous line (signifies a potential new segment)
        if(pdb_line_dict['segid'] != previous_segid):
            #Check to make sure the current segid is not blank. If it is blank it gets associated with the
            #segid before it
            if(pdb_line_dict['segid'] == '' or pdb_line_dict['segid'] == ''):
                pdb_line_dict['segid'] = previous_segid;
            previous_segid = pdb_line_dict['segid']
        #cross check the residue to see if it is a "good-het"
        #if it is a CHARMM-recognized  HETATM residue then it is a "good-het"
        if(resname == 'TIP3' or resname == 'ZN' or resname == 'FE' or resname == 'NA' \
        or resname == 'CA' or resname == 'MG' or resname == 'CS' \
        or resname == 'K' or resname == 'CL' or resname == 'HOH'):
            pdb_line_dict['tag'] = 'ATOM'
            #some PDB residues/atom types may not be CHARMM compliant, so the below method makes them
            #compliant
            renameGoodHetAtomsAndResidues(pdb_line_dict) 
            #if this is the first time a segid is found, store the segid as a key in a dictionary
            #and the filehandle will be the value
            fileHandleKey = pdb_line_dict['segid'] + '-goodhet'
            new_pdb_segment_filename = 'new_' + stripDotPDB(pdb_filename) + \
                                       '-' + pdb_line_dict['segid'] + '-goodhet' + charmm_format + '.pdb.tmp'
            if(fileHandleKey not in fileHandleDict.keys()):
                fileHandleDict[fileHandleKey] = open(fileLocation + new_pdb_segment_filename,'a')
            #write to PDB file
            fileHandleDict[fileHandleKey].write(make_pdb_line(pdb_line_dict))
        #if its not a goodhet, check to see if it is just a hetatm. If it is just a 
        #hetatm then it is a non-CHARMM recognized hetatm
        elif(pdb_line_dict['tag'] == 'HETATM'):
            #if this is the first time a segid is found, store the segid as a key in a dictionary
            #and the filehandle will be the value
            fileHandleKey = pdb_line_dict['segid'] + '-het'
            new_pdb_segment_filename = 'new_' + stripDotPDB(pdb_filename) + \
                                       '-' + pdb_line_dict['segid'] + '-het' + charmm_format + '.pdb.tmp'
            if(fileHandleKey not in fileHandleDict.keys()):
                fileHandleDict[fileHandleKey] = open(fileLocation + new_pdb_segment_filename,'a')
            #write to PDB file
            fileHandleDict[fileHandleKey].write(make_pdb_line(pdb_line_dict))
        #if it is not a hetatm or het then set this as the PDB filename
        elif(pdb_line_dict['tag'] == 'ATOM'):
            #some PDB residues/atom types may not be CHARMM compliant, so the below method makes them
            #compliant
            renameAtomsAndResidues(pdb_line_dict) 
            #if this is the first time a segid is found, store the segid as a key in a dictionary
            #and the filehandle will be the value
            fileHandleKey = pdb_line_dict['segid']
            new_pdb_segment_filename = 'new_' + stripDotPDB(pdb_filename) + \
                                       '-' + pdb_line_dict['segid'] + charmm_format + '.pdb.tmp'
            if(fileHandleKey not in fileHandleDict.keys()):
                fileHandleDict[fileHandleKey] = open(fileLocation + new_pdb_segment_filename,'a')
            #write to PDB file
            fileHandleDict[fileHandleKey].write(make_pdb_line(pdb_line_dict))
        #record filenames of segments created if it does not already exist
        if(new_pdb_segment_filename not in segmented_pdb_filenames):
            segmented_pdb_filenames.append(new_pdb_segment_filename);
    tempfile.close()
    #fileHandleDict.keys() is just segids
    print fileHandleDict.keys()
    for key in fileHandleDict.keys():
        #fileHandleDict[key].write("END")
        fileHandleDict[key].close()

#renumbers the atom columns and residue columns
def renumberPDBFiles():
    #variables
    #prev_orig_line_dict stores the previous line without any modifications
    #prev_mod_line_dict stores the previous line with modifications such as fixing of line numbering
    #curr_orig_line_dict stores the current line without any modifications
    #curr_mod_line_dict stores the current line with modifications such as fixing of line numbering
    #new_pdb_file_data stores all the lines of the PDB that will get written to the new
    #             PDB file at the end. This is done so the PDB file does not have
    #             to be opened during every iteration
    new_pdb_line = ''
    for pdb_filename in segmented_pdb_filenames:
        #print "processing " + pdb_filename
        prev_orig_line_dict = {}
        prev_mod_line_dict = {}
        curr_orig_line_dict = {}
        curr_mod_line_dict = {}
        #tempfile opens up the pdb file to read in
        tempfile = open(fileLocation + pdb_filename,'r')
        new_pdb_file_handle = open(fileLocation + stripDotTMP(pdb_filename),'a')

        for pdb_line in tempfile:
            curr_orig_line_dict = getPDBLineInfo(pdb_line)
            curr_mod_line_dict = copy.deepcopy(curr_orig_line_dict)
            #if there are no keys in the previous original line, aka if this is the
            #first line, then we have to give it some default values
            if(len(prev_orig_line_dict.keys()) == 0):
                #We have to ensure the numbering of certain fields like atom number and residue
                #number is correct (in this case starting at 1)
                curr_mod_line_dict['residue_number'] = 1
                curr_mod_line_dict['atom_number'] = 1

                #print "A Curr res number = %d" % curr_mod_line_dict['residue_number']
            #Basically, if this is not the first line in the pdb, do the following
            else:
                #print "B Prev res number = %d" % prev_mod_line_dict['residue_number']

                #If the atom number is not increasing by one in the current line, then make it the previous line's atom 
                #number plus 1
                if(curr_orig_line_dict['atom_number'] - prev_mod_line_dict['atom_number'] != 1):
                    curr_mod_line_dict['atom_number'] = prev_mod_line_dict['atom_number'] + 1

                #If the current residue and the previous residues are different and their resids are nonsequential (!= 1)
                #make them sequential
                orig_resnum_delta = curr_orig_line_dict['residue_number'] - prev_orig_line_dict['residue_number']
                new_resnum_delta = curr_orig_line_dict['residue_number'] - prev_mod_line_dict['residue_number']
                new_resname = curr_orig_line_dict['resname'] != prev_orig_line_dict['resname']

                # case 1: we have a new residue, but the number is non-sequential
                # solution: just add 1 to the previous residue number
                if orig_resnum_delta != 1 and new_resname:
                    curr_mod_line_dict['residue_number'] = prev_mod_line_dict['residue_number'] + 1

                # case 2: we have a new residue, and the original numbers are sequential 
                # solution: as above, add one to the previous midified residue number
                elif orig_resnum_delta == 1:
                    curr_mod_line_dict['residue_number'] = prev_mod_line_dict['residue_number'] + 1

                # case 3: the residue number is the same as the previous line
                # solution: keep it consistent
                elif orig_resnum_delta == 0 and not new_resname:
                    curr_mod_line_dict['residue_number'] = prev_mod_line_dict['residue_number']

                # otherwise, the residue number is good as-is


                #if(curr_orig_line_dict['residue_number'] - prev_mod_line_dict['residue_number'] != 1 and \
                #   curr_orig_line_dict['resname'] != prev_orig_line_dict['resname']):
                #    curr_mod_line_dict['residue_number'] = prev_mod_line_dict['residue_number'] + 1
                #    print "B.2 change currently set to %d" % curr_mod_line_dict['residue_number']
                #if the residue numbers in the original PDBs are the same (same residue), but the numbering is off
                #this will fix the numbering. For example if a PDB has one residue but its residue number starts at
                #150, this will make the entire PDB with one residue have a residue number of 1 ex: PDB.org ID 1O1O
                #elif(curr_orig_line_dict['residue_number'] == prev_orig_line_dict['residue_number'] and \
                #     curr_orig_line_dict['residue_number'] - prev_mod_line_dict['residue_number'] != 1):
                #    curr_mod_line_dict['residue_number'] = prev_mod_line_dict['residue_number']
                #    print "B.3 change now set to %d" % prev_mod_line_dict['residue_number']
                #elif(curr_orig_line_dict['residue_number'] - prev_mod_line_dict['residue_number'] == 1):
                #    curr_mod_line_dict['residue_number'] 

                #print "B Curr res number = %d" % curr_mod_line_dict['residue_number']

            #previously modified dictionary becomes the currently modified dictionary
            prev_mod_line_dict = copy.deepcopy(curr_mod_line_dict)
            #previous original (unmodified) dictionary becames the currently original (unmodified) dictionary
            prev_orig_line_dict = copy.deepcopy(curr_orig_line_dict)
            new_pdb_file_handle.write(make_pdb_line(curr_mod_line_dict))

        new_pdb_file_handle.write("TER")
        new_pdb_file_handle.close()
        tempfile.close()

    return listOfFileResnames





splitIntoSectionFiles();
renumberPDBFiles();
