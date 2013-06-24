
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
# Create your views here.
from django import forms
from django.template.loader import get_template
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from httplib import HTTPConnection
from account.models import *
from structure.models import Task 
from django.contrib.auth.models import User
from django.template import *
from account.views import checkPermissions
from lessons.models import LessonProblem
from mutation.models import mutateTask
import output, lesson1, lesson2, lesson3, lesson4, lessonaux
import structure.models, structure.views, input
from lesson_config import *
import os, copy, json, mimetypes, string, re
import pychm.io, charmming_config
from pychm.io.pdb import PDBFile
import cPickle
import openbabel

def clear_struct(request): #To write less repetitive code...
    return (None, charmming_config.user_home + "/" + request.user.username + "/")

def build_ligand(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')
#    logfp = open("/tmp/ligand_design.txt","w")
    segletter = "a"
    seggoodletter = "a"
    if request.POST and request.POST.has_key('attach_check'):
#        logfp.write("Found attach_check.\n")
        try:
            struct = structure.models.Structure.objects.get(owner=request.user,selected='y')
            logfp.write("Structure found.\n")
            filepath = struct.location
            sl = []
            sl.extend(structure.models.Segment.objects.filter(structure=struct))
            maxletter = ""
            maxgoodletter = "" #For GOOD segments
            for seg in sl: #This way we limit ourselves to the "bad" hetatms
                if seg.type == "bad" and seg.name[0] > maxletter:
                    maxletter = seg.name[0]
                if seg.type == "good" and seg.name[0] > maxgoodletter:
                    maxgoodletter = seg.name[0]
    #        logfp.write("name: " + maxletter)
            if maxletter == "z" or maxgoodletter == "z":
                return HttpResponse("Too many segments for this structure.") #This is only a minor setback! We shall figure out how to handle these later.
            else:
                if maxletter: #if there's an issue, default to 'a'
                    segletter = chr(ord(maxletter) + 1) #Can't handle anything but standard letters, but that's fine, neither can CHARMM.
                if maxgoodletter:
                    seggoodletter = chr(ord(maxgoodletter) + 1)
        except Exception as e: #If you get here, the user tried to attach to structure when there are no structures uploaded.
            logfp.write("No structures found.\n")
            struct, filepath = clear_struct(request)
        finally:
            logfp.close()
    else:
        struct, filepath = clear_struct(request)
    tdict = {}
#    logfp.write(str(struct) + "\n")
    tdict['super_user'] =request.user.is_superuser
    if request.POST: #Do stuff with the JSmol data if the POSTdata has some
        os.chdir(filepath)
        postdata = request.POST
        molname = postdata['LigName']
        moldata = postdata['molinfo']
        if 'attach_check' in postdata.keys():
            tdict['attach_check'] = True
        else:
            tdict['attach_check'] = False
        lines = moldata[moldata.index("HETATM"):moldata.index("]")]
#        logfp.write(str(lines) + "\n\t")
        lines = lines.split("\\n")
#        logfp.write(str(lines) + "\n\t")
        tmpfile = open("moldata.pdb","w") #We use a single tempfile to avoid encoding issues.
        end_data = ""
        for line in lines:
            if line.startswith("COMPND"):
                line = "COMPND    LIGAND: " + molname + "\n"
            end_data += line + "\n"
        tmpfile.write(end_data)
        tmpfile.close()
        obConversion = openbabel.OBConversion()
        obConversion.SetInFormat("pdb")
#        obConversion.SetInAndOutFormats("mol", "pdb") #So we can write to file for PDBFile if this is successful
        mol = openbabel.OBMol()
        obConversion.ReadFile(mol, "moldata.pdb")
        charge_model = openbabel.OBChargeModel.FindType("mmff94")
        charge_model.ComputeCharges(mol)
        partial_charges = charge_model.GetPartialCharges()

        molec_charges = 0.0
        core_charges = 0.0
        for atom in openbabel.OBMolAtomIter(mol):
            atom_idx = atom.GetIdx()
            molec_charges += atom.GetPartialCharge()
            core_charges += atom.GetAtomicNum()
#            logfp.write(str(atom_idx) + "\t" + str(atom.GetValence()) + "\n")
#        logfp.write("Molecular charge:\t" + str(molec_charges) + "\n")
#        logfp.write("Core charges:\t" + str(core_charges) + "\n")
#        logfp.write(str(charmming_config.user_home))
        tdict['filepath'] = '/charmming/pdbuploads/' + request.user.username + "/"
        if struct:
            tdict['filepath'] = tdict['filepath'] + struct.name + "/"
        tdict['moldata'] = moldata #pass them back and forth until we have valid data
        tdict['molname'] = molname
        if (core_charges - molec_charges) % 2 != 0 and postdata['force_charge'] == "false":
            tdict['charge_alert'] = True
            return render_to_response('html/ligand_design.html', tdict)
        #If force_charge is true we just pass the molecule on regardless.
        else:
            tdict['charge_alert'] = False
            resname = ""
            if len(molname) >= 3:
                resname = molname[0:3]
            else:
                resname = molname
            ligname = resname #This way you don't look for YZQ B_IDEAL
            resname = resname.upper() + (" " * (4-len(resname))) + segletter.upper()
            force_custom = postdata['force_custom'] == "true"
            #If set to false, we check whether the molecule exists in PDB.org,
            #and if it does, warn the user.
            #If set to true, the user has approved of the molecule, so we skip the if check.
            if not(force_custom): #i.e. if the user has not said he is sure this structure is correct
                conn = HTTPConnection("www.pdb.org")
                reqstring = "/pdb/files/ligand/%s_ideal.sdf" % ligname
                conn.request("GET", reqstring)
                resp = conn.getresponse()
                if resp.status == 200:
#                    logfp.write("I found a PDB for " + ligname)
#                    logfp.close()
                    tdict['sdf_link'] = "www.pdb.org" + reqstring
                    return render_to_response('html/ligand_design.html',tdict)
            #If either force_custom is True or it wasn't needed, the code just continues.
#            outstring = obConversion.WriteString(mol) #Consider changing these to ligand names
            #Now we need to add a chainid so that PDBFile doesnt' freak out.
            #Get ligand name, first 4 characters, add spaces.
            end_data = end_data.replace("UNK  ", resname)
            outfile = open("moldata.pdb","w")
            outfile.write(end_data)
            outfile.close()
            ligpdb = PDBFile("moldata.pdb")
            ligmol = ligpdb['model00']
            for seg in ligmol.iter_seg():
                if seg.segType == "pro": #If the builder caught it as GOOD/BAD, leave it, else change
                    seg.segType = "good" #If it became PRO, then we have params for it, so make it good
#            logfp.write(str(struct) + "\n")
            #Now we do a check for a structure being present
            if struct:
                #Pickle time.
                woof = open(struct.pickle)
                fullpdb = cPickle.load(woof)
                woof.close()
                thisMol = fullpdb['model00'] #Model 0 again...
#                logfp.write(str(thisMol) + "\n")
                #first a sanity check so that the user doesn't build more than one residue with the same name
                for seg in thisMol.iter_seg():
                    for res in seg.iter_res():
#                        logfp.write("residue: " + str(res.resName) + "\n")
#                        logfp.write("ligand residue: " + str(ligname.lower()) + "\n")
                        if res.resName == ligname.lower():
                            tdict['same_residue'] = True
                            return render_to_response('html/ligand_design.html',tdict)
#                for mol in thisMol:
#                    nyan = str(mol)
#                    logfp.write(nyan + "\n")
#                logfp.close()
                thisMol.extend(ligmol) #Snap them together and we have the real deal.
                #Now we need to save segments and dump the pickle.
                structure.views.getSegs(thisMol,struct,auto_append=False) #This saves all the Segments and stuff
                woof = open(struct.pickle,"w")
                cPickle.dump(fullpdb,woof) #we took the object, extended it, now we save
                woof.close()
                os.unlink("moldata.pdb")
#                os.unlink("moldata.mol")
                struct.save()
            else: #No structure, so we build one
                struct = structure.models.Structure()
                struct.owner = request.user
                struct.title = "Custom Ligand " + molname

                location = filepath
                dname = molname[0:3].lower()
               # tmpdname = dname
                #version = -1
                #THis stuff breaks. You shouldn't be building different ligands with same name anyway...
                #while os.path.exists(location + "/" + dname):
                #    version += 1
                #    dname = tmpdname + "-" + str(version)

                struct.name = dname
                struct.location = location + dname
                try:
                    os.mkdir(struct.location)
                except:
                    tdict['struct_error'] = True
                    return render_to_response('html/ligand_design.html',tdict)
                os.chmod(struct.location, 0775)

                fullpath = struct.location + "/" + dname + ".pdb"
                structure.views.getSegs(ligmol,struct,auto_append=False)
                pfname = location + dname + "/" + "pdbpickle.dat"
                picklefile = open(pfname, "w")
                cPickle.dump(ligpdb,picklefile)
                picklefile.close()
                struct.pickle = pfname

                try:
                    oldfile = structure.models.Structure.objects.get(owner=request.user,selected='y')
                    oldfile.selected = "n"
                    oldfile.save()
                except:
                    pass

                struct.selected = "y"
#                os.unlink("moldata.pdb")
#                os.unlink("moldata.mol")
                struct.save()
            return HttpResponseRedirect('/charmming/buildstruct/')
    else: #Render the page...
        return render_to_response('html/ligand_design.html', tdict)