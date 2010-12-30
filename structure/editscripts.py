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
from pdbinfo.models import PDBFile
from pdbinfo.aux import parseEnergy
from django.http import HttpResponseRedirect, HttpResponse
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
import minimization
import solvation
import normalmodes
import dynamics
import account

#This submits job to PBS/Condor
def runJob(request):
    if not account.views.isUserTrustworthy(request.user):
        return False
    file =  PDBFile.objects.filter(owner=request.user,selected='y')[0]
    scriptlist = request.POST['scriptlist'].split(',')
    jobtype = request.POST['jobtype']
    user_id = file.owner.id
    si = schedInterface()
    file.save()
    if jobtype in ['minimization']:
        newJobID = si.submitJob(user_id,file.location,scriptlist)
        #updates the parameter object's html
        minparam = minimization.models.minimizeParams.objects.filter(pdb=file,selected='y')[0]
        if newJobID < 0:
            minparam.statusHTML = "<font color=red>Failed</font>"
        else:
            file.minimization_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            minparam.statusHTML = statsDisplay(sstring,newJobID)
            file.save()
        minparam.save()
    elif jobtype in ['solvation']:
        newJobID = si.submitJob(user_id,file.location,scriptlist)
        #updates the parameter object's html
        solvparam = solvation.models.solvationParams.objects.filter(pdb=file,selected='y')[0]
        if newJobID < 0:
            solvparam.statusHTML = "<font color=red>Failed</font>"
	else:
            file.solvation_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            solvparam.statusHTML = statsDisplay(sstring,newJobID)
            file.save()
        solvparam.save()
    elif jobtype in ['nma']:
        newJobID = si.submitJob(user_id,file.location,scriptlist)
        nmparam = normalmodes.models.nmodeParams.objects.filter(pdb=file,selected='y')[0]
        if newJobID < 0:
            nmparam.statusHTML = "<font color=red>Failed</font>"
	else:
            file.nma_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            file.nma_status = statsDisplay(sstring,newJobID)
            nmparam.statusHTML = statsDisplay(sstring,newJobID)
            file.save()
        nmparam.save()
    elif jobtype in ['md']:
        newJobID = si.submitJob(user_id,file.location,scriptlist)
        mdparam = dynamics.models.mdParams.objects.filter(pdb=file,selected='y')[0]
        if newJobID < 0:
            mdparam.statusHTML = "<font color=red>Failed</font>"
	else:
            si = schedInterface()
            file.md_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            mdparam.statusHTML = statsDisplay(sstring,newJobID)
            mdparam.md_movie_status = None
            file.save()
        mdparam.save()
    elif jobtype in ['ld']:
        newJobID = si.submitJob(user_id,file.location,scriptlist)
        ldparam = dynamics.models.ldParams.objects.filter(pdb=file,selected='y')[0]
        if newJobID < 0:
            file.minimization_status = "<font color=red>Failed</font>"
	else:
            file.ld_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            ldparam.statusHTML = statsDisplay(sstring,newJobID)
            ldparam.md_movie_status = None
            file.save()
    elif jobtype in ['sgld']:
        newJobID = si.submitJob(user_id,file.location,scriptlist)
        sgldparam = dynamics.models.sgldParams.objects.filter(pdb=file,selected='y')[0]
        if newJobID < 0:
            file.minimization_status = "<font color=red>Failed</font>"
	else:
            file.sgld_jobID = newJobID
            sstring = si.checkStatus(newJobID)
            sgldparam.statusHTML = statsDisplay(sstring,newJobID)
            sgldparam.md_movie_status = None
            file.save()
    elif jobtype in ['energy']:
	file.lessons.onEnergySubmit(request.POST)
        enerjobID = si.submitJob(user_id,file.location,scriptlist)
        #Append/energy needs to be doen before energy can be parsed
        done = re.compile('Done')
        failed = re.compile('Failed')
        sstring = si.checkStatus(enerjobID)
        while(not done.search(statsDisplay(sstring,enerjobID)) and not failed.search(statsDisplay(sstring,enerjobID))):
            sstring = si.checkStatus(enerjobID)
        file.save()
	file.lessons.onEnergyDone(request.POST)
        return parseEnergy(file,'charmm-' + file.stripDotPDB(file.filename) + '-energy.out')

    
#Saves script changes
def changeScript(request):
    if not account.views.isUserTrustworthy(request.user):
        return False
    file =  PDBFile.objects.filter(owner=request.user,selected='y')[0]
    charmm_filename = request.POST['filename']
    input_text = request.POST['text']
    os.chdir(file.location)
    outfp = open(charmm_filename,'w+')
    outfp.write(input_text)
    outfp.close()
    sanitizeScript(file,charmm_filename)

#Retrieves input data
def getInputData(request):
    if not account.views.isUserTrustworthy(request.user):
        return False
    file = PDBFile.objects.filter(owner=request.user,selected='y')[0]
    file.checkRequestData(request)
    os.chdir(file.location)
    charmm_filename = request.POST['filename']
    infp = open(charmm_filename,'r')
    charmm_inp = infp.read()
    infp.close()
    infp = open(charmm_filename,'r')
    charmm_inplist = infp.readlines()
    infp.close()
    #charmm_inplist is used to find the number of rows for the textarea
    if len(charmm_inplist) > 70:
        rows = 70
    else:
        rows = len(charmm_inplist)
    html = '<textarea id="' + charmm_filename + '" name="' + charmm_filename + '" cols=100 rows=' + str(rows) + ' onBlur="sendScriptChanges(this);">\n'
    html += charmm_inp
    html += '</textarea>\n'
    return HttpResponse(html)
    
#generates the html to edit the script
def generateHTMLScriptEdit(charmm_inp,scriptlist,jobtype):
    html = '<div id="edit_scriptcontainer">\n'
    html += '<form id="editscript_form" method="post" action="." onsubmit='
    html += '"return send_form_editscript(this,\'/charmming/runjob/\',\'edit_scriptcontainer\')"> \n'
    html += "Choose a file to edit: <br> \n"
    last_script = ''
    #hidden line will be a hidden html form element that stores the scriptlist
    hiddenline = ''
    #stripfilepath will strip the file path if associated with a script
    stripfilepath = re.compile('.*/.*/')
    for filename in scriptlist:
        html += '<input type="radio" id=sel_"' + filename + '" name="edit_script"'
	#This will check if the file is the last item in the scriplist
	#Then store it. The last file in the script list should automatically
	#be selected
	if filename == scriptlist[len(scriptlist)-1]:
            infp = open(filename,'r')
            charmm_inplist = infp.readlines()
            infp.close()
            #charmm_inplist is used to find the number of rows for the textarea
            if len(charmm_inplist) > 70:
                rows = 70
            else:
                rows = len(charmm_inplist)
	    hiddenline += filename
	    html += 'value="' + filename + '" onclick="getInputScript(this,\'edit_scriptdiv\');" checked>' + stripfilepath.sub('',filename) + '<br> \n'
        else:
	    hiddenline += filename + ','
	    html += 'value="' + filename + '" onclick="getInputScript(this,\'edit_scriptdiv\');">' + stripfilepath.sub('',filename) + '<br> \n'
    html += '<input type="submit" value="Save files and Submit Job"><br>\n'
    html += '<div id="edit_scriptdiv">'
    html += '<textarea id="' + filename + '" name="' + filename + '" cols=100 rows=' + str(rows) + ' onBlur="sendScriptChanges(this);">\n'
    html += charmm_inp
    html += '</textarea>\n'
    html += '</div>'
    html += '<input type="hidden" id="jobtype" name="jobtype" value="' + jobtype + '"> \n'
    html += '<input type="hidden" name="scriptlist" value="' + hiddenline + '"> \n'
    html += '<input type="submit" value="Save files and Submit Job">\n'
    html += '</form>\n'
    html += '</div>\n'
    return html

#Removes malicious data from scripts
#Streams, opens, write, systems must all be handled strictly
def sanitizeScript(file,script_name):
    os.chdir(file.location)
    script = open(script_name,'r')
    scriptlines = script.readlines()
    script.close()
    script = open(script_name,'w+')
    syst = re.compile('syst', re.IGNORECASE)
    stream = re.compile('stream', re.IGNORECASE)
    name = re.compile('name', re.IGNORECASE)
    openword = re.compile('open', re.IGNORECASE)
    readword = re.compile('read', re.IGNORECASE)
    fileword = re.compile('file', re.IGNORECASE)
    writeword = re.compile('file', re.IGNORECASE)
    acceptable_syst_lines = ['syst "/usr/local/charmming/savechrg.py @TOTCH ' + file.location + file.stripDotPDB(file.filename) + '.charge"\n','syst "/usr/local/charmming/savegv.py @GREATERVALUE @SECONDVALUE rhdo ' + file.location + 'crystl_' + file.stripDotPDB(file.filename) + '.str"\n','syst "/usr/local/charmming/savegv.py @GREATERVALUE @SECONDVALUE cubic ' + file.location + 'crystl_' + file.stripDotPDB(file.filename) + '.str"\n','syst "/usr/local/charmming/savegv.py @GREATERVALUE @SECONDVALUE sphere ' + file.location + 'crystl_' + file.stripDotPDB(file.filename) + '.str"\n','syst "/usr/local/charmming/savegv.py @GREATERVALUE @SECONDVALUE hexa ' + file.location + 'crystl_' + file.stripDotPDB(file.filename) + '.str"\n','syst "/usr/local/charmming/savegv.py @GREATERVALUE @SECONDVALUE @XTALTYPE ' + file.location + 'crystl_' + file.stripDotPDB(file.filename) + '.str"\n','syst "/usr/localcharmming/calcewald.pl @GREATERVALUE ' + file.location + ' ' + file.stripDotPDB(file.filename) +' "\n']
    acceptable_stream_lines = ['stream ' + file.stripDotPDB(file.filename) + '-highnum.str\n','stream crystl_' + file.stripDotPDB(file.filename) + '.str\n','stream /usr/local/charmming/coarsegrain.str\n','stream /usr/local/charmming/enmrestraint-multi.str\n','stream ' + file.location + 'charmm-' + file.stripDotPDB(file.filename) + '-crystl.str\n','stream ' + file.location + 'new_' + file.stripDotPDB(file.filename) + '-addions.str\n']
    acceptable_open_lines = ['/usr/local/charmming/toppar/top_all27_prot_na.rtf','/usr/local/charmming/toppar/par_all27_prot_na.prm']
    for line in scriptlines:
        if syst.search(line):
	    if line not in acceptable_syst_lines and line[0] not in ['!']:
	        line = '!Line below emitted omitted: contains syst function.\n' + '!' + line + '\n'
	if stream.search(line):
	    if line not in acceptable_stream_lines and line[0] not in ['!']:
	        line = '!Line below emitted omitted: contains stream function.\n' + '!' + line + '\n'
	#Handles the "open" string in the scripts
	#Removes a filepath if there is one
	if openword.search(line):
	    if line[0] not in ['!']:
	        open_line = line.strip().split(' ')
		for i in range(len(open_line)-1):
		    if name.search(open_line[i]):
			#filepathindex is the index where the filepath should be
			#It is basically the index after "name"
			filepath_index = i + 1
			#This if test determines if that index even exists
			if(filepath_index >= len(open_line)):
			    #If it turns out that index does not exist, then the command is not
			    #formatted properly
			    filepath = "Path Error"
			else:
		            filepath = open_line[filepath_index]
			break
		    else:
	                filepath = "Path Error"
		     
                #stripfilepath will strip the file path if associated with a script
                stripfilepath = re.compile('\w*/.*/')
		#Checks to see if filepath is valid, it not replace it with a new filename with no filepath
		if stripfilepath.search(filepath) and filepath not in acceptable_open_lines:
		    filepath = stripfilepath.sub('',filepath)
		    open_line[filepath_index] = filepath
		    line = ' '.join(open_line)
	            line = '!Warning! The file path below was removed for security reasons.\n' + line + '\n'
		else:
		    if filepath == "Path Error":
 	                line = '!Warning! The filename should be right after the "name" argument when using open.\n' + line
		    else:
		        line = line
	elif (readword.search(line) or writeword.search(line)) and (fileword.search(line) or name.search(line)) and line[0] not in ['!']:
	    if line[0] not in ['!']:
	        open_line = line.strip().split(' ')
		for i in range(len(open_line)-1):
		    if name.search(open_line[i]) or fileword.search(open_line[i]):
			#filepathindex is the index where the filepath should be
			#It is basically the index after "name"
			filepath_index = i + 1
			#This if test determines if that index even exists
			if(filepath_index >= len(open_line)):
			    #If it turns out that index does not exist, then the command is not
			    #formatted properly
			    filepath = "Path Error"
			else:
		            filepath = open_line[filepath_index]
			break
		    else:
	                filepath = "Path Error"
                #stripfilepath will strip the file path if associated with a script
                stripfilepath = re.compile('\w*/.*/')
		#Checks to see if filepath is valid, it not replace it with a new filename with no filepath
		if stripfilepath.search(filepath) and filepath not in acceptable_open_lines:
		    filepath = stripfilepath.sub('',filepath)
		    open_line[filepath_index] = filepath
		    line = ' '.join(open_line)
	            line = '!Warning! The file path below was removed for security reasons.\n' + line + '\n'
		else:
		    if filepath == "Path Error":
 	                line = '!Warning! The filename should be right after the "name"/"file" argument when using read/write.\n' + line
		    else:
		        line = line

	script.write(line)
    script.close()
