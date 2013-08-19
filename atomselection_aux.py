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
import selection.models

#Takes a dictionary as input, then fills in Atom/OniomSelection keys as needed, and returns it.
#Centralizes the QM/MM atom selection procedure so that instead of writing 30 lines across n files, we write 1 line across n files and modify it here.
#Also takes a workingstructure as input to make the whole of it a bit easier to handle.
def getAtomSelections(tdict, ws):
    #Also get the atom selection, if any
    atomselection = None
    oniom_selections = None
    try:
        atomselection = selection.models.AtomSelection.objects.filter(workstruct=ws)[0] #Try getting the latest in the query...?
    except:
        pass
    try:
        oniom_selections = selection.models.OniomSelection.objects.filter(workstruct=ws)[0]
        oniom_selections = selection.models.OniomSelection.objects.filter(workstruct=ws)
    except:
        pass
    #Incoming block of generic code
    if oniom_selections != None:
        if atomselection.time_created < oniom_selections[0].time_created: #Timestamps shouldn't change by that much between when each OniomSelection was saved
            #Display the atomselection
            tdict['atomselection'] = atomselection
            lonepairs = selection.models.LonePair.objects.filter(selection=atomselection)
            tdict['lonepairs'] = lonepairs
        else:
            tdict['oniom_selections'] = oniom_selections
            lonepairs = []
            for i in range(0,len(oniom_selections)):
                lonepair_array = selection.models.LonePair.objects.filter(selection=oniom_selections[i])
                for lonepair in lonepair_array:
                    lonepairs.append(lonepair)
            tdict['lonepairs'] = lonepairs #These are generic. We can do post-processing on the form without issue, the DB doesn't have to care.
    else:
        tdict['atomselection'] = atomselection
        lonepairs = selection.models.LonePair.objects.filter(selection=atomselection)
        tdict['lonepairs'] = lonepairs
    #Might be worth something to make these more generic so that when we add more selection types the if statements don't get huge.
    #Now you just force these to render as a django template consisting entirely of javascript on the QM/MM dependent pages and you're done
    #In theory this updates tdict and mails it back out. But just in case...
    return tdict

"""
Takes a request (with attached POST data) as input,
and saves an Atom/OniomSelection to the database according
to the data in POST.
Will force-overwrite any selections already present.
When we implement naming selections, it would be a simple name-equality
check based on the POST data to determine whether we should replace
or just save.
If an error occurs, it will return the error message.
Otherwise, it will return "Done".
Check for "== 'Done'" to check whether it succeeded, and produce error
messages accordingly.
"""
def saveAtomSelections(request,ws,task): #We also need to pass in the working structure for the sake of queries, as well as the task associated to this atom selection
    #Set up qmhosts/mmhosts before we do anything, throw errors if something's wrong.
    save_type = "" #To make less branching
    graphical_select = False
    oldselect = False
    oldoniom = False
    oniom_selections = False
    modelSelected = ""
    if request.POST.has_key('model_selected'):
        modelSelected = request.POST['model_selected']
    else:
        modelSelected = request.POST['modelType'] #Scope creep...
        graphical_select = True
    layers = 0
    if modelSelected == "oniom":
        try:
            layers = int(request.POST['oniom_layers'])
        except:
            try:
                layers = int(request.POST['layers'])
            except:
                return "Invalid number of layers. If you are seeing this message in error, please report it."
    try:
        oldselect = selection.models.AtomSelection.objects.filter(workstruct=ws)[0]
        #There may be more than one atomselection associated to ws (see OniomSelection), so we just check if there are any.
        try: #Now check if there are any OniomSelections
            oldoniom = selection.models.OniomSelection.objects.filter(atomselection_ptr=oldselect.id)[0]
            #Seems so. Let's gather up the OniomSelections.
            oniom_selections = selection.models.OniomSelection.objects.filter(workstruct=ws)
            #Then we overwrite all their fields
            save_type = "old_oniom"
        except: #There are no OniomSelections
            #SO overwrite the atom selection that exists
            save_type = "old_qmmm"
    except: #There are no atom selections present
        #SO create a new one.
        if modelSelected in ["qmmm","oniom"]: #If not, something's very wrong.
            if modelSelected == "oniom": #Serious business. See the procedure in selection.views.
                save_type = "new_oniom"
            else: #Make one and attach fields to it.
                save_type = "new_qmmm"
        else:
            return "Bad model type. Please use one of the pre-set models. If you are seeing this message in error, please report it."
    selection_to_save = None
    if save_type == "old_oniom":
        #Delete the OniomSelections; we don't care about the old ones because the layers can be different. selection.views overwrites and it just makes things worse.
        if oniom_selections: #It should be set in the try block or it wouldn't get here
            oniom_selections.delete()
    elif save_type == "old_qmmm":
        #Delete the selection
        if oldselect: #You shouldn't get here if it's still False.
            oldselect.delete() #Just wipe it and make a new one. More generic that way.
    curr_return = True
    if modelSelected == "qmmm":
        qqh_count = 1
        curr_return = saveSelection(ws,task,modelSelected,request.POST,0,graphical_select,qqh_count)
    elif modelSelected == "oniom":
        qqh_count = 1 #TO update the qqh field in the database
        for current_layer in range(1,layers+1):
            curr_return = saveSelection(ws,task,modelSelected,request.POST,current_layer,graphical_select,qqh_count)
            try:
                qqh_count = int(curr_return)
            except:
                qqh_count = 1
                return curr_return
    else: #Oh dear.
        return "Bad save_type. Please report this message!"
    return curr_return


"""
Takes as input:
    WorkingSTructure
    Task
    Save_type (QM/MM or ONIOM)
    POST data
    current layer (if ONIOM), or 0
    GraphicalSelect (if from graphical selection page or False otherwise)
Attaches this data to an atom selection.
Makes this procedure nice and generic, avoiding issues with code overflow, and allows
more models to be easily built.
Returns True if successful, otherwise an error message.
"""
def saveSelection(ws,task,save_type,data,layer,graphical_select,qqh_count):
    model_string = ""
    total_layers = 1
    isQM = True
    if save_type == "oniom":
        model_string = "_layer_"+str(layer)
        if graphical_select:
            total_layers = int(data['layers'])
        else:
            total_layers = int(data['oniom_layers'])
        atomselect = selection.models.OniomSelection()
        try:
            if int(data['highest_qm_layer']) >= layer: #TODO: Modify this for QM/MM/QM calc...
                atomselect.isQM = isQM
            else:
                isQM = False
        except:
            return "highest_qm_layer is not an integer. Please report this message!"
        atomselect.total_layers = total_layers
        atomselect.layer_num = layer
    elif save_type == "qmmm": #We already checked for any other values, but best to use elif for future expansion
        atomselect = selection.models.AtomSelection()
    postdata_string = ""
    if not graphical_select:
        postdata_string = "qmmm_"
    atomselect.workstruct = ws
    atomselect.task = task
    atomselect.selection_type = save_type
    if layer < total_layers:
        if graphical_select:
            atomselect.selectionstring = data['atomselection'+model_string]
        else:
            atomselect.selectionstring = data['qmsele'+model_string]
    atomselect.save()
    logfp = open("/tmp/atomselect_test.txt","w")
    if isQM:
        atomselect.exchange = data[postdata_string+'exchange'+model_string]
        atomselect.charge = data[postdata_string+'charge'+model_string]
        atomselect.correlation = data[postdata_string+'correlation'+model_string]
        if graphical_select:
            atomselect.basis_set = data['basis_set'+model_string]
        else:
            atomselect.basis_set = data[postdata_string+'basisset'+model_string]
        atomselect.multiplicity = data[postdata_string+'multiplicity'+model_string]
        qmsearchstring = ""
        #Here we need to branch again because the two pages have different formats for keeping track of things.
        if graphical_select: #The difference in the search strings is my own fault - I did not expect the scope creep. But this is the standard, not the exception, here in CHARMMing... ~VS
            atomselect.linkatom_num = data['linkatom_num'+model_string]
            qmsearchstring = "qmhost"+model_string+"_" #e.g. linkqm_layer_1_
        else:
            atomselect.linkatom_num = data['num_linkatoms'+model_string]
            qmsearchstring = "linkqm"+model_string+"_"
        atomselect.save()
        if int(atomselect.linkatom_num) > 0:
            mmsearchstring = qmsearchstring.replace("qm","mm") #Should work alright, since "layer" and the ints do not contain qm. Watch for this for future models!
            woof = data[qmsearchstring+"0"]
            for key in data.keys():
                #Much more efficient than qmhosts/mmhosts
                if key.startswith(qmsearchstring):
                    divid = key.split("_")[-1]
                    if save_type == "oniom":
                        divid = "_layer_"+str(layer)+"_"+str(divid)
                    lonepair = selection.models.LonePair()
                    lonepair.selection = atomselect
                    lonepair.divid = divid
                    if graphical_select:
                        qm_data = data[key].split("\t")
                        lonepair.qmresid = qm_data[0]
                        lonepair.qmatomname = qm_data[1]
                        lonepair.qmsegid = qm_data[2]
                        mm_data = data[mmsearchstring+key.split("_")[-1]].split("\t") #e.g. linkmm_layer_1_0, the matching one to the linkqm
                        lonepair.mmresid = mm_data[0]
                        lonepair.mmatomname = mm_data[1]
                        lonepair.mmsegid = mm_data[2]
                        lonepair.save()
                    else:
                        lonepair.qmresid = data[key]
                        lonepair.qmatomname = data[key.replace("linkqm","qmatomtype")]
                        lonepair.qmsegid = data[key.replace("linkqm","linkqmsegid")]
                        lonepair.mmresid = data[key]
                        lonepair.mmatomname = data[key.replace("linkqm","mmatomtype")]
                        lonepair.mmsegid = data[key.replace("linkqm","linkmmsegid")]
                        lonepair.save()
                    lonepair.qqh = qqh_count #This way they're arranged by layer, and have labels we can pull on easily
                    qqh_count = qqh_count + 1
                    lonepair.save()
    atomselect.save()
    return qqh_count
"""
Validates a set of input data, spits out an error if something wrong is found, otherwise returns True
Layers is there for ONIOM simulations.
The error spat out is a string, which later gets turned into a messages.error call.
This input is validated BEFORE being fed in to the saveAtomSelections function, to make error reporting easy.
"""
def validateInputs(tdict,postdata,layers):
    specialchars = set("#$/;\n\\_+=[]{}()&^%")
    try:
        highest_qm_layer = int(postdata['highest_qm_layer'])
    except:
        return "Highest QM layer is not an integer."
    try:
        if layers > 1:
            for curr_layer in range(1,layers): #Exclusive end of range...last layer has no selection!
                foo = postdata['atomselection_layer_'+str(curr_layer)]
                if not(foo.startswith("bynum")) or (len(specialchars.intersection(foo)) > 0): #i.e. no special chars
                    return "Invalid selection string. Please use the automated selection system."
        else:
            foo = postdata['atomselection']
            if not(foo.startswith("bynum")) or (len(specialchars.intersection(foo)) > 0): #i.e. no special chars
                return "Invalid selection string. Please use the automated selection system."
    except: #I'm not sure what would break here, but let's do it anyway.
        return  "You have selected too many atoms, or performed an invalid selection. Please try again."
    try:
        if layers > 1 and highest_qm_layer >= 1:
            for curr_layer in range(1,layers):
                tdict['num_linkatoms_layer_'+str(curr_layer)] = int(postdata['linkatom_num_layer_'+str(curr_layer)])
        else:
            tdict['num_linkatoms'] = int(postdata['linkatom_num']) #This will fail if the user put in anything other than a number because of how int works.
    except:
        return "Number of link atoms is not an integer."
    return "Passed"
