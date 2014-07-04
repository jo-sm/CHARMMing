import selection.models
import copy, os, shutil, re #You can never be too safe in copying objects!
import output
import charmming_config
from django.template.loader import get_template
from django.template import *
import structure.qmmm
from structure.models import WorkingFile


"""
Subsystem objects have the following attributes:
    name: name of the subsystem (i.e. reg[0-9]+)
    coef: +/-1, depending on position in the equation, default 1
    out_file: Path to the output file for this subsystem (default: /home/schedd/user/structure/subsystems/; replace as needed)
    selection_string: atom selection associated to this subsystem
"""
class Subsystem:
    name = str()
    coef = float()
    ##out_file = charmming_config.user_home+ "/" +  + "/structure/subsys/"
    ##inp_file = charmming_config.user_home+"/structure/subsys/"
    layer = int() #This notes what layer it's associated to, for the sake of getting the right PSF/CRDs for it.
    selection_string = str()
    deletion_string = str() #Hack because reg1 etc. doesn't get transferred to the subsystem input files...
    level_of_theory = 0 #Fix
    lonepairs = []
    
    def __init__(self, nomen, coefficient, output, input, sele, lot, lay,lp):
        self.name = nomen
        self.coef = coefficient
        self.out_file = output
        self.inp_file = input
        self.selection_string = sele
        self.level_of_theory = lot
        self.layer = lay
        self.lonepairs = lp

    def __init__(self,username):
        self.name = ""
        self.coef = float()
        self.out_file = charmming_config.user_home + "/" + username + "/structure/subsys/"
        self.inp_file = charmming_config.user_home + "/" + username + "/structure/subsys/"
        self.selection_string = ""
        self.level_of_theory = 0
        self.layer = 1
        self.lonepairs = []

"""
Takes a list of OniomSelections as input.
Returns a list with the following structure:
    [(OniomSelection,[LonePair])]
i.e., a list of tuples, where the first half of the tuple is an OniomSelection
and the second half is the list of LonePair objects (if any) associated to that selection.
"""
def getLonePairs(oniom_selections):
    atomselections = []
    oniom_selections = oniom_selections.order_by("layer_num") #Sort by layer for easier processing later
    for current_selection in oniom_selections:
        selection_list = []
        selection_list.append(current_selection)
        current_lonepairs = []
        selection_list.append(current_lonepairs) #After this, we fill current_lonepairs up with things
        #The list at this point looks like this:
        # [OniomSelection, []]
        #In updating current_lonepairs, we update that to become:
        # [OniomSelection, [LonePair]]
        lonepairs = selection.models.LonePair.objects.filter(selection=current_selection)
        for lonepair in lonepairs: #If the array is empty, it doesn't matter
            current_lonepairs.append(lonepair) #This updates current_lonepairs, such that the array above gets filled.
        atomselections.append(selection_list) #This appends the tuple.
    return atomselections #We're done.

"""
Takes a list of OniomSelections as input, as well as a user's username and the name of the currently selected structure (to reduce DB queries in this method)
Returns subsystems, an array of Subsystem objects.
"""
def generateSubsystems(oniom_selections, username, structname):
    subsystems = []
    total_layers = oniom_selections[0].total_layers #This value is the same across the list, saves arguments for the function since this tells us how many subsystems
    num_subsystems = (2 * total_layers) - 1 #any MSCALE has 2n - 1 terms in the equation.
    oniom_selections = sorted(oniom_selections,key=lambda x:x.layer_num) #Sort for easier work
    #Region number is sorted by layer.
    for current_layer in reversed(xrange(1,total_layers+1)):
        #STart at layer of lowest level of theory (since that one gets added once and doesn't get subtracted by another level)
        #Add terms in pairs after.
        if current_layer == total_layers: #"Whole system" level.
            sub = Subsystem(username)
            sub.name = "alllow"
            sub.coef = 1.00
            sub.layer = current_layer
            sub.selection_string = ".not. ( type qqh* )" #Everything, but no QM link atoms. Irrelevant if this layer is MM or not.
            sub.out_file = sub.out_file.replace("structure",structname).replace("user",username) + sub.name + ".out"
            sub.inp_file = sub.out_file.replace(".out",".inp") #This is NOT a regular expression. No worries about the ., it's interpreted as a literal
            sub.level_of_theory = (total_layers - current_layer) + 1 #This will always equal 1. It is made with these variables for clarifying the pattern used in the code.
            subsystems.append(sub)
        else:
            #I am aware that this procedure is repetitive, but we need to take different standards for different pieces of the MSCALE equation.
            sub1 = Subsystem(username)
            sub1.name = "r"+str(current_layer)+"_plus" #For a fully generic experience, we need to just make these +/- since once we go past 3 layers, the terms "high/med/low" have no meaning anymore
            sub1.coef = 1.00
            sub1.layer = current_layer
            sub1.selection_string = "reg"+str(current_layer)
            sub1.deletion_string = oniom_selections[current_layer - 1].selectionstring
            sub_lonepairs = selection.models.LonePair.objects.filter(selection=oniom_selections[current_layer -1 ]) #oniom_selections are 0-indexed
            for lonepair in sub_lonepairs:
                sub1.selection_string = sub1.selection_string + " .or. type qqh" + str(lonepair.qqh)
                sub1.lonepairs.append(lonepair)
            sub1.out_file = sub1.out_file.replace("structure",structname).replace("user",username) + sub1.name + ".out"
            sub1.inp_file = sub1.out_file.replace(".out",".inp")
            sub1.level_of_theory = (total_layers - current_layer) + 1
            sub2 = Subsystem(username)
            sub2.name = sub1.name.replace("plus","minus")
            sub2.coef = -1.00
            sub2.layer = current_layer #Both of these live in the same place
            sub2.selection_string = copy.deepcopy(sub1.selection_string) #You can never be too careful.
            sub2.deletion_string = copy.deepcopy(sub1.deletion_string)
            sub2.out_file = sub1.out_file.replace("_plus","_minus") #Might be a problem with people with _plus in their names. Please check!
            sub2.inp_file = sub2.out_file.replace(".out",".inp")
            sub2.level_of_theory = (total_layers - current_layer)
            sub2.lonepairs = copy.deepcopy(sub1.lonepairs)
            subsystems.append(sub1)
            subsystems.append(sub2)
    return subsystems #This can now be used by django for the template.

"""
Makes the CHARMM inp files and directories for an MSCALE calculation,
then saves the results on the task and the dictionary.
"""
def make_mscale(template_dict, request, modelType, task):
    if modelType == "oniom":
        template_dict['mscale_path'] = charmming_config.charmm_mscale_exe
        oniom_selections = selection.models.OniomSelection.objects.filter(workstruct=task.workstruct) #TODO: Modify this for implementing name functionality
        atomselections = getLonePairs(oniom_selections)
        subsystems = generateSubsystems(oniom_selections,request.user.username,task.workstruct.structure.name)
        template_dict['nregion'] = oniom_selections[0].total_layers
        template_dict['num_subsystems'] = len(subsystems)
        template_dict['subsystems'] = subsystems
        template_dict['atomselections'] = atomselections
        template_dict['psf_path'] = task.workstruct.structure.location + "/subsys/"

        path_to_make = task.workstruct.structure.location + "/subsys/"
        if not os.path.isdir(path_to_make):
            os.mkdir(path_to_make) #Check permissions...

        os.chmod(path_to_make,0775)
        template_dict['atomselections_reversed'] = oniom_selections.order_by('-layer_num')
        #Now generate QChem scripts (see qmmm.py)
        os.chdir(task.workstruct.structure.location)
        #Normally we would place these in qcstuff/ but we can't because Q-Chem is really strange about paths.
        psf_path = template_dict['psf_path']
        template_dict['psf_path'] = template_dict['psf_path'].replace("subsys/","")
        current_level_of_theory = 1
        reversed_selections = oniom_selections.order_by('-layer_num') #Reverse so that the lowest level of theory is assigned to the layer with the "highest" number, e.g. layer 4 in a 4-layer system
        lowest_qm_level_of_theory = 0 #Once we hit a QM layer, evrerything below it is a QM level of theory (for now)
        current_layer = 0
        basepath = ""
        task.qmlayers = reversed_selections[0].layer_num #Very important for saving qchem inp/out files!
        task.save()
        while current_layer < reversed_selections[0].layer_num: #Stop at the max...
            if reversed_selections[current_layer].isQM:
                if lowest_qm_level_of_theory < 1:
                    lowest_qm_level_of_theory = current_level_of_theory #This makes the subsystem input scripts easier. Only need to check once, all layers "after" a QM layer are QM
                qmparams = structure.qmmm.makeQchem_val("oniom",reversed_selections[current_layer],layer=current_layer) #Since we start at, say, 4 and go down, this works out.
                qmparams['jobtype'] = 'Force' #????
                template_dict = structure.qmmm.makeQChem_tpl(template_dict,qmparams,"level-"+str(current_level_of_theory),task)  #This would make the level correct
                #We make input file WorkingFiles so that the user can see them...then in the finish() we fetch the .out versions.

                basepath = task.workstruct.structure.location +'/'
                path = basepath + "level-"+str(current_level_of_theory)+".inp"

            current_level_of_theory += 1
            current_layer += 1
        #Now generate subsystem scripts
        template_dict['lowest_qm_level_of_theory'] = lowest_qm_level_of_theory
        os.chdir("subsys") #Change to the right dir
        for subsystem in subsystems:
            if subsystem.level_of_theory < lowest_qm_level_of_theory:
                template_dict['QM_subsystem'] = False
            else:
                template_dict['QM_subsystem'] = True
            template_dict['sub'] = subsystem
            templ = get_template('%s/mytemplates/input_scripts/mscale_subsystem.inp' % charmming_config.charmming_root)
            charmm_input = output.tidyInp(templ.render(Context(template_dict))) #Since we already have atomselections and the other stuff, this should be easy.
            path = basepath + "subsys/" + subsystem.name + ".inp"
            inp_out = open(path,"w")
            inp_out.write(charmm_input)
            inp_out.close()
        template_dict['psf_path'] = psf_path #Just in case we need it again somewhere...

        # I have no idea if this is the right place to do things, but let's find out...
        template_dict['linkat_list'] = []
        template_dict['bynum_list'] = []
        template_dict['region_list'] = []

        for atomselection in oniom_selections:

            lonepairs = selection.models.LonePair.objects.filter(selection=atomselection)
            for lonepair in lonepairs:
                ldict = {}
                ldict['segq'] = lonepair.qmsegid
                ldict['segm'] = lonepair.mmsegid
                ldict['resq'] = lonepair.qmresid
                ldict['resm'] = lonepair.mmresid
                ldict['atypq'] = lonepair.qmatomname
                ldict['atypm'] = lonepair.mmatomname
                template_dict['linkat_list'].append(ldict)

            selar = atomselection.selectionstring.replace('.or.',"").split('bynum')
            for bnp in selar:
                bnp = bnp.strip()
                if bnp:
                    bynumdict = {}
                    # ToDo -- need some error checking in here
                    bynumdict['numa'] = bnp.split(':')[0]
                    try:
                        bynumdict['numb'] = bnp.split(':')[1]
                    except IndexError, e:
                        bynumdict['numb'] = -9999
                    template_dict['region_list'].append(atomselection.layer_num)
                    template_dict['bynum_list'].append(bynumdict)

        logfp = open('/tmp/mscale_renumber.txt','w')
        logfp.write('bynum_list: %s\n' % template_dict['bynum_list'])
        logfp.write('region_list: %s\n' % template_dict['region_list'])
        logfp.write('linkat_list: %s\n' % template_dict['linkat_list'])
        logfp.close()

        task.save()
        return template_dict
