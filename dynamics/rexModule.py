from django.template import *
from django.template.loader import get_template
import models

#pre: Requres request with post data and a dynamics model that has a "replica_exchange" foreign key field (md, ld, sgld objects as of
#     8/10/2009)
#post: Creates a rexParams object and links it to the dynamics model (md, ld, sgld objects as of 8/10/2009)
def createRexObjectFromPostData(request,dynamicsModel):
    #First check to make sure the user wanted replica exchange and this method was not called
    #due to a coding mistake
    try:
        useReplicaExchange = request.POST['use_replica_exchange']
    except:
        return;
    #Obtain number of replicas, how often to merge trajectories, npref value, and number of steps
    #from front-end POST data
    elog = open('/tmp/rlog.txt','w')
    rexmdl = models.rexParams()
    numOfReplicas = request.POST['num_rex']
    rexClean = request.POST['rex_clean']
    rexNpref = request.POST['rex_npref']
    nStep = dynamicsModel.nstep
    
    #Add parameters to rexParams model
    rexmdl.firstStep = 1
    rexmdl.lastStep = nStep
    rexmdl.cleanStep = rexClean
    rexmdl.npref = rexNpref
    rexmdl.nbath = numOfReplicas
    
    #Obtain temperatures from request Data and then store it into model
    rexList = []
    #This for loop goes through through the post data and each temperature has they key of rextemp+a number
    for i in range(1,int(numOfReplicas)+1):
        rexList.append(request.POST['rextemp' + `i`])
    #join list into a string using a space (" ") as a separator
    rexmdl.temperatures = " ".join(rexList)
    rexmdl.save()

    #Add rex model to the dynamics model sent
    dynamicsModel.replica_exchange = rexmdl
    dynamicsModel.save()

#pre: Requires the POST data from the front end which must contain the replica exchange data
#post: Creates a configuration file in the user directory with the filename as pdbfilename-rexconfig.cfg
def createConfigFileFromPostData(request,file,dynamicsModel):
    #First check to make sure the user wanted replica exchange and this method was not called
    #due to a coding mistake
    try:
        useReplicaExchange = request.POST['use_replica_exchange']
    except:
        return;
    #Obtain number of replicas, how often to merge trajectories, npref value, and number of steps
    #from front-end POST data
    numOfReplicas = request.POST['num_rex']
    rexClean = request.POST['rex_clean']
    rexNpref = request.POST['rex_npref']
    nStep = int(dynamicsModel.nstep)

    #Get replica temperatures and store it into the listOfReplicas array/list
    listOfReplicas = []
    for i in range(1,int(numOfReplicas)+1):
        listOfReplicas.append(request.POST['rextemp' + `i`])

    #Write data to config file
    newConfigFile = open(file.location + '/' + file.stripDotPDB(file.filename) + '-rexconfig' + '.cfg','w+')
    newConfigFile.write("FIRST 1 \n") 
    newConfigFile.write("LAST " + `nStep` + "\n")
    newConfigFile.write("CLEAN " + rexClean + "\n")
    newConfigFile.write("NPREF " + rexNpref + "\n")
    newConfigFile.write("NBATH " + numOfReplicas + "\n")
    #Write temperatures to data file
    for replica in listOfReplicas:
        newConfigFile.write(replica + "\n")

    #close config file and fin!
    newConfigFile.close()
    return

    
    
