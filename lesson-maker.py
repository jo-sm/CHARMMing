#!/usr/bin/python

#makes lessons in python
import sys
import os
import re
import MySQLdb
import settings

try:
    lesson_num = sys.argv[1]
except:
    lesson_num = ""

#If there are files that accompany the lesson, the developer should supply the location of such files 
#and the names of the files that will be used
if lesson_num == "" or lesson_num == "-h":
    print "Use this tool to create lessons.\n "
    print "python lesson-maker <<new_lesson_number>> <<location_of_files_the_user_must_download>> [name_of_files_user_must_download_in_brackets_seperated_by_a_comma]\n"
    sys.exit()

lesson_name = "lesson" + lesson_num

#Check and see if they supplied a directory where the files the user must download is located
try:
    lesson_file_location = sys.argv[2]
except:
    lesson_file_location = ""
if lesson_file_location == "":
    print "No file location was specified, assuming filenames are in current working directory\n"
    lesson_file_location = "."

#Check and see if the developer supplied any files the lesson user must download to continue with the elsson
try:
    lesson_filenames = sys.argv[3]
except:
    lesson_filenames = ""
if lesson_filenames == "":
    print "No files were given that should be used in this lesson. It is assumed that in this lesson, the user will not have to download any other files to upload.\n"
    print "If you would like to add files later, store them in mytemplates/lessons/" + lesson_name + "/\n"
else:
    print "Copying files that accompany lesson into mytemplates/lessons/" + lesson_name
    lesson_filenames = lesson_filenames.strip('[]')
    lesson_filenameslist = lesson_filenames.split(',')
    try:
        os.stat('mytemplates/lessons/' + lesson_name)
    except:
	os.mkdir('mytemplates/lessons/' + lesson_name)
    for filename in lesson_filenameslist:
        os.system("cp " + lesson_file_location + "/" + filename + " mytemplates/lessons/" + lesson_name + "/")
    

os.system("python manage.py startapp " + lesson_name)

print "Creating models.py ...\n"
modelfile = open(lesson_name + "/models.py",'w+')
modelfile_text = """
# lesson """ + lesson_num + """
# Enter Descriptions Here
from django.db import models
from django.contrib.auth.models import User
from lessons.models import LessonProblem
from solvation.models import solvationParams
from minimization.models import minimizeParams
from dynamics.models import mdParams, ldParams, sgldParams
from normalmodes.models import nmodeParams
import os
import pdbinfo
import lessonaux
import re

class Lesson""" + lesson_num + """(models.Model):
    class Admin:
        pass
    user = models.ForeignKey(User)
    #Please set the nstep value to the number of steps the lesson will have
    nSteps = models.PositiveIntegerField(default=4)
    curStep = models.FloatField(default=0,decimal_places=1,max_digits=3)
    
    #Fill in the below data. If you are not going to use one of these triggers
    #then just leave it blank

    def onFileUpload(self,postdata):
        return True

    def onEditPDBInfo(self,postdata):
        return True

    def onMinimizeSubmit(self,postdata):
        return True

    def onMinimizeDone(self,file):
        return True

    def onSolvationSubmit(self,postdata):
        return True

    def onSolvationDone(self,file):
        return True

    def onEnergySubmit(self,postdata):
        return True

    def onEnergyDone(self,file):
        return True

    def onNMASubmit(self,postdata):
        return True

    def onNMADone(self,file):
        return True

    def onMDSubmit(self,postdata):
        return True

    def onMDDone(self,file):
        return True

    def onLDSubmit(self,postdata):
        return True

    def onLDDone(self,file):
        return True

    def onSGLDSubmit(self,postdata):
        return True

    def onSGLDDone(self,file):
        return True
    
    #The code below is used to display lesson status. It generates the HTML for the lessons
    #The following is how the curstep should be numbered to remain consistent with the code:
    #An integer means the job is complete , for example 1 means step 1 was complete, 2 means step 2 was complete
    #A 0.5 means step 1 is running, a 1.5 means step 2 is running
    #If there is a lesson problem that has the same step number as the current step, that means there was a
    #user error in following the lessons and so the lesson step has failed
    #The below code can be used as an example
    def generateStatusHtml(self,file):
        step_status_list = []
        #First generate the html for each step and a short description 
        step_status_list.append("<tr class='status'><td class='status'>1. File Uploaded: ")
        step_status_list.append("<tr class='status'><td class='status'>2. Solvation: ")
        step_status_list.append("<tr class='status'><td class='status'>3. Minimization: ") 
        step_status_list.append("<tr class='status'><td class='status'>4. MD: ")
        
	#Check to see if a lessonproblem exists, and if so store it into lessonprob
	try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson""" + lesson_num + """',lesson_id=self.id)[0]
        except:
            lessonprob = None
	#Go through the step and see which are done,failed, still running. You can leave the code below this
	#comment untouched. Only edit it if you know what you are doing.
        for i in range(self.nSteps):
            if lessonprob and lessonprob.errorstep == (self.curStep+1) and self.curStep == i:
                step_status_list[i] += ("<font color='red'>Failed</font></td></tr>")
                continue
            elif (self.curStep-.5) == i and self.curStep%1 == 0.5:
                step_status_list[i] += ("<font color='blue'>Running</font></td></tr>")
                continue 
            elif i < self.curStep or len(step_status_list)+1 == self.curStep:
                step_status_list[i] += ("<font color='green'>Done</font></td></tr>")
                continue
            elif i+1 > self.curStep :
                step_status_list[i] += ("<font color='grey'>N/A</font></td></tr>")
                continue

        return step_status_list
    
    #This is used to display the lesson page. It tells the template to display which set of directions
    #for the next step
    #Each step is represented as an index in an array. The corresponding value tells the template
    #what is currently happening in regards to the lesson.
    #If the corresponding value is 2, the job is running
    #If it is 1 then the job is finished
    #0 means the job has not stared
    #-1 means there was an error
    #Only modify this code if your lesson will have greater than 4 steps...All you have to do if you want
    #to add more steps is to continue the if statements. For example if you want a lesson with 5 steps
    #just add the following:
    #    if self.curStep > 4:
    #        if self.curStep == 4.5:
    #            htmlcode_list[3] = 2
    #        else:
    #            htmlcode_list[3] = 1

    def getHtmlStepList(self):
        htmlcode_list = []
        for step in range(self.nSteps):
            htmlcode_list.append(0)
        if self.curStep > 0: 
            htmlcode_list[0] = 1
        if self.curStep > 1:
            if self.curStep == 1.5:
                htmlcode_list[1] = 2
            else:
                htmlcode_list[1] = 1
        if self.curStep > 2:
            if self.curStep == 2.5:
                htmlcode_list[2] = 2
            else:
                htmlcode_list[2] = 1
        if self.curStep > 3:
            if self.curStep == 3.5:
                htmlcode_list[3] = 2
            else:
                htmlcode_list[3] = 1
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson1',lesson_id=self.id)[0]
            htmlcode_list[lessonprob.errorstep-1] = -1
        except:
            lessonprob = None
        return htmlcode_list

"""
modelfile.write(modelfile_text)
modelfile.close()

print "Creating views.py ...\n"
viewsfile = open(lesson_name + "/views.py",'w+')
viewsfile_text = """
from django import forms
from django.http import HttpResponseRedirect, HttpResponse
from django.shortcuts import render_to_response
from account.views import isUserTrustworthy
from pdbinfo.models import PDBFile, PDBFileForm
from lessons.models import LessonProblem
from lesson""" + lesson_num + """.models import Lesson""" + lesson_num + """
from django.contrib.auth.models import User
from django.core import validators
from django import newforms as forms
from django.template import *
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
import re
import copy
import os

#Displays lesson """ + lesson_num + """ page
def lesson""" + lesson_num + """Display(request):
    PDBFile().checkRequestData(request)
    try:
        file = PDBFile.objects.filter(owner=request.user,selected='y')[0]
    except:
        return render_to_response('html/lesson""" + lesson_num + """.html')
    #If its a lesson""" + lesson_num + """ object, get the id by the file id
    if file.lesson_type == 'lesson""" + lesson_num + """':
        lesson_obj = Lesson""" + lesson_num + """.objects.filter(user=request.user,id=file.lesson_id)[0]
        html_step_list = lesson_obj.getHtmlStepList()
    else:
        lesson_obj = None
        html_step_list = None
    try:
        lessonprob_obj = LessonProblem.objects.filter(lesson_type='lesson""" + lesson_num + """',lesson_id=lesson_obj.id)[0]
    except:
        lessonprob_obj = None
    return render_to_response('html/lesson""" + lesson_num + """.html',{'lesson""" + lesson_num + """':lesson_obj,'lessonproblem':lessonprob_obj,'html_step_list':html_step_list})

"""
viewsfile.write(viewsfile_text)
viewsfile.close()

print "Creating mytemplates/html/" + lesson_name + ".html ...\n"
htmlfile = open('mytemplates/html/' + lesson_name + '.html','w+')
htmlfile_text = """
<center>
<h1>Lesson """ + lesson_num + """{% ifequal lesson""" + lesson_num + """.curStep 4 %}: <font color="red">Completed</font> {% endifequal %}</h1>
</center>

Lesson Objectives<br>
<li> List objectives here </li>
<br>
<p> Insert introductory text here </p>"""

if lesson_filenames != "":
    for filename in lesson_filenameslist:
        htmlfile_text += '<li><a href="/charmming/lessons/download/lesson' + lesson_num + '/' + filename + '">' + filename +'</a></li>\n'
    
htmlfile_text += """
<br>
<br>
<!-- 
The way this status display logic works is as follows:
html_step_list is a list/array where each index represents a different step in the lesson
if the value at each index (an index represents a step) is 1, then the step is complete
If the value at each index/step is 2 then the step is currently running
-->

{% if """ + lesson_name + """ %}
 {% ifequal html_step_list.0 1 %}
 <font color="red"> Step 1: Upload Done </font><br>
  <p>
 Insert text here
 </p>
{% endifequal %}

 {% ifequal html_step_list.1 2 %}

 <font color="red"> Step 2: Solvation Running </font><br>
 <p> Insert status text here </p>
 {% endifequal %}
 {% ifequal html_step_list.1 1 %}

 <font color="red"> Step 2: Solvation Done </font><br>
 <p>
 Insert status text here
 </p>
 {% endifequal %}

 {% ifequal html_step_list.2 2 %}
 <font color="red"> Step 3: Minimization Running </font><br>
 <p>
 Insert status text here
 </p>
 {% endifequal %}
 {% ifequal html_step_list.1 1 %}
 <font color="red"> Step 2: Solvation Done </font><br>
 <p>
 Insert text here
 </p>
 {% endifequal %}

<!-- This is the lesson problem object -->
 {% if lessonproblem %}
 <font color="red">LESSON ERROR: STEP - {{lessonproblem.errorstep}}, SEVERITY - {{lessonproblem.severity}}</font> <br>
 {{lessonproblem.description}}
 {% endif %}
{% endif %}

"""
htmlfile.write(htmlfile_text)
htmlfile.close()

print "Installing lesson into settings.py\n"
settingsfile = open('settings.py','r')
settingsfile2 = open('settings2.py','w+')
installedapps = re.compile('INSTALLED_APPS')
endparanth = re.compile('\)')
instapp_init = 0
instapp_text = ""
instapp_firstline = 0
for line in settingsfile:
    if installedapps.search(line):
        instapp_init = 1
        instapp_firstline = 1
        settingsfile2.write(line)
        continue
    if instapp_init and endparanth.search(line):
        settingsfile2.write("    '" + lesson_name + "',\n")
        settingsfile2.write(line)
	continue
    settingsfile2.write(line)
settingsfile2.close()
settingsfile.close()


print "Installing lesson into mytemplates/html/skeleton.html\n"
htmlfile = open('mytemplates/html/skeleton.html','r')
htmlfile2 = open('mytemplates/html/skeleton2.html','w+')
span = re.compile('<span>')
structure = re.compile('Structure')
account = re.compile('Account')
lessons = re.compile('Lessons')
enddiv = re.compile('</div>')
structure_exists = 0
lesson_exists = 0


for line in htmlfile:
    #the user has a lessons section and hasn't removed it
    if span.search(line) and lessons.search(line):
        lesson_exists = 1
	htmlfile2.write(line)
        continue
    #If the user does not have a lessons section, make one right before the accounts section
    if span.search(line) and account.search(line) and lesson_exists < 1:
        htmlfile2.write('<div>\n')
        htmlfile2.write('<span>Lessons</span>\n')
	htmlfile2.write('    <a title="Put lesson description here" href="javascript:openNewURL(\'/charmming/lessons/' + lesson_name + '/\', \'Lesson ' + lesson_num + '\');timeoutWikipage(\'Lesson ' + lesson_num + '\');">Lesson ' + lesson_num + '</a>\n')
        htmlfile2.write('</div>\n')
	htmlfile2.write(line)
	continue
    #if they do have a lessons section then write the url right before the </div> in the lessons section
    if enddiv.search(line) and lesson_exists == 1:
	htmlfile2.write('    <a title="Put lesson description here" href="javascript:openNewURL(\'/charmming/lessons/' + lesson_name + '/\',\'Lesson ' + lesson_num + '\');timeoutWikipage(\'Lesson ' + lesson_num + '\');">Lesson ' + lesson_num + '</a>\n')
	htmlfile2.write(line)
	lesson_exists = 2
	continue
    htmlfile2.write(line)

htmlfile2.close()
        
      

print "Adding the lesson into the Submit Structure page... mytemplates/html/fileuploadform.html \n"
fileuploadfile = open('mytemplates/html/fileuploadform.html','r')
fileuploadfile2 = open('mytemplates/html/fileuploadform2.html','w+')
whatlesson = re.compile("What Lesson is this PDB")
endselect = re.compile("/select")
submitstructure = re.compile("Submit Structure")
onlesson = 0
fileuploadfile_text = ""

for line in fileuploadfile:
    if whatlesson.search(line):
        onlesson = 1
        fileuploadfile2.write(line)
	continue
    #Puts the lesson at the end of the list on fileuploadform
    if endselect.search(line) and onlesson:
        fileuploadfile2.write('         <option value="' + lesson_name + '"> Lesson ' + lesson_num + ' </option>\n')
	fileuploadfile2.write(line)
	continue
    #This checks if the user removed the fileuploadform lesson box
    if submitstructure and onlesson == 0:
        fileuploadform2_text = '''
        <select name="lesson" id="lesson"> 
         <option value="nolesson" selected> No Lesson </option>
         <option value="''' + lesson_name + '''"> Lesson ''' + lesson_num + ''' </option>
        </select><br><br>
	'''
	fileuploadfile2.write(line)
	continue
    fileuploadfile2.write(line)
       
print "Adding lesson into lessonwrapper located in ./lessonaux.py doLessonAct\n"
lessonaux = open('lessonaux.py','r')
lessonaux2 = open('lessonaux2.py','w+')

importpdbinfo = re.compile('import pdbinfo')
onfileupload = re.compile("function == 'onFileUpload'")
on_dolessonact = 0
for line in lessonaux:
    #import the lessons
    if importpdbinfo.search(line):
        lessonaux2.write("import " + lesson_name + '\n') 
    #add the logic for the lesson code
    if onfileupload.search(line):
        lessonaux2.write("    elif lessontype == 'lesson" + lesson_num + "':\n")
	lessonaux2.write("        lesson_obj = lesson" + lesson_num + ".models.Lesson" + lesson_num + ".objects.filter(id=lessonid)[0]\n")
    lessonaux2.write(line)

lessonaux.close()
lessonaux2.close()
	


print "Adding lesson into the lessons/urls.py file\n"

urls = open('lessons/urls.py','r')
urls2 = open('lessons/urls2.py','w+')
admin = re.compile('admin/')
for line in urls:
    if admin.search(line):
        urls2.write("     (r'^lesson" + lesson_num + "/', include('lesson" + lesson_num + ".urls')),\n")    
    urls2.write(line)

urls.close()
urls2.close()





print "Creating custom urls.py for lesson. " + lesson_name + "/urls.py\n"

urls = open(lesson_name + "/urls.py",'w+')
urls_text = """
from django.template import Context, loader
from django.contrib.auth.views import login, logout
from django.conf.urls.defaults import *

urlpatterns = patterns('',
     (r'^$', 'lesson""" + lesson_num + """.views.lesson""" + lesson_num + """Display'),
)
"""

urls.write(urls_text)
urls.close()

print "Installing lesson into pdbinfo/views.py\n"
viewsfile = open('pdbinfo/views.py','r')
viewsfile2 = open('pdbinfo/views2.py','w+')
importos = re.compile('import os')
handlecrdpsf = re.compile('def handleCrdPsf')
editpdbinfo = re.compile('editpdbinfo')
newupload = re.compile('newupload')
elseword = re.compile('else')
on_newupload = 0
on_handlecrdpsf = 0
for line in viewsfile:
    if importos.match(line):
        viewsfile2.write('import ' + lesson_name + '\n')
	viewsfile2.write(line)
	continue
    if handlecrdpsf.search(line):
        on_handlecrdpsf = 1
	viewsfile2.write(line)
	continue
    #adds it into the handlecrdpsf area
    if on_handlecrdpsf == 1 and editpdbinfo.search(line):
        viewsfile2.write("    elif file.lesson_type == 'lesson" + lesson_num + "':\n")
        viewsfile2.write("        lesson_obj = lesson" + lesson_num + ".models.Lesson" + lesson_num + ".objects.filter(user = file.owner,id=file.lesson_id)[0].onFileUpload(request.POST)\n")
	viewsfile2.write(line)
	on_handlecrdpsf = 0
	continue
    if newupload.search(line):
        on_newupload = 1
	viewsfile2.write(line)
	continue
    #In newupload, the first else is right after the lesson lines
    if on_newupload and elseword.search(line):
        viewsfile2.write("        elif lnum == 'lesson" + lesson_num +"':\n")
        viewsfile2.write("            lessontype = 'lesson" + lesson_num + "'\n")
        viewsfile2.write("            lesson_obj = lesson" + lesson_num + ".models.Lesson" + lesson_num + "()\n")
	viewsfile2.write(line)
	on_newupload = 0
	continue
    viewsfile2.write(line)
viewsfile.close()
viewsfile2.close()

print "Adding wiki board for lesson \n"
boardname = "Lesson " + lesson_num
try:
    dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, db=settings.DATABASE_NAME)
    cra = dba.cursor()
except:
    print "Authentication to MySQL db failed!"
    sys.exit(-1)

try:
    cra.execute("INSERT INTO wiki_board (name) VALUES ('" + boardname + "')")
except:
    print "Could not insert row for wiki into db!"
    sys.exit(-1)

cra.close()
dba.close()

os.rename('mytemplates/html/skeleton.html','mytemplates/html/skeleton.html.prelesson')
os.rename('mytemplates/html/skeleton2.html','mytemplates/html/skeleton.html')
os.rename('mytemplates/html/fileuploadform.html','mytemplates/html/fileuploadform.html.prelesson')
os.rename('mytemplates/html/fileuploadform2.html','mytemplates/html/fileuploadform.html')
os.rename('settings.py','settings.py.prelesson')
os.rename('settings2.py','settings.py')
os.rename('lessonaux.py','lessonaux.py.prelesson')
os.rename('lessonaux2.py','lessonaux.py')
os.rename('lessons/urls.py','lessons/urls.py.prelesson')
os.rename('lessons/urls2.py','lessons/urls.py')
os.rename('pdbinfo/views.py','pdbinfo/views.py.prelesson')
os.rename('pdbinfo/views2.py','pdbinfo/views.py')

os.system('python manage.py syncdb')

print "Please restart the server to complete the process\n"
sys.exit()
