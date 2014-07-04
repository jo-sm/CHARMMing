#!/usr/bin/python

#makes lessons in python
import sys
import os
import re
import MySQLdb
import settings
import django
from django.template import Context, loader
from django.template.loader import get_template
from output import tidyTxt

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

#Check and see if the developer supplied any files the lesson user must download to continue with the lesson
try:
    lesson_filenames = sys.argv[3]
except:
    lesson_filenames = ""
lesson_filenameslist = []
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

maker_dict = {}
maker_dict['lesson_name'] = lesson_name 
maker_dict['lesson_num'] = lesson_num 
maker_dict['left'] = "{%"
maker_dict['right'] = "%}"
maker_dict['left2'] = "{{"
maker_dict['right2'] = "}}"
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/modelfile.txt')
modelsfile_text = t.render(Context(maker_dict))

modelfile.write(modelsfile_text)
modelfile.close()

print "Creating views.py ...\n"
viewsfile = open(lesson_name + "/views.py",'w+')
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/viewsfile.txt')
viewsfile_text = t.render(Context(maker_dict))

viewsfile.write(viewsfile_text)
viewsfile.close()

print "Creating mytemplates/html/" + lesson_name + ".html ...\n"
htmlfile = open('mytemplates/html/' + lesson_name + '.html','w+')

maker_dict['lesson_filenames'] = lesson_filenames
maker_dict['lesson_filenameslist'] = lesson_filenameslist
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/htmlfile.txt')
htmlfile_text = tidyTxt(t.render(Context(maker_dict)))
#htmlfile_text = t.render(Context(maker_dict))

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
maker_dict['settingsfile'] = []

for line in settingsfile:
    dict = {}
    dict['line'] = line
    if line.strip() == '':
        dict['line'] = "blank_line\n"
    dict['case1'] = ''
    dict['case2'] = ''
    if installedapps.search(line):
        instapp_init = 1
        instapp_firstline = 1
        dict['case1'] = 'y'
        maker_dict['settingsfile'].append(dict) 
        continue
    if instapp_init and endparanth.search(line):
        dict['case2'] = 'y'
        maker_dict['settingsfile'].append(dict) 
	continue
    maker_dict['settingsfile'].append(dict) 
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/settingsfile2.txt')

settingsfile2.write(tidyTxt(t.render(Context(maker_dict))))
#settingsfile2.write(t.render(Context(maker_dict)))
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
maker_dict['htmlfile'] = []
for line in htmlfile:
    dict = {}
    dict['line'] = line
    if line.strip() == '':
        dict['line'] = "blank_line\n"
    dict['case1'] = ''
    dict['case2'] = ''
    dict['case3'] = ''
    #the user has a lessons section and hasn't removed it
    if span.search(line) and lessons.search(line):
        lesson_exists = 1
        dict['case1'] = 'y'
        maker_dict['htmlfile'].append(dict) 
        continue
    #If the user does not have a lessons section, make one right before the accounts section
    if span.search(line) and account.search(line) and lesson_exists < 1:
        dict['case2'] = 'y'
        maker_dict['htmlfile'].append(dict) 
	continue
    #if they do have a lessons section then write the url right before the </div> in the lessons section
    if enddiv.search(line) and lesson_exists == 1:
        dict['case3'] = 'y'
	lesson_exists = 2
        maker_dict['htmlfile'].append(dict) 
	continue
    maker_dict['htmlfile'].append(dict) 
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/htmlfile2.txt')
htmlfile2.write(tidyTxt(t.render(Context(maker_dict))))
htmlfile2.close()   
      
print "Adding the lesson into the Submit Structure page... mytemplates/html/fileuploadform.html \n"
fileuploadfile = open('mytemplates/html/fileuploadform.html','r')
fileuploadfile2 = open('mytemplates/html/fileuploadform2.html','w+')
whatlesson = re.compile("What Lesson is this PDB")
endselect = re.compile("/select")
submitstructure = re.compile("Submit Structure")
onlesson = 0
fileuploadfile_text = ""
maker_dict['upload_lis'] = []

for line in fileuploadfile:
    dict = {}
    dict['line'] = line
    if line.strip() == '':
        dict['line'] = "blank_line\n"
    dict['case1'] = ''
    dict['case2'] = ''
    dict['case3'] = ''

    if whatlesson.search(line):
        onlesson = 1
        dict['case1'] = 'y'
        maker_dict['upload_lis'].append(dict)
	continue

    #Puts the lesson at the end of the list on fileuploadform
    if endselect.search(line) and onlesson:
        dict['case2'] = 'y'
        maker_dict['upload_lis'].append(dict)
 	continue
    #This checks if the user removed the fileuploadform lesson box

    if submitstructure and onlesson == 0:
        fileuploadform2_text = '''
        <select name="lesson" id="lesson"> 
         <option value="nolesson" selected> No Lesson </option>
         <option value="''' + lesson_name + '''"> Lesson ''' + lesson_num + ''' </option>
        </select><br><br>
	'''
        dict['case3'] = 'y'
        maker_dict['upload_lis'].append(dict)
 	continue
#    print line
    maker_dict['upload_lis'].append(dict)

#for dict in maker_dict['upload_lis']:
#    print dict['line']   
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/fileuploadfile2.txt')
fileuploadfile2.write(tidyTxt(t.render(Context(maker_dict))))
fileuploadfile.close()
fileuploadfile2.close()
       
print "Adding lesson into lessonwrapper located in ./lessonaux.py doLessonAct\n"
lessonaux = open('lessonaux.py','r')
lessonaux2 = open('lessonaux2.py','w+')

importpdbinfo = re.compile('import pdbinfo')
onfileupload = re.compile("function == 'onFileUpload'")
on_dolessonact = 0
maker_dict['lesson_aux'] = []
for line in lessonaux:
    dict = {}
    dict['line'] = line
    if line.strip() == '':
        dict['line'] = "blank_line\n"
    dict['case1'] = 0
    dict['case2'] = 0    
    #import the lessons
    if importpdbinfo.search(line):
        dict['case1'] = 1
    #add the logic for the lesson code
    if onfileupload.search(line):
        dict['case2'] = 1
    maker_dict['lesson_aux'].append(dict)
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/lessonaux2.txt')
lessonaux2.write(tidyTxt(t.render(Context(maker_dict))))
lessonaux.close()
lessonaux2.close()

print "Adding lesson into the lessons/urls.py file\n"
urls = open('lessons/urls.py','r')
urls2 = open('lessons/urls2.py','w+')
admin = re.compile('admin/')
maker_dict['urls_lis'] = []
for line in urls:
    dict = {}
    dict['line'] = line
    if line.strip() == '':
        dict['line'] = "blank_line\n"
    dict['case1'] = 0
    if admin.search(line):
        dict['case1'] = 1
    maker_dict['urls_lis'].append(dict)    
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/urls2.txt')
urls2.write(tidyTxt(t.render(Context(maker_dict))))
    
urls.close()
urls2.close()

print "Creating custom urls.py for lesson. " + lesson_name + "/urls.py\n"

urls = open(lesson_name + "/urls.py",'w+')
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/urls.txt')

urls.write(t.render(Context(maker_dict)))
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
maker_dict['viewsfile_lis'] = []
for line in viewsfile:
    dict = {}
    dict['line'] = line
    if line.strip() == '':
        dict['line'] = "blank_line\n"
    dict['case1'] = ''
    dict['case2'] = ''
    dict['case3'] = ''
    dict['case4'] = ''
    dict['case5'] = ''
    if importos.match(line):
        dict['case1'] = 'y'
        maker_dict['viewsfile_lis'].append(dict) 
	continue
    if handlecrdpsf.search(line):
        on_handlecrdpsf = 1
        dict['case2'] = 'y'
        maker_dict['viewsfile_lis'].append(dict) 
	continue
    #adds it into the handlecrdpsf area
    if on_handlecrdpsf == 1 and editpdbinfo.search(line):
        dict['case3'] = 'y'
	on_handlecrdpsf = 0
        maker_dict['viewsfile_lis'].append(dict) 
	continue
    if newupload.search(line):
        on_newupload = 1
        dict['case4'] = 'y'
        maker_dict['viewsfile_lis'].append(dict) 
	continue
    #In newupload, the first else is right after the lesson lines
    if on_newupload and elseword.search(line):
        dict['case5'] = 'y'
	on_newupload = 0
        maker_dict['viewsfile_lis'].append(dict) 
	continue
    maker_dict['viewsfile_lis'].append(dict) 
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/viewsfile2.txt')

viewsfile2.write(tidyTxt(t.render(Context(maker_dict))))
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
