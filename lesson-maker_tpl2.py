#!/usr/bin/python

#makes lessons in python
import sys, getopt
import os
import re
import MySQLdb
import settings
import django
from django.template import Context, loader
from django.template.loader import get_template
from output import tidyTxt
from lesson_config import *

# A brief documentation for testing only
# The final documentation will be in a text file called 'DOCUMENTATION'
# Then a command '' will desplay all documentation  
def usage():
    print "                       LESSON-MAKER DOCUMENTATION" 
    print "The program 'lesson-maker_tpl2.py' is used to create a new lesson by users.\nWith the lesson users can uploade files, implement minimization, solvation and so on through 'charmming'.\n"
    print "Below is a brief description of command line option usage."
    print "python lesson-maker_tpl2.py -n(-name) <<new_lesson_name>> -d(--directory) <<location_of_files_the_user_must_download>> -f(--files) [name_of_files_user_must_download_in_brackets_seperated_by_a_comma]" 
    print "-r(--description) <<a_text_file_of_lesson_description>> -c(--config) <<lesson_configuration_location>> -h(--help) <<documentation_of_help>>"
    print "NOTE: THIS IS VERY IMPORTANT. ANY OPTION VALUE MUST BE WITHOUT ANY SPACE!!! OPTION-VALUR MUST BE PAIRED EXCEPT FOR OPTION -h/--help"
    print "Option '-h/--help' will allow users to find how to run 'lesson-maker_tpl2.py'. A detailed example is whown through this option."
    print "Option '-n/--name' specifies the name of the lesson such as 'Lesson1', 'Lesson2' and so on."
    print "Option '-r/--description is a 'txt' file specifying what the lesson will do. If the file is in a directory other than the working directory,"
    print "Option '-f/--files' allow users to put file names needed for the lesson. NOTE the file names must be put in brackets and seperated by comma. "
    print "NOTE again no spaces are allowed in brackets!!! Otherwise the rest files after first space will be lost."
    print "the right path should be preceding the file name. User could be able to modify it later by editting 'lesson_config.py'."
    print "'lesson_config.py' is placed in the working directory or the directory specified by option -c(--config)."
    print "The option of '-f/--files' and '-d/--directory' are necessary only if users need to upload pdb, crd, psf or any other files in order to create the lesson.\n"
    print "Following is a detailed example demonstrating how to create a new lesson with lesson-maker_tpl2.py."

if len(sys.argv) == 1 or not sys.argv[1].startswith('-'):
    usage()
    sys.exit(1)

try:
    opts, args = getopt.getopt(sys.argv[1:], 'n:hd:f:r:c:', ['name=', 'help', 'description=', 'files=', 'directory=', 'config='])

except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)

# as default 
new_lesson_type = 'Put the File Type Here'
new_lesson_desc = "Put the Lesson Description Here"
lesson_filenameslist = []
lesson_file_location = '.'
lesson_config_location = '.'
# maker_dict will pass the information needed by templates
maker_dict = {}

for o, a in opts:
    if a.startswith('-') or not o.startswith('-'):
        print "Command line option error. Options must start with '-' while values must not start with '-'."
        print "Also options and values must be paired except for -h/--help. Use -h/--help for detailed information."
        sys.exit() 
         
    if o in ("-h", "--help"):
        usage()
        sys.exit(2)
    elif o in ("-n", "--name") and a[6:]:
        try:
            new_lesson_num = a[6:]
            if not new_lesson_num.isdigit():
                print "Lesson name must be a format of 'lesson' plus a digit number. Use -h/--help for detailed information."
                sys.exit() 
            if new_lesson_num in lesson_num_lis:
                print "This lesson is already created. The lessons you have created are in the list: "
                print lesson_num_lis
                print "Please create a new lesson.\n"
                sys.exit(3)
            lesson_name = "lesson" + new_lesson_num
            # as default
            new_lesson_txt = 'Lesson '+ new_lesson_num
            # if wish a text format of the lesson modify the corresponding term in 'lesson_config.py'
        except:
            print "Lesson name must be provided. Use 'python lesson-maker_tpl2.py -h/--help' for detailed information\n"
            sys.exit(3)

    #Check and see if users supplied a directory where the files must be downloaded 
    elif o in ("-d", "--directory") and a:
        try: 
            if os.path.isdir(a):
                lesson_file_location = a
            else:
                print "The directory is not specified correctly. Please check and try again."
                sys.exit() 
        except:
            print "No file location was specified, assuming all files are in current working directory\n"

    #Check and see if users supplied any files which must be downloaded to continue with the lesson
    #If there are files that accompany the lesson, a location of such files is created to store the files
    elif o in ("-f", "--files") and a:
        try: 
            print "Copying the files that accompany lesson into mytemplates/lessons/" + lesson_name
            for type in ('.crd','.psf','.CRD','.PSF'):
                if a.find(type) != -1:
                    new_lesson_type = 'CRD/PSF'
                    break
            lesson_filenameslist = a.strip('[]').split(',')
            for filename in lesson_filenameslist:
                if not os.path.isfile(filename):
                    print "%s does not exist or the path is not specified correctly. Please check and try again." % filename
                    sys.exit() 
            try:
                os.stat('mytemplates/lessons/' + lesson_name)
            except:
                os.mkdir('mytemplates/lessons/' + lesson_name)
            for filename in lesson_filenameslist:
                if os.path.isfile(filename):
                    os.system("cp " + lesson_file_location + "/" + filename + " mytemplates/lessons/" + lesson_name + "/")
        except:
            print "No files were given that should be used in this lesson. It is assumed that in this lesson, the user will not have to download any other files to upload.\n"
            print "If you would like to add files later, store them in mytemplates/lessons/" + lesson_name + "/\n"
    elif o in ("-r", "--description") and a:
        try:             
            if os.path.isfile(a):
                f = open(a, 'r')
                for l in f:     
                    # make only one line
                    new_lesson_desc += l.strip()+" "
            else:
                 print "%s does not exist or the path is not specified correctly. Please check and try again." % a
                 sys.exit()
        except:
            new_lesson_desc = 'Put the lesson description here.'  
    elif o in ("-c", "--config") and a:
        try: 
            lesson_config_location = a
        except:
            pass
    else:
        print "'%s' unhandled option!!!" % o
        usage()
        sys.exit(3)

os.system("python manage.py startapp " + lesson_name)

maker_dict['lesson_name'] = lesson_name 
maker_dict['lesson_num'] = new_lesson_num 
maker_dict['left'] = "{%"
maker_dict['right'] = "%}"
maker_dict['left2'] = "{{"
maker_dict['right2'] = "}}"
# a list containing lesson information of lesson number, description and its text-format such as 'Lesson One', 'Lesson Two' 
lesson_info_lis = []

# add lesson information into lesson_info_lis and pass to templates
for num in lesson_num_lis:
    lesson_info_lis.append({'lesson_num':num, 'lesson_txt':lesson_txt[num], 'lesson_desc':lesson_desc[num], 'file_type':file_type[num]})
lesson_info_lis.append({'lesson_num':new_lesson_num, 'lesson_txt':new_lesson_txt, 'lesson_desc':new_lesson_desc, 'file_type':file_type[num]})
maker_dict['lesson_info_lis'] = lesson_info_lis

# those variables are passed to templates too
maker_dict['new_lesson_txt'] = new_lesson_txt
maker_dict['new_lesson_desc'] = new_lesson_desc
maker_dict['new_lesson_type'] = new_lesson_type
#maker_dict['lesson_num_lis'] = lesson_num_lis.append(new_lesson_num)

# create 'lesson_config.py' so other 'py' programs can import it 
lesson_config = open(lesson_config_location+'/lesson_config.py', 'w')
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/lesson_config.txt')
lesson_config.write(tidyTxt(t.render(Context(maker_dict))))
lesson_config.close()

print "Creating models.py ...\n"
modelfile = open(lesson_name + "/models.py",'w+')
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

maker_dict['lesson_filenameslist'] = lesson_filenameslist
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/htmlfile.txt')
htmlfile_text = tidyTxt(t.render(Context(maker_dict)))
htmlfile.write(htmlfile_text)
htmlfile.close()

print "Installing lesson into mytemplates/html/skeleton.html\n"
htmlfile2 = open('mytemplates/html/skeleton2.html','w+')
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/htmlfile2_tpl.txt')
htmlfile2.write(tidyTxt(t.render(Context(maker_dict))))
htmlfile2.close()   
      
print "Adding the lesson into the Submit Structure page... mytemplates/html/fileuploadform.html \n"
fileuploadfile2 = open('/var/www/html/charmming/mytemplates/html/fileuploadform2.html','w+')
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/fileuploadfile2_tpl.txt')
fileuploadfile2.write(tidyTxt(t.render(Context(maker_dict))))
fileuploadfile2.close()
       
print "Creating custom urls.py for lesson. " + lesson_name + "/urls.py\n"
urls = open(lesson_name + "/urls.py",'w+')
t = get_template('/var/www/html/charmming/mytemplates/lesson_maker/urls.txt')
urls.write(t.render(Context(maker_dict)))
urls.close()

print "Adding wiki board for lesson \n"
boardname = "Lesson " + new_lesson_num
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

os.system('python manage.py syncdb')

print "Please restart the server to complete the process\n"
sys.exit()
