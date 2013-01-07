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

import lesson1
import lesson2
import lesson3
import lesson4
#import lesson6
#import lesson96
# by importing lesson_config to import all lessons created
from lesson_config import *
import sys, traceback

def doLessonAct(file,function,task=None,filename=None,finale=None):
    lessontype = file.lesson_type
    lessonid = file.lesson_id

    estr = lessontype+'.models.'+lessontype.capitalize()+'()'
    lesson_num_obj = eval(estr)
    lesson_num_class=lesson_num_obj.__class__
    lesson_obj = lesson_num_class.objects.get(id=lessonid)

    if function == 'onFileUpload':
        lesson_obj.onFileUpload()
        return True
    elif function == 'onMinimizeSubmit':
        lesson_obj.onMinimizeSubmit(task,filename)
        return True
    elif function == 'onMinimizeDone':
        lesson_obj.onMinimizeDone(task)
        return True
    elif function == 'onSolvationSubmit':
        lesson_obj.onSolvationSubmit(task)
        return True
    elif function == 'onSolvationDone':
        lesson_obj.onSolvationDone(task)
        return True
    elif function == 'onMDSubmit':
        lesson_obj.onMDSubmit(task,filename)
        return True
    elif function == 'onMDDone':
        lesson_obj.onMDDone(task)
        return True
    elif function == 'onLDSubmit':
        lesson_obj.onLDSubmit(task)
        return True
    elif function == 'onLDDone':
        lesson_obj.onLDDone(file)
        return True
    elif function == 'onSGLDSubmit':
        lesson_obj.onSGLDSubmit(task)
        return True
    elif function == 'onSGLDDone':
        lesson_obj.onSGLDDone(task)
        return True
    elif function == 'onEnergySubmit':
        lesson_obj.onEnergySubmit(task)
        return True
    elif function == 'onEnergyDone':
        lesson_obj.onEnergyDone(task)
        return True
    elif function == 'onRMSDSubmit':
        lesson_obj.onRMSDSubmit(task)
        return True
    elif function == 'onRedoxSubmit':
        lesson_obj.onRedoxSubmit(task)
        return True
    elif function == 'onRedoxDone':
        lesson_obj.onRedoxDone(task)
        return True
    elif function == 'getObject':
        return lesson_obj

def getPDBListFromPostdata(file,postdata):
    filename_list = file.getLimitedFileList('blank')
    pdblist = []
    for i in range(len(filename_list)):
        if postdata.has_key(filename_list[i]):
            pdblist.append(filename_list[i])
        elif postdata.has_key('min'):
            pdblist.append(postdata['min'])
            break
        elif postdata.has_key('solv_or_min'):
            pdblist.append(postdata['solv_or_min'])
            break
        else:
            continue
    return pdblist


def diffPDBs(file,filename1,filename2):
    file1list = []
    file2list = []

    file1 = open(filename1,'r')
    for line in file1:
        file1list.append(line.strip())
    file1.close()
    file2 = open(filename2,'r')
    for line in file2:
        file2list.append(line.strip())
    file2.close()
    if file1list != file2list:
        return False
    return True
