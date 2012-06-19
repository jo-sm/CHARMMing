
# lesson 1, upload 1YJP, minimize in vacuum, solvate/neutralize, minimize again, and
# run dynamics
from django.db import models
from django.contrib.auth.models import User
from lessons.models import LessonProblem
from minimization.models import minimizeTask
import os, re
import structure, lessonaux, charmming_config


class Lesson4(models.Model):
    # data for lessons (should not be overridden by subclasses)
    # specifying both user and Structure is redundant (since the Structure references the user),
    # but makes it easier to look up all lessons being done by a particular user.
    user = models.ForeignKey(User)
    nSteps = models.PositiveIntegerField(default=4)
    curStep = models.DecimalField(default=0,decimal_places=1,max_digits=3)

    
    def onFileUpload(self):
        try:
            LessonProblem.objects.filter(lesson_type='lesson4',lesson_id=self.id)[0].delete()
        except:
            pass
        file = structure.models.Structure.objects.get(selected='y',owner=self.user,lesson_id=self.id)
        #all_segids = file.segids.split() + file.good_het.split() + file.nongood_het.split()
        try:
            filename1 = file.location + '/butane.crd'
            os.stat(filename1)
        #    filename2 = file.location + 'rtf-' + file.stripDotPDB(file.filename) + '.crd.rtf'
        #    os.stat(filename2)
        except:
            lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=1,severity=9,description='The coordinate file was not uploaded properly.')
            lessonprob.save()
            return False
        try:
            #filename1 = '%s/mytemplates/lessons/lesson4/butane.psf' % charmming_config.charmming_root
            #os.stat(filename1)
            #filename2 = file.location + 'psf-' + file.stripDotPDB(file.filename) + '.psf'
            #os.stat(filename2)
            #if not lessonaux.diffPDBs(file,filename1,filename2):
            #    lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=1,severity=9,description='The PDB you uploaded was not the correct PDB.')
            #    lessonprob.save()
            #    return False
            filename1 = '%s/mytemplates/lessons/lesson4/butane.crd' % charmming_config.charmming_root
            os.stat(filename1)
            filename2 = file.location + '/butane.crd'
            os.stat(filename2)
            if not lessonaux.diffPDBs(file,filename1,filename2):
                lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=1,severity=9,description='The CRD you uploaded was not the correct PDB.')
                lessonprob.save()
                return False
        except:
            lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=1,severity=9,description='The CRD/PSF you submitted did not upload properly.')
            lessonprob.save()
            return False
        self.curStep = '1.0'
        self.save()
        return True

    def onEditPDBInfo(self,postdata):
        return True

    def onMinimizeSubmit(self,mp,filename):
        try:
            LessonProblem.objects.filter(lesson_type='lesson4',lesson_id=self.id)[0].delete()
        except:
            pass
        #file = structure.models.Structure.objects.filter(selected='y',owner=self.user,lesson_id=self.id)[0]
        #mp = minimizeTask.objects.filter(pdb=file,selected='y')[0]
        if mp.sdsteps != 100:
            lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=4,severity=2,description='SD steps were not set to 100.')
            lessonprob.save()
            return False
        if mp.abnrsteps != 1000:
            lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=4,severity=2,description='ABNR steps were not set to 1000.')
            lessonprob.save()
            return False
        if float(mp.tolg) != .05:
            lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=4,severity=2,description='TOLG was not set not 0.05.')
            lessonprob.save()
            return False
        # 3.5 Means it is running
        self.curStep = '3.5'
        self.save()
        return True

    def onMinimizeDone(self,mp):
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson4',lesson_id=self.id)[0]
        except:
            lessonprob = None
        #mp = minimizeTask.objects.filter(pdb=file,selected='y')[0]
        fail = re.compile('Failed')
        if lessonprob:
            return False
        if mp.status=='F':
        #if fail.search(mp.statusHTML):
            lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=3,severity=9,description='The Job did not complete correctly.')
            lessonprob.save()
            return False
        else:
            self.curStep = '4'
            self.save()
        return True

    def onSolvationSubmit(self,postdata):
        return True

    def onSolvationDone(self,file):
        return True

    def onNMASubmit(self,postdata):
        return True

    def onNMADone(self,file):
        return True

    def onMDSubmit(self,postdata,filename):
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

    def onEnergyDone(self,et):
        l4log=open("/tmp/l4.log",'aOA')
        
        try:
            lessonprob = LessonProblem.objects.get(lesson_type='lesson4',lesson_id=self.id)
        except:
            lessonprob = None
        #mp = minimizeTask.objects.filter(pdb=file,selected='y')[0]
        #fail = re.compile('Failed')
        
        if float(et.finale)==0.0: return True
        if lessonprob:
            l4log.write("there is a problem: %s\n" % (str(lessonprob.id) + " " +lessonprob.description))
            if float(self.curStep)==1.5: self.curStep = '1'
            if float(self.curStep)==2.5: self.curStep = '2'
            self.save()
            return False
        #l4log.write("curstep is:%s\n" % str(float(self.curStep))) 
        if float(self.curStep) == 1.5:
            #l4log.write("finale: %s\n" % (et.finale))
            # finale should be 6.84842 but we'll give a slight margin for error
            if float(et.finale) < float(6.8) or float(et.finale) > float(6.9):
                #t4log.write("saving error finale: %s\n" % str(float(et.finale)))
                lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=2,severity=9,description="Your energy was wrong blah(should be ~ 6.864, you got %s). Please check that the RTF file you uploaded was correct. If you can't figure out the problem, the correct RTF is <a href=/charmming/lessons/lesson4/butane.rtf>here</a>." % str(float(et.finale)))
                lessonprob.save()
                return False
            self.curStep = '2'
            self.save()
            return True
        elif float(self.curStep) == 2.5:
            # finale should be -49122.53689 but we'll give a slight margin for error
            #l4log.write("qmfinale: %s\ntask is: %s" % (str(float(et.finale)),et.id))
            if float(et.finale) < float(-49037.5) or float(et.finale) > float(-49034.5):
            #if float(et.finale) < -49124.5 or float(et.finale) > -49121.5:
                lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=3,severity=9,description='Your energy was wrong (should be ~ -49035.98). Please check that you set up the QM/MM calculation correctly.')
                lessonprob.save()
                return False
            self.curStep = '3'
            self.save()
            return True
        else:
            return False

    def onRMSDSubmit(self,file):
        return True

    #generates html for the lesson status page
    def onEnergySubmit(self,et):
        #l4submitlog=("/tmp/l4submit.log",'w')
        try:
            LessonProblem.objects.get(lesson_type='lesson4',lesson_id=self.id).delete()
        except:
            pass
        if float(self.curStep) == 1.0:
            # nothing really to validate here...
            self.curStep='1.5'
            self.save()
            return True
        elif float(self.curStep) == 2.0:
            #if not postdata.has_key('useqmmm'):
            if et.useqmmm!='y':
                lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=3,severity=9,description='You did not use QM/MM in your energy calculation.')
                lessonprob.save()
                return False
            if et.qmmmsel !='type c1 .or. type h11 .or. type h12 .or. type h13 .or. type c2 .or. type h21 .or. type h22':
            #if et.qmmmsel != 'atom buta 1 c1 .or. atom buta 1 c2 .or. atom buta 1 h11 .or. atom buta 1 h12 .or. atom buta 1 h13 .or. atom buta 1 h21 .or. atom buta 1 h22':
                lessonprob = LessonProblem(lesson_type='lesson4',lesson_id=self.id,errorstep=3,severity=9,description='You did not make the correct QM/MM selection. Please re-read the instructions.')
                return False
            self.curStep='2.5'
            self.save()
            return True
        else:
            return True

    def onRMSDSubmit(self,file):
        return True

    #generates html for the lesson status page
    def generateStatusHtml(self,file):
        step_status_list = []
        step_status_list.append("<tr class='status'><td class='status'>1. File Uploaded: ") 
        step_status_list.append("<tr class='status'><td class='status'>2. Classical Energy: ") 
        step_status_list.append("<tr class='status'><td class='status'>3. QM/MM Energy: ") 
        step_status_list.append("<tr class='status'><td class='status'>4. QM/MM Minimization: ") 
        #This will store all the status and the steps, clearing the template of logic
        #And only displaying the status
        try:
            lessonprob = LessonProblem.objects.get(lesson_type='lesson4',lesson_id=self.id)
        except:
            lessonprob = None
        for i in range(self.nSteps):
            if lessonprob and lessonprob.errorstep == (self.curStep+1) and self.curStep == i:
                step_status_list[i] += ("<font color='red'>Failed</font></td></tr>")
                continue
            elif (float(self.curStep)-0.5) == i and float(self.curStep) % 1 == 0.5:
                step_status_list[i] += ("<font color='blue'>Running</font></td></tr>")
                continue
            elif i < float(self.curStep):
                step_status_list[i] += ("<font color='green'>Done</font></td></tr>")
                continue
            elif i >= float(self.curStep):
                step_status_list[i] += ("<font color='grey'>N/A</font></td></tr>")
                continue
        return step_status_list 

    #Returns a list where each index corresponds to lesson progress
    #on the display lesson page
    def getHtmlStepList(self):
        #2 is running
        #1 is successfully done
        #0 is not started
        #-1 is error
        htmlcode_list = []
        for step in range(self.nSteps):
            htmlcode_list.append(0)
        if float(self.curStep) > 0:
            htmlcode_list[0] = 1
        if float(self.curStep) > 1:
            if float(self.curStep) == 1.5:
                htmlcode_list[1] = 2
            else:
                htmlcode_list[1] = 1
        if float(self.curStep) > 2:
            if float(self.curStep) == 2.5:
                htmlcode_list[2] = 2
            else:
                htmlcode_list[2] = 1
        if float(self.curStep) > 3:
            if float(self.curStep) == 3.5:
                htmlcode_list[3] = 2
            else:
                htmlcode_list[3] = 1
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson4',lesson_id=self.id)[0]
            htmlcode_list[lessonprob.errorstep-1] = -1
        except:
            lessonprob = None
        return htmlcode_list

