import structure
from django.db import models
from django.contrib.auth.models import User
from lessons.models import LessonProblem
from dynamics.models import mdTask
import os, math
import structure, lessonaux, charmming_config
import re

class Lesson5(models.Model):
    user = models.ForeignKey(User)
    nSteps = models.PositiveIntegerField(default=5)
    curStep = models.DecimalField(default=0,decimal_places=1,max_digits=3)

    def onFileUpload(self):
        try:
            LessonProblem.objects.get(lesson_type='lesson5',lesson_id=self.id).delete()
        except:
            pass
        file = structure.models.Structure.objects.filter(selected='y',lesson_id=self.id)[0]
        try:
            filename3 = '%s/mytemplates/lessons/lesson5/1prb.pdb' % charmming_config.charmming_root
            os.stat(filename3)
            filename4 = file.location + '/1prb.pdb'
            os.stat(filename4)
        except:
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=1,severity=9,description='The structure you submitted did not upload properly. Check to make sure the structure files are valid.')
            lessonprob.save()
            return False
        if not lessonaux.diffPDBs(file,filename3,filename4):
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=1,severity=9,description='The structure you uploaded was not the correct structure.')
            lessonprob.save()
            return False
        self.curStep = '1'
        self.save()
        return True

    def onEditPDBInfo(self,postdata):
        return True

    def onMinimizeSubmit(self,mp,filename):
        try:
            LessonProblem.objects.filter(lesson_type='lesson5',lesson_id=self.id)[0].delete()
        except:
            pass

        if mp.sdsteps != 10:
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=3,severity=2,description='SD steps were not set to 10.')
            lessonprob.save()
            return False
        if mp.abnrsteps != 100:
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=3,severity=2,description='ABNR steps were not set to 100.')
            lessonprob.save()
            return False
        if float(mp.tolg) != .05:
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=3,severity=2,description='TOLG was not set not 0.05.')
            lessonprob.save()
            return False
        if mp.usepbc == 't':
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=3,severity=2,description='Somehow your model is using periodic boundary conditions. This should not be possible and is a mistake.')
            lessonprob.save()
            return False

        self.curStep = '2.5'
        self.save()
        return True

    def onMinimizeDone(self,mp):

        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson5',lesson_id=self.id)[0]
        except:
            lessonprob = None

        logfp = open('/tmp/l5mindone.txt','w')
        logfp.write('called\n')
        logfp.flush() 

        #mp = minimizeTask.objects.filter(pdb=file,selected='y')[0]
        if lessonprob:
            self.curStep = '2'
            self.save()
            logfp.write('Got lessonprob\n')
            logfp.close()
            return False

        if mp.status=='F':
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=3,severity=9,description='The Job did not complete correctly.')
            lessonprob.save()
            self.curStep = '2'
            self.save()
            logfp.write('min failed\n')
            logfp.close()
            return False
        else:
            self.curStep = '3'
            self.save()
            logfp.write('all good\n')
            logfp.close()
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

    def onLDSubmit(self,ldp):
        logfp = open('/tmp/l5-sub.txt','w')
        logfp.write('lesson 5 submits.\n')
        logfp.close()

        try:
            LessonProblem.objects.filter(lesson_type='lesson5',lesson_id=self.id)[0].delete()
        except:
            pass
        if ldp.parent.action != 'minimization':
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=4,severity=2,description='Please run :angevin dynamics on the PDB from the minimization.')
            lessonprob.save()
            return False
        if int(ldp.nstep) != 1000:
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=4,severity=2,description='Please set the number of steps to 1000 to continue.')
            lessonprob.save()
            return False
        if float(ldp.fbeta) != 1.0:
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=4,severity=2,description='Please set the collision frequency (FBETA).')
            lessonprob.save()
            return False

        self.curStep = '3.5'
        self.save()
        return True

    def onLDDone(self,ldp):
        logfp = open('/tmp/l5-done.txt','w')
        logfp.write('lesson 5 done.\n')
        logfp.close()

        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson5',lesson_id=self.id)[0]
        except:
            lessonprob = None
        if lessonprob:
            return False
        if  ldp.status == 'F':
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=4,severity=9,description='The job did not complete correctly.')
            lessonprob.save()
            return False
        else:
            self.curStep = '4'
            self.save()
        return True

    def onMDSubmit(self,postdata,filename):
        return True

    def onMDDone(self,file):
        return True

    def onSGLDSubmit(self,postdata):
        return True

    def onSGLDDone(self,file):
        return True

    def onBuildStructureSubmit(self,postdata):
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson5',lesson_id=self.id)[0].delete()
        except:
            lessonprob = None

        if postdata['buildtype'] != 'go':
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=2,severity=9,description='You did not build a Go model. Please be sure you select a Go model from the CG working structure page.')
            lessonprob.save()
            return False

        if float(postdata['gm_nscale']) != 0.91:
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=2,severity=9,description='Wrong nScale set. Please set nScale to 0.91')
            lessonprob.save()
            return False

        if postdata['gm_contact_type'] != 'mj':
            lessonprob = LessonProblem(lesson_type='lesson5',lesson_id=self.id,errorstep=2,severity=9,description='You did not use the MJ contact set. Please use this contact set.')
            lessonprob.save()
            return False

        self.curStep = '2'
        self.save()
        return True

    def onBuildStructureDone(self,file):
        self.curStep = '2'
        self.save()
        return True

    def onRedoxSubmit(self,postdata):
        return True

    def onRedoxDone(self,file):
        return True

    def onNATQSubmit(self,postdata):
        # no checking
        try:
            LessonProblem.objects.filter(lesson_type='lesson5',lesson_id=self.id)[0].delete()
        except:
            pass
        self.curStep = '5'
        self.save()
        return True

    def onRMSDSubmit(self,postdata):
        return True

    #generates html for the lesson status page
    def generateStatusHtml(self,file):
        step_status_list = []
        step_status_list.append("<tr class='status'><td class='status'>1. File Uploaded: ")
        step_status_list.append("<tr class='status'><td class='status'>2. Working Structure Built: ")
        step_status_list.append("<tr class='status'><td class='status'>3. Light minimization: ")
        step_status_list.append("<tr class='status'><td class='status'>4. Langevin dynamics: ")
        step_status_list.append("<tr class='status'><td class='status'>5. Analysis: ")
        #This will store all the status and the steps, clearing the template of logic
        #And only displaying the status
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson5',lesson_id=self.id)[0]
        except:
            lessonprob = None
        for i in range(self.nSteps):
            if lessonprob and lessonprob.errorstep == math.floor(self.curStep+1) and math.floor(self.curStep) == i:
                step_status_list[i] += ("<font color='red'>Failed</font></td></tr>")
                continue
            elif (float(self.curStep)-0.5) == i and float(self.curStep) % 1.0 == 0.5:
                step_status_list[i] += ("<font color='blue'>Running</font></td></tr>")
                continue
            elif i < float(self.curStep) or len(step_status_list) == float(self.curStep):
                step_status_list[i] += ("<font color='green'>Done</font></td></tr>")
                continue
            elif i + 1 > float(self.curStep):
                step_status_list[i] += ("<font color='grey'>N/A</font></td></tr>")
                continue

        return step_status_list

    def getHtmlStepList(self):
        #2 is running
        #1 is successfully done
        #0 is not started
        #-1 is error
        htmlcode_list = []
        for step in range(self.nSteps):
            htmlcode_list.append(0)
        htmlcode_list.append(self.curStep)
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
        if float(self.curStep) > 4:
            if float(self.curStep) == 4.5:
                htmlcode_list[4] = 2
            else:
                htmlcode_list[4] = 1
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson5',lesson_id=self.id)[0]
            htmlcode_list[lessonprob.errorstep-1] = -1
        except:
            lessonprob = None
        return htmlcode_list

