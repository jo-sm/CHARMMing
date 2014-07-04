
# lesson 6, upload 1CKU, calculate energy, calculate redox potential, change HIS protonation state,
# calculate redox potential, mutate residue, calculate redox potential, compare
from django.db import models
from django.contrib.auth.models import User
from lessons.models import LessonProblem
from apbs.models import redoxTask
from dynamics.models import mdTask
import os, math
import structure, lessonaux, charmming_config
import re

class Lesson6(models.Model):
    # data for lessons (should not be overridden by subclasses)
    # specifying both user and Structure is redundant (since the Structure references the user),
    # but makes it easier to look up all lessons being done by a particular user.
    user = models.ForeignKey(User)
    nSteps = models.PositiveIntegerField(default=5)
    curStep = models.DecimalField(default=0,decimal_places=1,max_digits=5)

    def LogStep(self,other='Done'):
        log=open("/tmp/lesson6_step%s.txt" % self.curStep,'w')
        log.write("lesson_id: %s\n" % (self.id))
        log.write("nSteps: %s\n" % (self.nSteps))
        log.write("curSteps: %s\n" % (self.curStep))
        log.write("%s\n" % other)
        log.close()

    def onFileUpload(self):
        self.LogStep()
        try:
            LessonProblem.objects.get(lesson_type='lesson6',lesson_id=self.id).delete()
        except:
            pass
        file = structure.models.Structure.objects.filter(selected='y',lesson_id=self.id)[0]
        try:
            filename3 = '%s/mytemplates/lessons/lesson6/1CKU.pdb' % charmming_config.charmming_root
            os.stat(filename3)
            filename4 = file.location+ '/1cku.pdb'
            os.stat(filename4)
        except:
            lessonprob = LessonProblem(lesson_type='lesson6',lesson_id=self.id,errorstep=1,severity=9,description='The structure you submitted did not upload properly. Check to make sure the structure files are valid.')
            lessonprob.save()
            return False
        if not lessonaux.diffPDBs(file,filename3,filename4):
            lessonprob = LessonProblem(lesson_type='lesson6',lesson_id=self.id,errorstep=1,severity=9,description='The structure you uploaded was not the correct structure.')
            lessonprob.save()
            return False
        self.curStep = '0.0'
        self.save()
        return True

    def onEditPDBInfo(self,postdata):
        return True

    def onBuildStructureSubmit(self,postdata):
        '''
            I combined this step with the energy calculations. - BSP
        '''
        #self.LogStep()
        #try:
        #    LessonProblem.objects.get(lesson_type='lesson6',lesson_id=self.id).delete()
        #except:
        #    pass
        #self.curStep = '1.5'
        #self.save()
        return True

    def onBuildStructureDone(self,postdata):
        '''
            I combined this step with the energy calculations. - BSP
        '''
        #self.LogStep()
        #self.curStep = '2'
        #self.save()
        return True

    def onEnergySubmit(self,postdata):
        if float(self.curStep) < 2.0:
            self.LogStep()
            try:
                LessonProblem.objects.get(lesson_type='lesson6',lesson_id=self.id).delete()
            except:
                pass
            self.curStep = '0.5'
            self.save()
        return True

    def onEnergyDone(self,postdata):
        if float(self.curStep) < 2.0:
            self.LogStep()
            try:
                LessonProblem.objects.get(lesson_type='lesson6',lesson_id=self.id).delete()
            except:
                pass
            self.curStep = '1'
            self.save()
        return True

    def onRedoxSubmit(self,postdata):
        self.LogStep()
        try:
            LessonProblem.objects.get(lesson_type='lesson6',lesson_id=self.id).delete()
        except:
            pass
        #file = structure.models.Structure.objects.get(selected='y',owner=self.user,lesson_id=self.id)[0]
        '''
            Just ratchet up the values; no not check anything.
        '''
        if float(self.curStep) == 1.0:   # SOMETHING WRONG WITH ON ENERGY DONE
            self.curStep = '1.5'
        elif float(self.curStep) == 2.0:
            self.curStep = '2.5'
        elif float(self.curStep) == 3.0:
            self.curStep = '3.5'
        elif float(self.curStep) == 4.0:
            self.curStep = '4.5'
        self.save()
        return True

    def onRedoxDone(self,file):
        try:
            LessonProblem.objects.get(lesson_type='lesson6',lesson_id=self.id).delete()
        except:
            pass

        '''
            Here we check the results fo the redox calculations and if they match the precalculated values, 
            the counter will move onto the next step.
        '''
        try:
            struct = structure.models.Structure.objects.get(selected='y',owner=self.user,lesson_id=self.id)
            ws = structure.models.WorkingStructure.objects.filter(structure=struct,selected='y')[0]
            fp = open(struct.location + '/redox-' + ws.identifier + '-modpot.txt', 'r')
            modpot = float(fp.readline())
            fp.close()
            self.LogStep(modpot)
        except:
            lessonprob = LessonProblem(lesson_type='lesson6',lesson_id=self.id,errorstep=3,severity=2,description='There was an error with your redox calculation')
            lessonprob.save()
            return False

#        try:
#            os.stat(struct.location + '/redox-' + ws.identifier + '-modpotref.txt')
#        except OSError, oe:
#            calc_final = False
#        else:
#            print_result = True
#            fp = open(struct.location + '/redox-' + ws.identifier + '-modpotref.txt', 'r')
#            try:
#                modpotref = float(fp.readline())
#            except:
#                print_result = False
#            fp.close()

        self.LogStep(modpot)

        if float(self.curStep) == 1.5 or float(self.curStep) == 1.0 :
            if modpot == 13074.  :
                self.curStep = '2.0'
            else:
                lessonprob = LessonProblem(lesson_type='lesson6',lesson_id=self.id,errorstep=3,severity=2,description='There was an error with your redox calculation')
                lessonprob.save()
                return False
        elif float(self.curStep) == 2.5 or float(self.curStep) == 2.0:
            if modpot == 13089.7 :
                self.curStep = '3.0'
            else:
                lessonprob = LessonProblem(lesson_type='lesson6',lesson_id=self.id,errorstep=3,severity=2,description='There was an error with your redox calculation')
                lessonprob.save()
                return False
        elif float(self.curStep) == 3.5 or float(self.curStep) == 3.0:
            if modpot == 12925.1 :
                self.curStep = '4.0'
            else:
                lessonprob = LessonProblem(lesson_type='lesson6',lesson_id=self.id,errorstep=3,severity=2,description='There was an error with your redox calculation')
                lessonprob.save()
                return False
        elif float(self.curStep) == 4.5 or float(self.curStep) == 4.0:
            if modpot == 12991.1 :
                self.curStep = '5.0'
            else:
                lessonprob = LessonProblem(lesson_type='lesson6',lesson_id=self.id,errorstep=3,severity=2,description='There was an error with your redox calculation')
                lessonprob.save()
                return False

        # Don't check values, just ratchet up step values. 
#        if float(self.curStep) == 2.5 or float(self.curStep) == 2.0 :
#            self.curStep = '3.0'
#        elif float(self.curStep) == 3.5 or float(self.curStep) == 3.0:
#            self.curStep = '4.0'
#        elif float(self.curStep) == 4.5 or float(self.curStep) == 4.0:
#            self.curStep = '5.0'
#        elif float(self.curStep) == 4.5 or float(self.curStep) == 4.0:
#            self.curStep = '5.0'
        self.save()
        self.LogStep('Returning True')
        return True


    def onMDSubmit(self,mdp,filename):
        return True

    def onMDDone(self,mdp):
        return True

    def onLDSubmit(self,postdata):
        return True

    def onLDDone(self,file):
        return True

    def onSGLDSubmit(self,postdata):
        return True

    def onSGLDDone(self,file):
        return True

    def onRMSDSubmit(self,file):
        return True

    def onNATQSubmit(self,postdata):
        return True

    #generates html for the lesson status page
    def generateStatusHtml(self,file):
        step_status_list = []
        step_status_list.append("<tr class='status'><td class='status'>1. Upload File, Build Working Structure, and Calculate Energy: ") 
        step_status_list.append("<tr class='status'><td class='status'>2. Eo for Cv HiPIP: ") 
        step_status_list.append("<tr class='status'><td class='status'>3. Eo for Cv HiPIP with Charged H42: ")
        step_status_list.append("<tr class='status'><td class='status'>4. Eo for Cv HiPIP with H42 Knock Out: ")
        step_status_list.append("<tr class='status'><td class='status'>5. Eo for Cv HiPIP with H42A: ")
        #This will store all the status and the steps, clearing the template of logic
        #And only displaying the status
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson6',lesson_id=self.id)[0]
        except:
            lessonprob = None
        for i in range(self.nSteps):
            if lessonprob and lessonprob.errorstep == math.floor(self.curStep+1) and math.floor(self.curStep) == i:
                step_status_list[i] += ("<font color='red'>Failed</font></td></tr>")
                continue
            elif (float(self.curStep)-0.5) == i and float(self.curStep) % 1.0 == 0.5:
                step_status_list[i] += ("<font color='blue'>Running</font></td></tr>")
                continue
            elif i < float(self.curStep) or len(step_status_list)+1 == float(self.curStep):
                step_status_list[i] += ("<font color='green'>Done</font></td></tr>")
                continue
            elif i+1 > float(self.curStep):
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
        if float(self.curStep) > 5:
            if float(self.curStep) == 5.5:
                htmlcode_list[5] = 2
            else:
                htmlcode_list[5] = 1
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson6',lesson_id=self.id)[0]
            htmlcode_list[lessonprob.errorstep-1] = -1
        except:
            lessonprob = None
        return htmlcode_list
                
