# lesson 
# Enter Descriptions Here
from django.db import models
from django.contrib.auth.models import User
from lessons.models import LessonProblem
from solvation.models import solvationParams
from minimization.models import minimizeParams
from dynamics.models import mdParams, ldParams, sgldParams
from normalmodes.models import nmodeParams
import os
import structure
import lessonaux
import re
from django.template import Context, loader
from django.template.loader import get_template

class Lesson(models.Model):
    class Admin:
        pass
    user = models.ForeignKey(User)
    #Please set the nstep value to the number of steps the lesson will have
    nSteps = models.PositiveIntegerField(default=4)
    curStep = models.DecimalField(default=0,decimal_places=1,max_digits=3)
    
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
        # a template dictionary passing information needed by the template 'generateStatus_Html' 
        template_dict = {}   
        template_dict['status_list'] = []
        
	#Check to see if a lessonproblem exists, and if so store it into lessonprob
	try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson',lesson_id=self.id)[0]
        except:
            lessonprob = None
	#Go through the steps and see which are done,failed, still running. You can leave the code below this
	#comment untouched. Only edit it if you know what you are doing.
        for i in range(self.nSteps):
            # a dictionary to pass the values of function, color and status
            dict = {'num':i+1, 'fun':"function %i" % (i+1), 'color':"", 'status':""}
            # Spesify what functions you want your lesson to do with dict['fun']
            # for example first step "File Uploaded", sencond step "Solvation" and so on
            # if there are more than 4 steps, add lines such as
            # elif i == 5:
            #     dict['fun'] = "xxxx" # and go on 
 
            # Here you specify your own functions
            #if i == 1:
            #    dict['fun'] = "File Uploaded"
            #elif i == 2:
            #    dict['fun'] = "Solvation"
            #elif i == 3:
            #    dict['fun'] = "Minimization"
            #elif i == 4:
            #    dict['fun'] = "MD"

            if lessonprob and lessonprob.errorstep == (self.curStep+1) and self.curStep == i:
                dict['color'] = 'red'
                dict['status'] = 'Failed'  
                template_dict['status_list'].append(dict)
                continue
            elif (self.curStep-.5) == i and self.curStep%1 == 0.5:
                dict['color'] = 'blue'
                dict['status'] = 'Running'  
                template_dict['status_list'].append(dict)
                continue 
            elif i < self.curStep: # or len(step_status_list)+1 == self.curStep:
                dict['color'] = 'green'
                dict['status'] = 'Done'  
                template_dict['status_list'].append(dict)
                continue
            elif i+1 > self.curStep :
                dict['color'] = 'grey'
                dict['status'] = 'N/A'                  
                template_dict['status_list'].append(dict)
                continue
            template_dict['status_list'].append(dict)
        t = get_template('%s/mytemplates/lesson_maker/generateStatus_Html' % charmming_config.charmming_root)        
        return t.render(Context(template_dict))
    
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
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson',lesson_id=self.id)[0]
            htmlcode_list[lessonprob.errorstep-1] = -1
        except:
            lessonprob = None
        return htmlcode_list

