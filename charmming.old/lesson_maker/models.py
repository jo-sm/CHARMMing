#this is the lesson maker
#it makes lessons by analyzing usage patterns
from django.db import models
#from solvation.models import solvationTask
#from minimization.models import minimizeTask
#from dynamics.models import mdTask
#we don't know which of these we need yet.
import structure

class Lesson(models.Model):
    # data for lessons (should not be overridden by subclasses)
    # specifying both user and Structure is redundant (since the Structure references the user),
    # but makes it easier to look up all lessons being done by a particular user.
    #presumably you could calculate this out from the count of Steps it currently has, but it doesn't allow for deletio
    finished = models.BooleanField(default=False,null=False) #update this field at the end
    structure = models.ForeignKey(structure.models.Structure,null=False) #set this at the start...
    filename = models.CharField(max_length=25,null=True) #this keeps track of the filename of the thing the user uploaded
    name = models.CharField(max_length=25,null=True) #this keeps track of the name, it's for the create script
    intro_text = models.TextField(null=True) #this tracks the intro text. TextField cause this can get huge.
    title = models.CharField(max_length=75,null=True) #this tracks the title. We shouldn't need more than 75 chars, and should do input sanitization on entry.
    #we need this to compare with the file later on...

    def onTaskDone(self,inp_task):
        #we call the exact same code over and over for every task, so we're making it generic
        task = inp_task
        #we can pass the task directly in structure.models...we don't have to care
        if task.status=="F":
            return False
        else:
            file_step = Step()
            file_step.lesson = Lesson.objects.get(structure=task.workstruct.structure)
            #now we do have to care, so we fetch the object since we are calling this from
            #an empty Lesson()
            file_step.step = len(Step.objects.filter(lesson=self)) + 1 #make this generic, so we start at 0
            file_step.task = task
            file_step.type = task.action
            file_step.structure = None #avoid data redundancy!
            file_step.workingstructure = None
            file_step.save()
            self.save() #just in case...
            return True

    def onFileUpload(self):
        #this one is called in structure.views so we don't have to worry about circulat import crap
        file_step = Step()
        file_step.lesson = self
        file_step.step = len(Step.objects.filter(lesson=self)) + 1 #make this generic, so we start at 0
        #we only need whole number steps, but we will leave them in the format that is compatible with lessons currently available
        file_step.task = None
        file_step.type = "fileupload"
        file_step.save()
        self.save() #structure is defined by this point so we don't have to care
        self.structure.lessonmaker_active = True
        self.structure.save()
        return True
        #we'll pass this step the structure from the Lesson object

    def onBuildStructureDone(self,postdata):
        return True
#this one isn't called either

    def onBuildStructureSubmit(self):
        #this one actually gets the structure built before it gets called, so we don't have to worry about it
        #we sadly don't actually have the task object otherwise I'd save that
        file_step = Step()
        file_step.lesson = self
        #we also do this on the page directly so we don't have to care just yet
        file_step.step = len(Step.objects.filter(lesson=self)) + 1 #make this generic, so we start at 0
        file_step.task = None
        file_step.type = "buildstruct"
        #we'll pass this step a WorkingStructure later.
        file_step.save() #highly doubt this does anything but why not
        self.save()
        return True

#    def onMinimizeSubmit(self,mp,filename):
#        return True
#    #we don't care about minimization until it succeeds
#
#    def onMinimizeDone(self,mp):
#        return self.onTaskDone(mp)
#
#    def onSolvationSubmit(self,sp):
#        return True
#    #again we don't care
#
#    def onSolvationDone(self,sp):
#        return self.onTaskDone(sp)
#
#    def onNMASubmit(self,postdata):
#        return True
#
#    def onNMADone(self,nmodetask):
#        return self.onTaskDone(nmodetask)
#
#    def onMDSubmit(self,mdp,filename):
#        return True
#    #don't care about this one
#
#    def onMDDone(self,mdp):
#        return self.onTaskDone(mdp)
#
#    def onLDSubmit(self,postdata):
#        return True
#    #don't care here either
#
#    def onLDDone(self,ldtask):
#        return self.onTaskDone(ldtask)
#
#    def onSGLDSubmit(self,postdata):
#        return True
#
#    def onSGLDDone(self,sgldtask):
#        return self.onTaskDone(sgldtask)
#
#    def onRMSDSubmit(self,file):
#        return True
#    #can't do much with this one yet
#
#    def onEnergySubmit(self,postdata):
#        return True
#
#    def onEnergyDone(self,finale):
#        return self.onTaskDone(finale)
#
#    def onRedoxDone(self,redoxtask):
#        return self.onTaskDone(redoxtask)

#    #generates html for the lesson status page
#use this later for the finished builder
#    def generateStatusHtml(self,file):
#        step_status_list = []
#        step_status_list.append("<tr class='status'><td class='status'>1. File Uploaded: ") 
#        step_status_list.append("<tr class='status'><td class='status'>2. Solvation: ") 
#        step_status_list.append("<tr class='status'><td class='status'>3. Minimization: ") 
#        step_status_list.append("<tr class='status'><td class='status'>4. MD: ") 
#        #This will store all the status and the steps, clearing the template of logic
#        #And only displaying the status
#        try:
#            lessonprob = LessonProblem.objects.filter(lesson_type='lesson1',lesson_id=self.id)[0]
#        except:
#            lessonprob = None
#        for i in range(self.nSteps):
#            if lessonprob and lessonprob.errorstep == math.floor(self.curStep+1) and math.floor(self.curStep) == i:
#                step_status_list[i] += ("<a class='failed' href='javascript:open_failure();'>Failed</a></td></tr>")
#                continue
#            elif (float(self.curStep)-0.5) == i and float(self.curStep) % 1.0 == 0.5:
#                step_status_list[i] += ("<span class='running'>Running</span></td></tr>")
#                continue
#            elif i < float(self.curStep) or len(step_status_list)+1 == float(self.curStep):
#                step_status_list[i] += ("<span class='done'>Done</span></td></tr>")
#                continue
#            elif i+1 > float(self.curStep):
#                step_status_list[i] += ("<span class='inactive'>N/A</span></td></tr>")
#                continue
#            
#        return step_status_list
#
#    #Returns a list where each index corresponds to lesson progress
#    #on the display lesson page
#    def getHtmlStepList(self):
#        #2 is running
#        #1 is successfully done
#        #0 is not started
#        #-1 is error
#        htmlcode_list = []
#        for step in range(self.nSteps):
#            htmlcode_list.append(0)
#        if float(self.curStep) > 0:
#            htmlcode_list[0] = 1
#        if float(self.curStep) > 1:
#            if float(self.curStep) == 1.5:
#                htmlcode_list[1] = 2
#            else:
#                htmlcode_list[1] = 1
#        if float(self.curStep) > 2:
#            if float(self.curStep) == 2.5:
#                htmlcode_list[2] = 2
#            else:
#                htmlcode_list[2] = 1
#        if float(self.curStep) > 3:
#            if float(self.curStep) == 3.5:
#                htmlcode_list[3] = 2
#            else:
#                htmlcode_list[3] = 1
#        try:
#            lessonprob = LessonProblem.objects.filter(lesson_type='lesson1',lesson_id=self.id)[0]
#            htmlcode_list[lessonprob.errorstep-1] = -1
#        except:
#            lessonprob = None
#        return htmlcode_list
#                
class Step(models.Model):
    lesson = models.ForeignKey(Lesson)
    step = models.DecimalField(default=3,decimal_places=1,max_digits=3)
    task = models.ForeignKey(structure.models.Task,null=True) #Task will be interpreted depending on the type...
    #We will have to update Task stuff for the job sched rebuild
    type = models.CharField(max_length=20,null=False) #no default value cause this will get modified immediately
    running_text = models.TextField(null=True) #text to display while it's running
    done_text = models.TextField(null=True) #text to display when it's done, both of these are longwinded.

class Objective(models.Model):
    lesson = models.ForeignKey(Lesson)
    obj_num = models.IntegerField(null=False) #this gives the objective number because we don't want them to be out of order
    obj_text = models.CharField(max_length=240,null=False) #this should be spawned on creation of an Objective Object
