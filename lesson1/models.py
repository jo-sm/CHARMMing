
# lesson 1, upload 1YJP, minimize in vacuum, solvate/neutralize, minimize again, and
# run dynamics
from django.db import models
from django.contrib.auth.models import User
from lessons.models import LessonProblem
from solvation.models import solvationTask
from minimization.models import minimizeTask
from dynamics.models import mdTask
import os, math
import structure, lessonaux, charmming_config
import re

class Lesson1(models.Model):
    # data for lessons (should not be overridden by subclasses)
    # specifying both user and Structure is redundant (since the Structure references the user),
    # but makes it easier to look up all lessons being done by a particular user.
    user = models.ForeignKey(User)
    nSteps = models.PositiveIntegerField(default=4)
    curStep = models.DecimalField(default=0,decimal_places=1,max_digits=3)

    
    def onFileUpload(self):
        try:
            LessonProblem.objects.get(lesson_type='lesson1',lesson_id=self.id).delete()
        except:
            pass
        uploadlog=open("/tmp/uploadlog.txt",'w')
        uploadlog.write("lesson_id: %s\n" % (self.id))
        file = structure.models.Structure.objects.filter(selected='y',lesson_id=self.id)[0]
        #all_segids = file.segids.split() + file.good_het.split() + file.nongood_het.split()
        try:
        #if 2==1:
            #filename1 = '%s/mytemplates/lessons/lesson1/1water.psf' % charmming_config.charmming_root 
            #os.stat(filename1)
            #filename2 = file.location + '1water.psf' + file.stripDotPDB(file.filename) + '.psf'
            #os.stat(filename2)
            filename3 = '%s/mytemplates/lessons/lesson1/1water.crd' % charmming_config.charmming_root
            os.stat(filename3)
            filename4 = file.location+ '/1water.crd'
            os.stat(filename4)
        except:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=1,severity=9,description='The structure you submitted did not upload properly. Check to make sure the structure files are valid.')
            lessonprob.save()
            return False
        if not lessonaux.diffPDBs(file,filename3,filename4):
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=1,severity=9,description='The structure you uploaded was not the correct structure.')
            lessonprob.save()
            return False
        self.curStep = '1'
        self.save()
        return True

    def onEditPDBInfo(self,postdata):
        return True

    def onBuildStructureDone(self,postdata):
        return True

    def onBuildStructureSubmit(self,postdata):
        return True

    def onMinimizeSubmit(self,mp,filename):
        try:
            LessonProblem.objects.get(lesson_type='lesson1',lesson_id=self.id).delete()
        except:
            pass
        #file = structure.models.Structure.objects.filter(selected='y',owner=self.user,lesson_id=self.id)[0]
        #mp = minimizeTask.objects.filter(pdb=file,selected='y')[0]
        ##if filename not in ['new_' + file.stripDotPDB(file.filename) + '-solv']:
        ##    lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=3,severity=2,description='Please minimize the solvated PDB.')
        ##    lessonprob.save()
        ##    return False
        #Now check the sp (solvation parameter) object to make sure they used an RHDO structure
        #with a 15 angstrom distance to the edge of the protein
        
        task=structure.models.Task.objects.get(id=mp.task_ptr_id)
        parent_task=structure.models.Task.objects.get(id=task.parent_id)
        if parent_task.action!='solvation':
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=3,severity=2,description='Please minimize the solvated PDB.')                                                             
            lessonprob.save()
            return False
        if mp.sdsteps != 1000:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=3,severity=2,description='SD steps were not set to 1000.')
            lessonprob.save()
            return False
        if mp.abnrsteps != 1000:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=3,severity=2,description='ABNR steps were not set to 1000.')
            lessonprob.save()
            return False
        if float(mp.tolg) != 0.01:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=3,severity=2,description='TOLG was set to %s, and not 0.01.' % mp.tolg)
            lessonprob.save()
            return False
        if mp.usepbc != 't':
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=3,severity=2,description='Minimization did not use Periodic Boundary Conditions.')
            lessonprob.save()
            return False

        #2.5 Means it is running
        self.curStep = '2.5'
        self.save()
        return True

    def onMinimizeDone(self,mp):
        try:
            lessonprob = LessonProblem.objects.get(lesson_type='lesson1',lesson_id=self.id)
        except:
            lessonprob = None
        #mp = minimizeTask.objects.filter(pdb=file,selected='y')[0]
        #fail = re.compile('Failed')
        if lessonprob:
            self.curStep = '2'
            self.save()
            return False
        #if fail.search(mp.statusHTML):
        if mp.status=='F':
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=3,severity=9,description='The Job did not complete correctly.')
            lessonprob.save()
            self.curStep = '2'
            self.save()
            return False
        else:
            self.curStep = '3'
            self.save()
        return True

    def onSolvationSubmit(self,sp):
        try:
            LessonProblem.objects.filter(lesson_type='lesson1',lesson_id=self.id)[0].delete()
        except:
            pass
        #file = structure.models.Structure.objects.get(selected='y',owner=self.user,lesson_id=self.id)
        #workingstruct = structure.models.WorkingStructure.objects.get(structure=file,selected='y')
        #task=structure.models.Task(workstruct=workingstruct,
        #sp = solvationTask.objects.filter(task_ptr_id=,active='y')[0]
        #pdb_list = lessonaux.getPDBListFromPostdata(file,postdata)

        
        #Now check the sp (solvation parameter) object to make sure they used an RHDO structure
        #with a 15 angstrom distance to the edge of the protein
        if sp.solvation_structure != 'rhdo':
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=2,severity=2,description='You used the wrong solvation structure. To go onto the next step you must use the RHDO structure.')
            lessonprob.save()
            return False
        if (float(sp.xtl_x)<30 or float(sp.xtl_x)>33) or (float(sp.xtl_y)<30 or float(sp.xtl_y)>33) or (float(sp.xtl_z)<30 or float(sp.xtl_z)>33):
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=2,severity=2,description='The wrong radius size was set. Use a value fo 15 to move on in the lesson.')
            lessonprob.save()
            return False
        self.curStep = '1.5'
        self.save()
        return True
       
         #Verify patching
        try:
        #    patch_text = open(file.location + 'new_' + file.stripDotPDB(file.filename) + '-segpatch.txt','r')
        #    patchline = patch_text.readlines()
        #    patch_text.close()
            #if "SEGMENT wat :: generate setu wat first NONE last NONE\n" not in patchline:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=2,severity=2,description='The patching was not done correctly.')
            lessonprob.save()
            return False
        except:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=2,severity=2,description='The patching was not done correctly.')
            lessonprob.save()
            return False

    def onSolvationDone(self,sp):
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson1',lesson_id=self.id)[0]
        except:
            lessonprob = None
        #sp = solvationTask.objects.filter(pdb=file,active='y')[0]
        #fail = re.compile('Failed')
        if lessonprob:
            self.curStep = '1'
            self.save()
            return False
        #if fail.search(sp.statusHTML):
        if sp.status=='F':
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=2,severity=9,description='The job did not complete correctly.')
            lessonprob.save()
            self.curStep = '1'
            self.save()
            return False
        else:
            self.curStep = '2'
            self.save()
        return True

    def onNMASubmit(self,postdata):
        return True

    def onNMADone(self,file):
        return True

    def onMDSubmit(self,mdp,filename):
        #Clear any old lessonproblems
        try:
            LessonProblem.objects.get(lesson_type='lesson1',lesson_id=self.id).delete()
        except:
            pass
        #file = structure.models.Structure.objects.get(selected='y',owner=self.user,lesson_id=self.id)[0]
        #mdp = mdTask.objects.filter(pdb=file,selected='y')[0]
        ##if filename not in ['new_' + file.stripDotPDB(file.filename) + '-min.pdb']:
        ##    lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=2,description='Please run dynamics on the minimized PDB (-min).')
        ##    lessonprob.save()
        ##    return False
        task=structure.models.Task.objects.get(id=mdp.task_ptr_id)
        parent_task=structure.models.Task.objects.get(id=task.parent_id)
        if parent_task.action!='minimization':
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=4,description='Please use the minimized and solvated PDB coordinates.')                                                             
            lessonprob.save()
            return False
        if mdp.ensemble != 'heat':
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=2,description='You used an equilibration calculation instead of a heating one. Please use heating to continue.')
            lessonprob.save()
            return False
        if mdp.nstep != 1000:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=2,description='Please set the number of steps to 1000 to continue.')
            lessonprob.save()
            return False
        if float(mdp.firstt) != 210.15:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=2,description='Please set the starting temperature to 210.15 K to continue.')
            lessonprob.save()
            return False
        if float(mdp.teminc) != 10:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=2,description='Please set the temperature increment to 10 K to continue.')
            lessonprob.save()
            return False
        if float(mdp.ihtfrq) != 100:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=2,description='Please set the heating frequency to 100 steps.')
            lessonprob.save()
            return False
        if float(mdp.finalt) != 310.15:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=2,description='Please set the final temperature to physiological temperature (310.15 K).')
            lessonprob.save()
            return False
        if float(mdp.tbath) != 310.15:
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=2,description='Please set the temperature bath physiological temperature (310.15 K).')
            lessonprob.save()
            return False
        self.curStep = '3.5'
        self.save()
        return True

    def onMDDone(self,mdp):
        try:
            lessonprob = LessonProblem.objects.get(lesson_type='lesson1',lesson_id=self.id)
        except:
            lessonprob = None
        #mdp = mdTask.objects.filter(pdb=file,selected='y')[0]
        #fail = re.compile('Failed')
        if lessonprob:
            return False
        #if fail.search(mdp.statusHTML):
        if mdp.status=='F':
            lessonprob = LessonProblem(lesson_type='lesson1',lesson_id=self.id,errorstep=4,severity=9,description='The job did not complete correctly.')
            lessonprob.save()
            return False
        else:
            self.curStep = '4'
            self.save()
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

    def onEnergySubmit(self,postdata):
        return True

    def onEnergyDone(self,finale):
        return True
    
    #generates html for the lesson status page
    def generateStatusHtml(self,file):
        step_status_list = []
        step_status_list.append("<tr class='status'><td class='status'>1. File Uploaded: ") 
        step_status_list.append("<tr class='status'><td class='status'>2. Solvation: ") 
        step_status_list.append("<tr class='status'><td class='status'>3. Minimization: ") 
        step_status_list.append("<tr class='status'><td class='status'>4. MD: ") 
        #This will store all the status and the steps, clearing the template of logic
        #And only displaying the status
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson1',lesson_id=self.id)[0]
        except:
            lessonprob = None
        for i in range(self.nSteps):
            if lessonprob and lessonprob.errorstep == math.floor(self.curStep+1) and math.floor(self.curStep) == i:
                step_status_list[i] += ("<a class='failed' href='javascript:open_failure();'>Failed</a></td></tr>")
                continue
            elif (float(self.curStep)-0.5) == i and float(self.curStep) % 1.0 == 0.5:
                step_status_list[i] += ("<span class='running'>Running</span></td></tr>")
                continue
            elif i < float(self.curStep) or len(step_status_list)+1 == float(self.curStep):
                step_status_list[i] += ("<span class='done'>Done</span></td></tr>")
                continue
            elif i+1 > float(self.curStep):
                step_status_list[i] += ("<span class='inactive'>N/A</span></td></tr>")
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
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson1',lesson_id=self.id)[0]
            htmlcode_list[lessonprob.errorstep-1] = -1
        except:
            lessonprob = None
        return htmlcode_list
                
