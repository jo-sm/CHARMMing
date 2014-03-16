
# lesson 1, upload 1YJP, minimize in vacuum, solvate/neutralize, minimize again, and
# run dynamics
from django.db import models
from django.contrib.auth.models import User
from lessons.models import LessonProblem
from solvation.models import solvationTask
from minimization.models import minimizeTask
from dynamics.models import mdTask, dynamicsTask
import os, re, math
import structure, lessonaux, charmming_config

class Lesson2(models.Model):
    # data for lessons (should not be overridden by subclasses)
    # specifying both user and PDBFile is redundant (since the PDBFile references the user),
    # but makes it easier to look up all lessons being done by a particular user.
    user = models.ForeignKey(User)
    nSteps = models.PositiveIntegerField(default=6)
    curStep = models.DecimalField(default=0,decimal_places=1,max_digits=3)


    def onFileUpload(self):
        try:
            LessonProblem.objects.get(lesson_type='lesson2',lesson_id=self.id).delete()
        except:
            pass
        file = structure.models.Structure.objects.get(selected='y',owner=self.user,lesson_id=self.id)
        #all_segids = file.segids.split() + file.good_het.split() + file.nongood_het.split()
        try:
        #if 2==1:
            #filename1 = '%s/mytemplates/lessons/lesson2/par_all27_cmap_chol.prm' % charmming_config.charmming_root
            #os.stat(filename1)
            #filename2 = file.location + 'prm-' + file.stripDotPDB(file.filename) + '.prm'
            #os.stat(filename2)
            #filename3 = '%s/mytemplates/lessons/lesson2/top_all27_prot_lipid.rtf' % charmming_config.charmming_root
            #os.stat(filename3)
            #filename4 = file.location + 'rtf-' + file.stripDotPDB(file.filename) + '.rtf'
            #os.stat(filename4)
            filename1 = '%s/mytemplates/lessons/lesson2/1zhy.pdb' % charmming_config.charmming_root                                                                                                  
            os.stat(filename1)
            filename2 = file.location + '/1zhy.pdb'
            os.stat(filename2)
        except:
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=1,severity=9,description='The PDB you submitted did not upload properly. Check to make sure the PDB.org ID was valid.')
            lessonprob.save()
            return False
        #if not lessonaux.diffPDBs(file,filename1,filename2) or not lessonaux.diffPDBs(file,filename3,filename4):
        #    lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=1,severity=9,description='The topology/parameter files you uploaded were not correct.')
        #    lessonprob.save()
        #    return False
        #try:
        #    for segid in all_segids:
        #        filename1 = '%s/mytemplates/lessons/lesson2/1zhy-%s.pdb' % (charmming_config.charmming_root,segid)
        #        os.stat(filename1)
        #        filename2 = file.location + 'new_' + file.stripDotPDB(file.filename) + '-' + segid + '.pdb'
        #        os.stat(filename2)
        if not lessonaux.diffPDBs(file,filename1,filename2):
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=1,severity=9,description='The PDB you uploaded was not the correct PDB.')
            lessonprob.save()
            return False
        #except:
        #    lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=1,severity=9,description='The PDB you submitted did not upload properly. Check to make sure the PDB.org ID was valid.')
        #    lessonprob.save()
            return False
        self.curStep = '1'
        self.save()
        return True

    def onEditPDBInfo(self,postdata):
        return True

    def onMinimizeSubmit(self,mp,filename):
        try:
            LessonProblem.objects.filter(lesson_type='lesson2',lesson_id=self.id)[0].delete()
        except:
            pass
        #file = structure.models.Structure.objects.filter(selected='y',owner=self.user,lesson_id=self.id)[0]
        #mp = minimizeTask.objects.filter(pdb=file,selected='y')[0]
        #if filename not in ['new_' + file.stripDotPDB(file.filename) + '-neutralized']:
        #    lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=3,severity=2,description='Please minimize the neutralized PDB.')
        #    lessonprob.save()
        #    return False
        #Now check the sp (solvation parameter) object to make sure they used an RHDO structure
        #with a 15 angstrom distance to the edge of the protein
        if mp.sdsteps != 1000:
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=4,severity=2,description='SD steps were not set to 1000.')
            lessonprob.save()
            return False
        if mp.abnrsteps != 1000:
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=4,severity=2,description='ABNR steps were not set to 1000.')
            lessonprob.save()
            return False
        if float(mp.tolg) != .01:
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=4,severity=2,description='TOLG was not set not 0.01.')
            lessonprob.save()
            return False
        if mp.usepbc != 't':
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=4,severity=2,description='Minimization did not use Periodic Boundary Conditions.')
            lessonprob.save()
            return False

        #2.5 Means it is running
        self.curStep = '3.5'
        self.save()
        return True

    def onMinimizeDone(self,mp):
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson2',lesson_id=self.id)[0]
        except:
            lessonprob = None
        #mp = minimizeTask.objects.filter(pdb=file,selected='y')[0]
        fail = re.compile('Failed')
        
        if lessonprob:
            self.curStep = '2'
            self.save()
            return False
        
        #if fail.search(mp.statusHTML):
        if mp.status=='F':
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=4,severity=9,description='The Job did not complete correctly.')
            lessonprob.save()
            self.curStep = '3'
            self.save()
            return False
        else:
            self.curStep = '4'
            self.save()
        return True

    def onSolvationSubmit(self,sp):
        #l2log=open("/tmp/l2.log",'w')
        try:
            LessonProblem.objects.get(lesson_type='lesson2',lesson_id=self.id).delete()
        except:
            pass
        #file = structure.models.Structure.objects.filter(selected='y',owner=self.user,lesson_id=self.id)[0]
        #sp = solvationTask.objects.filter(pdb=file,active='y')[0]

        #This ensures the user selected the right PDBs
        ##pdb_list = lessonaux.getPDBListFromPostdata(file,postdata)
        ##acceptable_pdb_list = ['new_' + file.stripDotPDB(file.filename) + '-a-pro.pdb','new_' + file.stripDotPDB(file.filename) + '-a-pro-final.pdb','new_' + file.stripDotPDB(file.filename) + '-a-goodhet.pdb', \
        ##                       'new_' + file.stripDotPDB(file.filename) + '-a-goodhet-final.pdb','new_' + file.stripDotPDB(file.filename) + '-a-het.pdb','new_' + file.stripDotPDB(file.filename) + '-a-het-final.pdb']
        ##for pdb in pdb_list:
        ##    if pdb not in acceptable_pdb_list:
        ##        lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=3,severity=2,description='Please select all of the initial segments. Do not use a segment that has had a calculation done on it.')
        ##        lessonprob.save()
        ##        return False
        #Now check the sp (solvation parameter) object to make sure they used an RHDO structure
        #with a 15 angstrom distance to the edge of the protein
        #l2log.write("solvation structure is: %s\n" % (sp.solvation_structure))
        if sp.solvation_structure != 'rhdo':
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=3,severity=2,description='You used the wrong solvation structure. To go onto the next step you must use the RHDO structure.')
            lessonprob.save()
            return False
        #if sp.no_pref_radius != 10:
        if (float(sp.xtl_x)<80 or float(sp.xtl_x)>83) or (float(sp.xtl_y)<80 or float(sp.xtl_y)>83) or (float(sp.xtl_z)<80 or float(sp.xtl_z)>83):
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=3,severity=2,description='The wrong radius size was set. Use a value fo 15 to move on in the lesson.')
            lessonprob.save()
            return False
        if sp.salt != "POT":
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=3,severity=2,description='The wrong salt was used. Please use potassium chloride.')
            lessonprob.save()
            return False
        if float(sp.concentration) != 0.15:
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=3,severity=2,description='The wrong concentration was used. Please use a concentration of .15.')
            lessonprob.save()
            return False
        if float(sp.ntrials) != 1:
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=3,severity=2,description='The wrong number of trials were used. Please use 1 trial.')
            lessonprob.save()
            return False
        self.curStep = '2.5'
        self.save()
        return True

    def onSolvationDone(self,sp):
        try:
            lessonprob = LessonProblem.objects.get(lesson_type='lesson2',lesson_id=self.id)
        except:
            lessonprob = None
        #sp = solvationTask.objects.filter(pdb=file,active='y')[0]
        #fail = re.compile('Failed')
        if lessonprob:
            self.curStep = '1'
            self.save()
            return False
        #if fail.search(sp.statusHTML):
        if sp.status == 'F':
            lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=3,severity=9,description='The job did not complete correctly.')
            lessonprob.save()
            self.curStep = '2'
            self.save()
            return False
        else:
            self.curStep = '3'
            self.save()
        return True

    def onNMASubmit(self,postdata):
        return True

    def onNMADone(self,file):
        return True

    def onMDSubmit(self,mdp,filename):
        l2mdlog=open("/tmp/l2md.log","w")
        l2mdlog.write("starting md submit with mdtask %s\n" % (mdp))
        #Clear any old lessonproblems
        try:
            LessonProblem.objects.get(lesson_type='lesson2',lesson_id=self.id).delete()
        except:
            pass
        #file = structure.models.Structure.objects.filter(selected='y',owner=self.user,lesson_id=self.id)[0]
        #mdp = mdTask.objects.filter(pdb=file,selected='y')[0]

        if float(self.curStep) == float(4.0):
            l2mdlog.write("step is 4.0\n")
            #if filename not in ['new_' + file.stripDotPDB(file.filename) + '-min.pdb']:
            #    lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=2,description='Please run dynamics on the minimized PDB (-min).')
            #    lessonprob.save()
            #    return False
            parent_task = structure.models.Task.objects.get(id=structure.models.Task.objects.get(id=dynamicsTask.objects.get(id=mdp.dynamicstask_ptr_id).task_ptr_id).parent_id)
            l2mdlog.write("parent task id and action are: %s and %s\n" % (parent_task.id, parent_task.action))
            l2mdlog.close()
            if parent_task.action !='minimization':
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=2,description='Please run heating dynamics on the PDB from the minimization.')
                lessonprob.save()
                return False
            if mdp.ensemble != 'heat':
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=2,description='You used an equilibration calculation instead of a heating one. Please use heating to continue.')
                lessonprob.save()
                return False
            if mdp.nstep != 1000:
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=2,description='Please set the number of steps to 1000 to continue.')
                lessonprob.save()
                return False
            if float(mdp.firstt) != 210.15:
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=2,description='Please set the starting temperature to 100 K to continue.')
                lessonprob.save()
                return False
            if float(mdp.teminc) != 10:
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=2,description='Please set the temperature increment to 10K to continue.')
                lessonprob.save()
                return False
            if float(mdp.ihtfrq) != 100:
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=2,description='Please set the heating frequency to 20 steps.')
                lessonprob.save()
                return False
            if float(mdp.finalt) != 310.15:
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=2,description='Please set the final temperature to physiological temperature (310.15 K).')
                lessonprob.save()
                return False
            if float(mdp.tbath) != 310.15:
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=2,description='Please set the temperature bath physiological temperature (310.15 K).')
                lessonprob.save()
                return False
            self.curStep = '4.5'
            self.save()
            return True
        elif float(self.curStep) == float(5.0):
            parent_task = structure.models.Task.objects.get(id=structure.models.Task.objects.get(id=dynamicsTask.objects.get(id=mdp.dynamicstask_ptr_id).task_ptr_id).parent_id)
            l2mdlog.write("parent task id and action are: %s and %s\n" % (parent_task.id, parent_task.action))
            l2mdlog.close()
            if parent_task.action != 'md':
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=6,severity=2,description='Please run equilibration dynamics on the PDB from the MD heating run (-md).')
                lessonprob.save()
                return False 
            if mdp.ensemble != 'npt':
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=6,severity=2,description='You used an incorrect type of run. Please use NPT to continue.')
                lessonprob.save()
                return False
            if mdp.nstep != 1000:
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=6,severity=2,description='Please set the number of steps to 1000 to continue.')
                lessonprob.save()
                return False
            if float(mdp.temp) != 310.15:
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=6,severity=2,description='Please set the simulation temperature  to 310.15 K.')
                lessonprob.save()
                return False
            self.curStep = '5.5'
            self.save()
            return True
        else:
            # to-do, MD done at the wrong time
            return False

    def onMDDone(self,mdp):
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson2',lesson_id=self.id)[0]
        except:
            lessonprob = None
        #mdp = mdTask.objects.filter(pdb=file,selected='y')[0]
        #fail = re.compile('Failed')
        if lessonprob:
            return False
        #if fail.search(mdp.statusHTML):
        if  mdp.status == 'F':
            if float(self.curStep) == 4.5:
               lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=5,severity=9,description='The job did not complete correctly.')
            else:
               lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=6,severity=9,description='The job did not complete correctly.')
            lessonprob.save()
            return False
        else:
            if float(self.curStep) == 4.5:
                self.curStep = '5'
            else:
                self.curStep = '6'
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
        try:
            LessonProblem.objects.get(lesson_type='lesson2',lesson_id=self.id).delete()
        except:
            pass

        if float(self.curStep) == 1.0:
            # nothing really to validate here...
            self.curStep='1.5'
            self.save()

        return True

    def onEnergyDone(self,et):
        if float(self.curStep) == 1.5:
            if float(et.finale) < float(40690) or float(et.finale) > float(40695):
                #t4log.write("saving error finale: %s\n" % str(float(et.finale)))
                lessonprob = LessonProblem(lesson_type='lesson2',lesson_id=self.id,errorstep=2,severity=9,description="Your energy was wrong (should be ~ 40692, you got %s). Make sure to use all default settings on the energy page, and to use the custom parameter and topology files provided in this lesson" % et.finale)
                lessonprob.save()
                self.curStep = '1'
                self.save()
                return False
        self.curStep = '2'
        self.save()
        return True

    #generates html for the lesson status page
    def generateStatusHtml(self,file):
        step_status_list = []
        step_status_list.append("<tr class='status'><td class='status'>1. File Uploaded: ") 
        step_status_list.append("<tr class='status'><td class='status'>2. Energy: ") 
        step_status_list.append("<tr class='status'><td class='status'>3. Solvation: ") 
        step_status_list.append("<tr class='status'><td class='status'>4. Minimization: ") 
        step_status_list.append("<tr class='status'><td class='status'>5. MD Heating: ") 
        step_status_list.append("<tr class='status'><td class='status'>6. MD Equilibration: ") 
        #This will store all the status and the steps, clearing the template of logic
        #And only displaying the status
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson2',lesson_id=self.id)[0]
        except:
            lessonprob = None
        for i in range(self.nSteps):
            if lessonprob and lessonprob.errorstep == math.floor(self.curStep+1) and math.floor(self.curStep) == i:
                step_status_list[i] += ("<a class='failed' href='javascript:open_failure()'>Failed</span></td></tr>")
                continue
            elif (float(self.curStep)-0.5) == i and float(self.curStep) % 1 == 0.5:
                step_status_list[i] += ("<span class='running'>Running</span></td></tr>")
                continue
            elif i < float(self.curStep):
                step_status_list[i] += ("<span class='done'>Done</span></td></tr>")
                continue
            elif i + 1 > float(self.curStep):
                step_status_list[i] += ("<span class='inactive'>N/A</span></td></tr>")
                continue
        return step_status_list

    #Returns a list where each index corresponds to lesson progress
    #on the display lesson page
    def getHtmlStepList(self):
        # 2 is running
        # 1 is successfully done
        # 0 is not started
        # -1 is error

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
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson2',lesson_id=self.id)[0]
            htmlcode_list[lessonprob.errorstep-1] = -1
        except:
            lessonprob = None
        return htmlcode_list
#Exists so that building the structure doesn't explode
    def onBuildStructureDone(self,postdata):
        return True

    def onBuildStructureSubmit(self,postdata):
        return True
