# lesson {{lesson.name}} # {{lesson.title}}
from django.db import models
from django.contrib.auth.models import User
from lessons.models import LessonProblem
from solvation.models import solvationTask
from minimization.models import minimizeTask
from dynamics.models import mdTask, ldTask, sgldTask
from normalmodes.models import nmodeTask
import os
import math
import structure
import lessonaux
import re
from django.template import Context, loader
from django.template.loader import get_template
import charmming_config
{%comment%}blankme
#Formatting notes:
#    lesson.name is just an int with the lesson number
#    lesson.title is the lesson descriptor
#    stepcount is len(steps)
#    lesson.filename is the filename (with extension!) of the structure file we
#                       copy from the lesson maker upload
#  
#IMPORTANT NOTE!
#Django doesn't know what indentation is. Keep the indentation in any includes you have
#and then put the include at the lowest indentation level
#that way django will just render the tabs, and you stop worrying.

{%endcomment%}blankme
class Lesson{{lesson.name}}(models.Model):
    user = models.ForeignKey(User)
    nSteps = models.PositiveIntegerField(default={{stepcount}})
    curStep = models.DecimalField(default=0,decimal_places=1,max_digits=3)
    
    def onFileUpload(self):
{%include 'new_lesson_maker/lessonprob_delete_tpl.py' %}
        file = structure.models.Structure.objects.get(selected='y',owner=self.user,lesson_id=self.id)
        #this should never return more than 1...if it does we've got some serious DB trouble
        try:
            filename1 = '%s/mytemplates/lessons/lesson{{lesson.name}}/{{lesson.filename}}' % charmming_config.charmming_root
            os.stat(filename1)
            filename2 = file.location + '/{{lesson.filename}}'
            os.stat(filename2)
        except:
            lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep=1,severity=9,description='The PDB you submitted did not upload properly. Check to make sure the PDB.org ID was valid.')
            lessonprob.save()
            return False

        if not lessonaux.diffPDBs(file,filename1,filename2):
            lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep=1,severity=9,description='The PDB you uploaded was not the correct PDB.')
            lessonprob.save()
            return False
        self.curStep = '1' {%comment%}If this is ever NOT 1, we're in trouble.{%endcomment%}
        self.save()
        return True

    def onBuildStructureSubmit(self,task):
        self.curStep = '2' {%comment%}Same deal as in file upload.{%endcomment%}
        self.save()
        return True

    def onBuildStructureDone(self,task):
        return True

    def onMinimizeSubmit(self,mp,filename):
{%comment%}In order to make this truly generic, we'll check for multiple steps even if there aren't any{%endcomment%}
{%comment%}blankme
#task_dict.minimization_tasks just sends a bunch of minimizationTask objects into here, so you can check each of those
#for the parameters we care about. However, that's slot 0. Slot 1 carries the step number, so we
#don't lose track of that stuff.
#Because Django is silly, slot 2 carries slot 1, with 0.5 added to it, to make Django's filters 
#not crash.
{%endcomment%}blankme
    {%if task_dict.minimization_tasks%}blankme
{%include 'new_lesson_maker/lessonprob_delete_tpl.py' %}
    {%for minimization_task in task_dict.minimization_tasks%}blankme
        if float(self.curStep) == float({{minimization_task.1|add:"1"}}):
            {% if minimization_task.0.parent.action %}blankme
            parent_task_action = mp.parent.action
            if parent_task_action != '{{minimization_task.0.parent.action}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{minimization_task.1|add:"1"}},severity=2,description='Please run {{minimization_task.0.action}} on the coordinates from {{minimization_task.parent.action}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
            if mp.sdsteps != {{minimization_task.0.sdsteps}}
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{minimization_task.1|add:"1"}},severity=2,description='SD steps were not set to {{minimization_task.0.sdsteps}}.')
                lessonprob.save()
                return False
            if mp.abnrsteps != {{minimization_task.0.abnrsteps}}:
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{minimization_task.1|add:"1"}},severity=2,description='ABNR steps were not set to {{minimization_task.0.abnrsteps}}.')
                lessonprob.save()
                return False
            if float(mp.tolg) != float({{minimization_task.0.tolg}}):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{minimization_task.1|add:"1"}},severity=2,description='TOLG was not set to {{minimization_task.0.tolg}}.')
                lessonprob.save()
                return False
            if mp.usepbc != '{{minimization_task.0.usepbc}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{minimization_task.1|add:"1"}},severity=2,description='PBC use was not set to {{minimization_task.0.usepbc}}.')
                lessonprob.save()
                return False
            if mp.useqmmm != '{{minimization_task.0.useqmmm}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{minimization_task.1|add:"1"}},severity=2,description='QM/MM use was not set to {{minimization_task.0.useqmmm}}.')
                lessonprob.save()
                return False

            self.curStep = '{{minimization_task.2}}'
            self.save()
            return True
    {%endfor%}blankme
    {%endif%}blankme
        return True

    def onMinimizeDone(self,mdp):
    {%if task_dict.minimization_tasks %}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        if mdp.status == 'F':
            {%for task in task_dict.minimization_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=9,description='The job did not complete correctly.')
                lessonprob.save()
                self.curStep = '{{task.1}}'
                self.save()
                return False
            {%endfor%}blankme
        else:
            {%for task in task_dict.minimization_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                self.curStep = '{{task.1|add:"1"}}'
                self.save()
                return True
            {%endfor%}blankme
            {%endif%}blankme
        return True

    def onSolvationSubmit(self,sp):
    {%if task_dict.solvation_tasks%}blankme
{%include 'new_lesson_maker/lessonprob_delete_tpl.py' %}
    {%for solvation_task in task_dict.solvation_tasks%}blankme
        if float(self.curStep) == float({{solvation_task.1|add:"1"}}):
            {% if solvation_task.0.parent %}blankme
            parent_task_action = sp.parent.action
            if parent_task_action != '{{solvation_task.0.parent.action}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{solvation_task.1|add:"1"}},severity=2,description='Please run {{solvation_task.0.action}} on the coordinates from {{solvation_task.parent.action}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
            if sp.solvation_structure != '{{solvation_task.0.solvation_structure}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{solvation_task.1|add:"1"}},severity=2,description='Solvation Structure was not set to {{solvation_task.0.solvation_structure}}.')
                lessonprob.save()
                return False
            #I'm using +/- 2 because I don't really know what this value means, but also Django can't multiply.
            if (float(sp.xtl_x)<float({{solvation_task.0.xtl_x}} - 2) or float(sp.xtl_x) > float({{solvation_task.0.xtl_x}} + 2)) or (float(sp.xtl_y)<float({{solvation_task.0.xtl_y}}-2) or float(sp.xtl_y)>float({{solvation_task.0.xtl_y}}+2)) or (float(sp.xtl_z)<float({{solvation_task.0.xtl_z}}-2) or float(sp.xtl_z)>float({{solvation_task.0.xtl_z}}+2)):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{solvation_task.1|add:"1"}},severity=2,description'The wrong radius size was set. Use a value of {{solvation_task.0.sp_radius}} instead.')
                lessonprob.save()
                return False
            if sp.salt != '{{solvation_task.0.salt}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{solvation_task.1|add:"1"}},severity=2,description='The wrong salt was used.') #Modify this line for your salt...
                lessonprob.save()
                return False
            if float(sp.concentration) != {{solvation_task.0.concentration}}:
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{solvation_task.1|add:"1"}},severity=2,description='PBC use was not set to {{solvation_task.0.usepbc}}.')
                lessonprob.save()
                return False
            if float(sp.ntrials) != {{solvation_task.0.ntrials}}:
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{solvation_task.1|add:"1"}},severity=2,description='QM/MM use was not set to {{solvation_task.0.useqmmm}}.')
                lessonprob.save()
                return False

            self.curStep = '{{solvation_task.2}}'
            self.save()
            return True
    {%endfor%}blankme
    {%endif%}blankme
        return True

    def onSolvationDone(self,sp):
    {%if task_dict.solvation_tasks %}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        if sp.status == 'F':
            {%for task in task_dict.solvation_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=9,description='The job did not complete correctly.')
                lessonprob.save()
                self.curStep = '{{task.1}}'
                self.save()
                return False
            {%endfor%}blankme
        else:
            {%for task in task_dict.solvation_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                self.curStep = '{{task.1|add:"1"}}'
                self.save()
                return True
            {%endfor%}blankme
            {%endif%}blankme
        return True

    def onEnergySubmit(self,et):
        {%if task_dict.energy_tasks%}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        #there's a couple things we can verify here
    {%for energy_task in task_dict.energy_tasks%}blankme
        if float(self.curStep) == float({{energy_task.1|add:"1"}}):
            {% if energy_task.0.parent %}blankme
            parent_task_action = et.parent.action
            if parent_task_action != '{{energy_task.0.parent.action}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{energy_task.1|add:"1"}},severity=2,description='Please run {{energy_task.0.action}} on the coordinates from {{energy_task.parent.action}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
            if et.useqmmm != '{{energy_task.0.useqmmm}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{energy_task.1|add:"1"}},severity=2,description='You {%ifequal energy_task.0.useqmmm 'n'%}used{%else%}did not use{%endifequal%} QM/MM methods in your energy calculation.')
                lessonprob.save()
                return False
            if et.usepbc != '{{energy_task.0.usepbc}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{energy_task.1|add:"1"}},severity=2,description='You {%ifequal energy_task.0.useqmmm 'n'%}used{%else%}did not use{%endifequal%} Periodic Boundary Conditions (PBC) in your energy calculation.')
                lessonprob.save()
                return False
            self.curStep = '{{energy_task.2}}'
            self.save()
            return True
    {%endfor%}blankme
    {%endif%}blankme
        return True

    def onEnergyDone(self,et):
    {%if task_dict.energy_tasks %}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        if et.status == 'F':
            {%for task in task_dict.energy_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=9,description='The job did not complete correctly.')
                lessonprob.save()
                self.curStep = '{{task.1}}'
                self.save()
                return False
            {%endfor%}blankme
        else:
            {%for task in task_dict.energy_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                if float(et.finale) < (float({{energy_task.0.finale}}) - 5.0) or float(et.finale) > (float({{energy_task.0.finale}}) + 5.0):
                    lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{energy_task.1|add:"1"}},severity=2,description='Your calculated energy did not match. Make sure you are using the correct parameter and topology files for this lesson.')
                    lessonprob.save()
                    self.curStep = '{{task.1|add:"1"}}'
                    self.save()
                    return False
                self.curStep = '{{task.1|add:"1"}}'
                self.save()
                return True
            {%endfor%}blankme
        {%endif%}blankme
        return True

    def onNMASubmit(self,nmt):
        {%if task_dict.nmode_tasks%}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        #there's a couple things we can verify here
    {%for nmode_task in task_dict.nmode_tasks%}blankme
        if float(self.curStep) == float({{nmode_task.1|add:"1"}}):
            {%if nmode_task.0.parent%}blankme
            parent_task_action = nmt.parent.action
            if parent_task_action != '{{nmode_task.0.parent.action}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{nmode_task.1|add:"1"}},severity=2,description='Please run {{nmode_task.0.action}} on the coordinates from {{nmode_task.parent.action}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
            if nmt.useqmmm != '{{nmode_task.0.useqmmm}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{nmode_task.1|add:"1"}},severity=2,description='You {%ifequal nmode_task.0.useqmmm 'n'%}used{%else%}did not use{%endifequal%} QM/MM methods in your nmode calculation.')
                lessonprob.save()
                return False
            {% if nmode_task.0.make_nma_movie %}blankme
            {%comment%}blankme
            if you made a movie and you didn't need to, that's fine
            if you DIDN'T make a movie and you needed to, we make a lessonprob
            so, if make_nma_movie was TRue, and you set your task to false, raise
            exception
            {%endcomment%}blankme
            if not nmt.make_nma_movie:
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{nmode_task.1|add:"1"}},severity=2,description='You did not make a movie of your normal mode calculation.')
                lessonprob.save()
                return False
            {%endif%}blankme
            if nmt.useenm != '{{nmode_task.0.useenm}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{nmode_task.1|add:"1"}},severity=2,description='You {%ifequal nmode_task.0.useenm 'n'%}used{%else%}did not use{%endifequal%} an Elastic Network Model (ENM) in your nmode calculation.')
                lessonprob.save()
                return False
            if nmt.nmodes != {{nmode_task.0.nmodes}}:
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{nmode_task.1|add:"1"}},severity=2,description='Number of normal modes was not set to {{nmode_task.0.nmodes}}.')
                lessonprob.save()
                return False
            if nmt.num_trjs != {{nmode_task.0.num_trjs}}:
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{nmode_task.1|add:"1"}},severity=2,description='Number of trajectories to be generated was not set to {{nmode_task.0.num_trjs}}.')
                lessonprob.save()
                return False
            {% if nmode_task.0.useenm == "y" %}blankme
            if float(nmt.rcut) != float({{nmode_task.0.rcut}}):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{nmode_task.1|add:"1"}},severity=2,description='RCUT was not set to {{nmode_task.0.rcut}}')
                lessonprob.save()
                return False
            if float(nmt.kshort) != float({{nmode_task.0.kshort}}):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{nmode_task.1|add:"1"}},severity=2,description='KSHORT was not set to {{nmode_task.0.kshort}}.')
                lessonprob.save()
                return False
            if float(nmt.klong) != float({{nmode_task.0.klong}}):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{nmode_task.1|add:"1"}},severity=2,description='KLONG was not set to {{nmode_task.0.klong}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
            {%comment%}blankme
            vibran doesn't need extra params so we stop here, finally
            {%endcomment%}blankme
            self.curStep = '{{nmode_task.2}}'
            self.save()
            return True
    {%endfor%}blankme
    {%endif%}blankme
        return True

    def onNMADone(self,nmt):
    {%if task_dict.nmode_tasks %}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        if nmt.status == 'F':
            {%for task in task_dict.nmode_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=9,description='The job did not complete correctly.')
                lessonprob.save()
                self.curStep = '{{task.1}}'
                self.save()
                return False
            {%endfor%}blankme
        else:
            {%for task in task_dict.nmode_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                self.curStep = '{{task.1|add:"1"}}'
                self.save()
                return True
            {%endfor%}blankme
        {%endif%}blankme
        return True

    def onMDSubmit(self,mdt):
        {%if task_dict.md_tasks%}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        #there's a couple things we can verify here
    {%for task in task_dict.md_tasks%}blankme
        if float(self.curStep) == float({{task.1|add:"1"}}):
            {%if task.0.parent %}blankme
            parent_task_action = mdt.parent.action
            if parent_task_action != '{{task.0.parent.action}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='Please run {{task.0.action}} on the coordinates from {{task.parent.action}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
{%include 'new_lesson_maker/lessonprob_md.py'%}
            if mdt.ensemble != '{{task.0.ensemble}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='{{task.0.ensemble}} dynamics were not selected.')
                lessonprob.save()
                return False
            {%if task.0.ensemble == "heat"%}blankme
            if mdt.firstt != '{{task.0.firstt}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='FIRSTT was not set to {{task.0.firstt}}.')
                lessonprob.save()
                return False
            if mdt.finalt != '{{task.0.finalt}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='FINALT was not set to {{task.0.finalt}}.')
                lessonprob.save()
                return False
            if mdt.teminc != '{{task.0.teminc}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='teminc was not set to {{task.0.teminc}}.')
                lessonprob.save()
                return False
            if mdt.ihtfrq != '{{task.0.ihtfrq}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='ihtfrq was not set to {{task.0.ihtfrq}}.')
                lessonprob.save()
                return False
            if mdt.tbath != '{{task.0.tbath}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='tbath was not set to {{task.0.tbath}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
            {%comment%}blankme
            ensemble can only ever be one of three or so values
            so we don't have to do elif-style structures, thankfully
            nve has no new params for some reason?
            we need to add these params but we currently don't care a lot
            as much as we care for this module to work.
            {%endcomment%}blankme
            {%if task.0.ensemble == "nve" %}
            {%comment%}blankme
            nve apparently takes nothing? I'm keeping this here because 
            something feels wrong about that
            {%endcomment%}blankme
            {%endif%}blankme
            {%if task.0.ensemble == "nvt" %}blankme
            if mdt.temp != '{{task.0.temp}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='temp was not set to {{task.0.temp}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
            {%if task.0.ensemble == "npt" %}blankme
            if mdt.temp != '{{task.0.temp}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='temp was not set to {{task.0.temp}}.')
                lessonprob.save()
                return False
           {%endif%}blankme
           self.curStep = '{{task.2}}'
            self.save()
            return True
    {%endfor%}blankme
    {%endif%}blankme
        return True

    def onMDDone(self,mdt):
    {%if task_dict.md_tasks %}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        if mdt.status == 'F':
            {%for task in task_dict.md_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=9,description='The job did not complete correctly.')
                lessonprob.save()
                self.curStep = '{{task.1}}'
                self.save()
                return False
            {%endfor%}blankme
        else:
            {%for task in task_dict.md_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                self.curStep = '{{task.1|add:"1"}}'
                self.save()
                return True
            {%endfor%}blankme
        {%endif%}blankme
        return True

    def onLDSubmit(self,mdt):
        {%if task_dict.ld_tasks%}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        #there's a couple things we can verify here
    {%for task in task_dict.ld_tasks%}blankme
        if float(self.curStep) == float({{task.1|add:"1"}}):
            {%if task.0.parent%}blankme
            parent_task_action = mdt.parent.action
            if parent_task_action != '{{task.0.parent.action}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='Please run {{task.0.action}} on the coordinates from {{task.parent.action}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
{%include 'new_lesson_maker/lessonprob_md.py'%}
            if float(mdt.fbeta) != float({{task.0.fbeta}}):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='FBETA was not set to {{task.0.fbeta}}')
                lessonprob.save()
                return False
            self.curStep = '{{task.2}}'
            self.save()
            return True
    {%endfor%}blankme
    {%endif%}blankme
        return True

    def onLDDone(self,ldt):
    {%if task_dict.ld_tasks %}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        if ldt.status == 'F':
            {%for task in task_dict.ld_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=9,description='The job did not complete correctly.')
                lessonprob.save()
                self.curStep = '{{task.1}}'
                self.save()
                return False
            {%endfor%}blankme
        else:
            {%for task in task_dict.ld_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                self.curStep = '{{task.1|add:"1"}}'
                self.save()
                return True
            {%endfor%}blankme
            {%endif%}blankme
        return True

    def onSGLDSubmit(self,mdt):
        {%if task_dict.sgld_tasks%}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        #there's a couple things we can verify here
    {%for task in task_dict.sgld_tasks%}blankme
        if float(self.curStep) == float({{task.1|add:"1"}}):
            {%if task.0.parent%}blankme
            parent_task_action = mdt.parent.action
            if parent_task_action != '{{task.0.parent.action}}':
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='Please run {{task.0.action}} on the coordinates from {{task.parent.action}}.')
                lessonprob.save()
                return False
            {%endif%}blankme
{%include 'new_lesson_maker/lessonprob_md.py'%}
            if float(mdt.tsgavg) != float({{task.0.tsgavg}}):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='TSGAVG was not set to {{task.0.tsgavg}}')
                lessonprob.save()
                return False
            if float(mdt.tempsg) != float({{task.0.tempsg}}):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=2,description='TEMPSG was not set to {{task.0.tempsg}}')
                lessonprob.save()
                return False
            self.curStep = '{{task.2}}'
            self.save()
            return True
    {%endfor%}blankme
    {%endif%}blankme
        return True

    def onSGLDDone(self,sgldt):
    {%if task_dict.sgld_tasks %}blankme
{%include 'new_lesson_maker/lessonprob_tpl.py' %}
        if sgldt.status == 'F':
            {%for task in task_dict.sgld_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                lessonprob = LessonProblem(lesson_type='lesson{{lesson.name}}',lesson_id=self.id,errorstep={{task.1|add:"1"}},severity=9,description='The job did not complete correctly.')
                lessonprob.save()
                self.curStep = '{{task.1}}'
                self.save()
                return False
            {%endfor%}blankme
        else:
            {%for task in task_dict.sgld_tasks %}blankme
            if float(self.curStep) == float('{{task.2}}'):
                self.curStep = '{{task.1|add:"1"}}'
                self.save()
                return True
            {%endfor%}blankme
            {%endif%}blankme
        return True

    #An integer means the job is complete , for example 1 means step 1 was complete, 2 means step 2 was complete
    #A 0.5 means step 1 is running, a 1.5 means step 2 is running
    #If there is a lesson problem that has the same step number as the current step, that means there was a
    #user error in following the lessons and so the lesson step has failed
    #The below code can be used as an example
    def generateStatusHtml(self,file):
        # a template dictionary passing information needed by the template 'generateStatus_Html'
        step_status_list = []
        {%for step in steps %}blankme
        step_status_list.append("<tr class='status'><td class='status'>"+str(int({{step.step}}))+" {%if step.type == "fileupload"%}File Uploaded"{%elif step.type == "buildstruct"%}Working Structure Built"{%else%}"+"{{step.type}}".capitalize(){%endif%}+": ")
        {%endfor%}blankme
        template_dict = {}
        template_dict['status_list'] = []
	#Check to see if a lessonproblem exists, and if so store it into lessonprob
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson{{lesson.name}}',lesson_id=self.id)[0]
        except:
            lessonprob = None
        for i in range(self.nSteps):
            if lessonprob and lessonprob.errorstep == math.floor(self.curStep+1) and math.floor(self.curStep) == i:
                step_status_list[i] += ("<a class='failed' onclick='javascript:open_failure();'>Failed</a></td></tr>")
                continue
            elif (float(self.curStep)-0.5) == i and float(self.curStep) % 1 == 0.5:
                step_status_list[i] += ("<a class='running'>Running</a></td></tr>")
                continue
            elif i < float(self.curStep):
                step_status_list[i] += ("<a class='done'>Done</a></td></tr>")
                continue
            elif i >= float(self.curStep):
                step_status_list[i] += ("<a class='inactive'>N/A</a></td></tr>")
                continue
        return step_status_list
    #This is used to display the lesson page. It tells the template to display which set of directions
    #for the next step
    #Each step is represented as an index in an array. The corresponding value tells the template
    #what is currently happening in regards to the lesson.
    #If the corresponding value is 2, the job is running
    #If it is 1 then the job is finished
    #0 means the job has not stared
    #-1 means there was an error

    def getHtmlStepList(self):
        htmlcode_list = []
        for step in range(self.nSteps):
            htmlcode_list.append(0)
        {%for step in steps_zero%}blankme
        {%comment%}blankme
        steps_zero is just the steps, zero-indexed. A simple but efficient way of getting around django restrictions.
        {%endcomment%}blankme
        {%if forloop.first%}blankme
        if self.curStep > 0:
            htmlcode_list[0] = 1
        {%else%}blankme
        if self.curStep > {{step}}:
            {%if not forloop.first%}blankme
            if self.curStep == {{step}}.5:
                htmlcode_list[{{step}}] = 2
            else:
                htmlcode_list[{{step}}] = 1
            {%else%}blankme
            htmlcode_list[step] = 1
        {%endif%}blankme
            {%endif%}blankme
        {%endfor%}blankme
        try:
            lessonprob = LessonProblem.objects.filter(lesson_type='lesson{{lesson.name}}',lesson_id=self.id)[0]
            htmlcode_list[lessonprob.errorstep-1] = -1
        except:
            lessonprob = None
        return htmlcode_list
