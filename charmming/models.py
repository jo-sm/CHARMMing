# This is there the database models should reside

from django.db import models

# A user of CHARMMing
class User(models.Model):
  username = models.CharField(max_length=20)
  password_digest = models.CharField(max_length=200)
  email_address = models.CharField(max_length=65)
  permissions = models.IntegerField()

  @property
  def is_admin(self):
    return self.permissions == 2

# A specific group of users, used for permissions and allowing only certain groups to have access to certain programs (for beta testing, e.g)
class Group(models.Model):
  name = models.CharField(max_length=45)
  description = models.CharField(max_length=150)

# A user session
class Session(models.Model):
  user = models.ForeignKey(User) # user_id
  key = models.CharField(max_length=64)
  user_agent = models.CharField(max_length=100)
  ip_address = models.CharField(max_length=32)
  last_accessed = models.DateTimeField()

# A user's uploaded structure
class Structure(models.Model):
  user = models.ForeignKey(User)
  directory = models.CharField(max_length=25) # directory within user's CHARMMing directory
  file_name = models.CharField(max_length=25) # filename within the directory of the uploaded structure
  name = models.CharField(max_length=25) # display name of structure
  description = models.CharField(max_length=100) # description of structure
  uploaded_at = models.DateTimeField()
  type = models.CharField(max_length=15) # uploaded PDB, uploaded PSF, RCSB downloaded PDB, etc.

class WorkingStructure(models.Model):
  structure_id = models.ForeignKey(Structure)

# A Mustache template of an input script for a given task in a program set
class InputScriptTemplate(models.Model):
  template = models.TextField()

# A set of programs that belong to the same group (CHARMM single vs multi core vs MPI)
class ProgramSet(models.Model):
  name = models.CharField(max_length=30)
  description = models.CharField(max_length=150)
  program_set_template = models.ForeignKey(InputScriptTemplate)

# A single program within a ProgramSet
class Program(models.Model):
  name = models.CharField(max_length=30)
  description = models.CharField(max_length=150)
  program_set = models.ForeignKey(ProgramSet)
  path = models.CharField(max_length=255)
  enabled = models.BooleanField() # for group
  
  # These are always going to be present regardless of the program
  # Therefore, instead of making them generic in the ProgramParameter
  # table below, they are defined on the Program model directly
  input_flag = models.CharField(max_length=10)
  output_flag = models.CharField(max_length=10)

# Parameters that a program can accept
class ProgramParameter(models.Model):
  program = models.ForeignKey(Program)
  name = models.CharField(max_length=25)
  flag = models.CharField(max_length=10)
  use_equals_sign = models.BooleanField() # Should the parameter be separated by a space or an equals sign for the parameter value

# A task for a specific program set, with a given input script template
class Task(models.Model):
  program_set = models.ForeignKey(ProgramSet)
  name = models.CharField(max_length=45)
  input_script_template = models.ForeignKey(InputScriptTemplate)
  type = models.CharField(max_length=25) # Calculation, Analysis? 
  position = models.IntegerField() # Order in the menu
  visible = models.BooleanField(default=True) # is it visible in the menu? Useful for disabling failing tasks and also for tasks such as Build that aren't called by a user explicitly

  # It's possible that a task can have many steps that have to happen in a row, so
  # this is the flag that checks for this and, if true, goes through each of the steps
  # and gets the information
  # Should this be used, and defined a multi step task table, or should multiple tasks be used?
  is_multi_step = models.BooleanField(default=False)

# A parameter for a task, used for both the input script and generic form creation
class TaskParameter(models.Model):
  task = models.ForeignKey(Task)
  name = models.CharField(max_length=30) # variable name for the input script template
  display_name = models.CharField(max_length=30) # name that will be displayed to the user when submitting a task
  description = models.CharField(max_length=150) # longer description of the parameter's use
  required = models.BooleanField()
  type = models.CharField(max_length=15) # Text, radio, checkbox

# If a field has more than one option, such as a checkbox or radio
class TaskParameterOption(models.Model):
  task_parameters = models.ForeignKey(TaskParameter)
  label = models.CharField(max_length=25) # display label for the option
  value = models.CharField(max_length=25) # value of the option


# A job that has been created by a user
class Job(models.Model):
  name = models.CharField(max_length=45)

# Task parameters that are inputted by a user for a job
class JobTaskParameter(models.Model):
  job = models.ForeignKey(Job)
  task_parameter = models.ForeignKey(TaskParameter)
  value = models.CharField(max_length=25) # longer

# A defined lesson
class Lesson(models.Model):
  name = models.CharField(max_length=45)
  description = models.CharField(max_length=150)

# A single step within a lesson
class LessonStep(models.Model):
  lesson = models.ForeignKey(Lesson)
