# This is there the database models should reside

from django.db import models

# A user of CHARMMing
class User(models.Model):
  username = models.CharField(max_length=20)
  password_digest = models.CharField(max_length=200)
  permissions = models.IntegerField()
  
  @property
  def is_admin(self):
    return self.permissions == 2

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

  # lesson_maker
  lesson_maker_active = models.BooleanField(default=False)

class WorkingStructure(models.Model):
  structure_id = models.ForeignKey(Structure)
