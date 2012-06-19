

from django.db import models
from django.contrib.auth.models import User


class Lesson5(models.Model):
    user = models.ForeignKey(User)
    nSteps = models.PositiveIntegerField(default=4)
    curStep = models.DecimalField(default=0,decimal_places=1,max_digits=3)
