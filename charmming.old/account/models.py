#
#                            PUBLIC DOMAIN NOTICE
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the authors' official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software is freely available
#  to the public for use.  There is no restriction on its use or
#  reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, NIH and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. NIH, NHLBI, and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any
#  particular purpose.
from django.db import models
from django.contrib.auth.models import User
from django.contrib.localflavor.us.forms import USPhoneNumberField
import django.forms

class teacherProfile(models.Model):
    institute  = models.CharField(max_length=200)
    phone_number = USPhoneNumberField()
    teacher = models.ForeignKey(User, unique=True)

class classroom(models.Model):
    name =  models.CharField(max_length=200)
    teacher = models.ForeignKey(teacherProfile)
    num_of_students = models.IntegerField()
    
    #pre: Needs a valid teacherProfile
    #This will get the list of classrooms based on a teacher
    def getClassList(self, teacher):
        return classroom.objects.filter(teacher=teacher)

class studentProfile(models.Model):
    institute  = models.CharField(max_length=200)
    student = models.ForeignKey(User, unique=True)
    charmm_license = models.BooleanField()

    #Gets the studentProfile based on the request's user
    def getStudentFromUser(request):
        try:
            student_profile = studentProfile.objects.filter(student = request.user)[0]
	    return student_profile
	except:
	    return None
