from django.contrib import admin
from account.models import teacherProfile, studentProfile, classroom

admin.site.register(teacherProfile)
admin.site.register(studentProfile)
admin.site.register(classroom)
