from django.contrib import admin
from minimization.models import minimizeTask

class minimizeTaskAdmin(admin.ModelAdmin):
    pass

admin.site.register(minimizeTask, minimizeTaskAdmin)
