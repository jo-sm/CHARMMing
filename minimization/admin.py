from django.contrib import admin
from minimization.models import minimizeParams

class minimizeParamsAdmin(admin.ModelAdmin):
    pass

admin.site.register(minimizeParams, minimizeParamsAdmin)
