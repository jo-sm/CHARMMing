from django.contrib import admin
from pdbinfo.models import PDBFile, energyParams, goModel, blnModel

class pdbfileAdmin(admin.ModelAdmin):
    pass

admin.site.register(PDBFile, pdbfileAdmin)
admin.site.register(energyParams)
admin.site.register(goModel)
admin.site.register(blnModel)
