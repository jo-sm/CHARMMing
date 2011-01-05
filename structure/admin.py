from django.contrib import admin
from structure.models import Structure, energyParams, goModel, blnModel

class structureAdmin(admin.ModelAdmin):
    pass

admin.site.register(Structure, structureAdmin)
admin.site.register(energyParams)
admin.site.register(goModel)
admin.site.register(blnModel)
