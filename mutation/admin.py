from django.contrib import admin
from mutation.models import mutateTask

class mutateTaskAdmin(admin.ModelAdmin):
    pass

admin.site.register(mutateTask, mutateTaskAdmin)
