from django.views import generic
from django.http import HttpResponse
from charmming.models import Session, ProgramSet, Program, Task
from views.mixins import PermissionsMixin
from django.db.models import Count, Max

class AdminIndexView(PermissionsMixin, generic.TemplateView):
  template_name = 'admin/index.html'
  permissions = 2

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    # This doesn't account for a user who logged out
    context['recent_users'] = Session.objects.values('user__id', 'user__username').annotate(last_accessed=Max('last_accessed'), session_count=Count('id'))[:5]
    context['program_sets'] = ProgramSet.objects.values('name', 'description').annotate(programs_count=Count('program')).order_by('id')[:5]
    context['programs'] = Program.objects.values('name', 'description', 'program_set__name', 'enabled')[:5]
    context['tasks'] = Task.objects.order_by('id')[:5]
    return self.render_to_response(context)

#  def post(self, request, *args, **kwargs):
#    # update a server setting
