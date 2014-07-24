from django.views import generic
from django.http import HttpResponse
from charmming.models import ProgramSet, InputScriptTemplate
from views.mixins import PermissionsMixin
from views.helpers import most_recent_sessions
import json

class AdminTaskNewView(PermissionsMixin, generic.TemplateView):
  template_name = 'admin/task/new.html'
  permissions = 2

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    context['program_sets'] = ProgramSet.objects.all()[:5]
    return self.render_to_response(context)

  def post(self, request, *args, **kwargs):
    # Create a new program set
    return HttpResponse(json.dumps({"success": True}), content_type='application/json', status=200)
