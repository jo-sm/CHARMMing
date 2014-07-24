from django.views import generic
from django.http import HttpResponse
from charmming.models import ProgramSet, InputScriptTemplate
from views.mixins import PermissionsMixin
from views.helpers import most_recent_sessions
import json

class AdminProgramSetNewView(PermissionsMixin, generic.TemplateView):
  template_name = 'admin/program_set/new.html'
  permissions = 2

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    return self.render_to_response(context)

  def post(self, request, *args, **kwargs):
    # Create a new program set
    name = request.POST.get('name')
    description = request.POST.get('description')
    template = request.POST.get('template')

    program_set_template = InputScriptTemplate.objects.create(template=template)
    program_set = ProgramSet.objects.create(name=name, description=description, program_set_template=program_set_template)
    return HttpResponse(json.dumps({"success": True}), content_type='application/json', status=200)
