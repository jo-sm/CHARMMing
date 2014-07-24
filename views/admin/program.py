from django.views import generic
from django.http import HttpResponse
from charmming.models import ProgramSet, Program
from views.mixins import PermissionsMixin
import json
import os

class AdminProgramNewView(PermissionsMixin, generic.TemplateView):
  template_name = 'admin/program/new.html'
  permissions = 2

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    context['program_sets'] = ProgramSet.objects.all()[:5]
    return self.render_to_response(context)

  def post(self, request, *args, **kwargs):
    name = request.POST.get('name')
    description = request.POST.get('description')
    enabled = request.POST.get('enabled')
    input_flag = request.POST.get('input_flag')
    output_flag = request.POST.get('output_flag')
    path = request.POST.get('program_path')
    program_set = request.POST.get('program_set')

    program = Program()
    program.name = name
    program.description = description
    program.enabled = enabled
    program.input_flag = input_flag
    program.output_flag = output_flag
    program.path = path
    program.program_set = ProgramSet.objects.filter(id=program_set)[0]
    program.save()

    return HttpResponse(json.dumps({'success': True}), content_type='application/json', status=200)

class AdminProgramVerifyView(PermissionsMixin, generic.TemplateView):
  permissions = 2

  def post(self, request, *args, **kwargs):
    path = request.POST.get('path')
    if os.path.isfile(path) and os.access(path, os.X_OK):
      return HttpResponse(json.dumps({'success': True}), content_type='application/json', status=200)
    else:
      return HttpResponse(status=422)
