from django.views import generic
from django.http import HttpResponse
from django.core import serializers
from django.forms.models import model_to_dict
from charmming.models import ProgramSet, InputScriptTemplate
from views.mixins import PermissionsMixin, APIMixin
from views.helpers import most_recent_sessions
import json
import re

class AdminProgramSetIndexView(PermissionsMixin, APIMixin, generic.TemplateView):
  template_name = 'admin/program_set/index.html'
  permissions = 2

  class JSON:
    def get(self, request, *args, **kwargs):
      context = self.get_context_data()
      program_sets = ProgramSet.objects.all()
      if re.search('application/json', request.META.get('HTTP_ACCEPT')):
        _list = []
        for _s in program_sets:
          _list.append(model_to_dict(_s))
        return HttpResponse(json.dumps(_list), content_type='application/json', status=200)

  def get(self, request, *args, **kwargs):
    return HttpResponse(json.dumps({'hello': 'world'}), content_type='application/json', status=200)

class AdminProgramSetNewView(PermissionsMixin, generic.TemplateView):
  template_name = 'admin/program_set/new.html'
  permissions = 2

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    return self.render_to_response(context)

  def post(self, request, *args, **kwargs):
    # Create a new program set
    name = request.POST.get('name')
    slug = slugify(name)
    description = request.POST.get('description')
    template = request.POST.get('template')

    program_set_template = InputScriptTemplate.objects.create(template=template)
    program_set = ProgramSet()
    program_set.name = name
    program_set.slug = slug
    program_set.description = description
    program_set.program_set_template = program_set_template
    program_set.save()

    return HttpResponse(json.dumps({"success": True}), content_type='application/json', status=200)
