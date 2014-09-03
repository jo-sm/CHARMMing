from django.views import generic
from django.http import HttpResponse, HttpResponseRedirect
from charmming.models import Task
from views.mixins import PermissionsMixin
import datetime
import json
import os

class ProgramSetTaskIndexView(PermissionsMixin, generic.TemplateView):
  def get(self, request, *args, **kwargs):
    return HttpResponse(status=200)

class ProgramSetIndexView(PermissionsMixin, generic.TemplateView):
  def get(self, request, *args, **kwargs):
    return HttpResponse(status=200)

class ProgramSetShowView(PermissionsMixin, generic.TemplateView):
  def get(self, request, *args, **kwargs):
    return HttpResponse(status=200)

class ProgramSetTaskShowView(PermissionsMixin, generic.TemplateView):
  template_name = 'program_set/show.html'

  def get(self, request, *args, **kwargs):
    data = {
      'groups': []
    }
    slug = kwargs.get('slug')
    task_slug = kwargs.get('task_slug')
    task = Task.objects.filter(program_set__slug=slug, slug=task_slug)
    
    if not task:
      return HttpResponse(status=404)
 
    task = task[0]
    data['name'] = task.name
    data['description'] = task.description

    for _group in task.taskgroup_set.all():
      group = {}
      group['name'] = _group.name
      group['order'] = _group.order
      group['id'] = _group.internal_id
      group['conditionals'] = []

      for _conditional in _group.conditionals:
        conditional = {
          'condition': _conditional.condition,
          'field': str(_conditional.group) + '-' + str(_conditional.field),
          'value': _conditional.value,
          'thenDo': _conditional.thenDo
        }

        group['conditionals'].append(conditional)

      group['fields'] = []

      for _field in _group.taskfield_set.all():
        field = {}
        field['name'] = _field.name
        field['description'] = _field.description
        field['displayName'] = _field.display_name
        field['id'] = _field.internal_id
        field['order'] = _field.order
        field['required'] = _field.required
        field['conditionals'] = []

        for _conditional in _field.conditionals:
          conditional = {
            'condition': _conditional.condition,
            'field': str(_conditional.group) + '-' + str(_conditional.field),
            'value': _conditional.value,
            'thenDo': _conditional.thenDo
          }
          field['conditionals'].append(conditional)

        field['type'] = {}
        field['type']['name'] = _field.type

        if _field.type == 'radio' or _field.type == 'dropdown':
          field['type']['values'] = []
          for _option in _field.taskfieldoption_set.all():
            field['type']['values'].append({
              'name': _option.label,
              'value': _option.value
            })

        group['fields'].append(field)

      data['groups'].append(group)

    context = self.get_context_data()
    context['data'] = json.dumps(data)
    context['task'] = task

    return self.render_to_response(context)
