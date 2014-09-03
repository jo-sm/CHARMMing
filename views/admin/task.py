from django.views import generic
from django.http import HttpResponse
from django.db import transaction
from django.template.defaultfilters import slugify
from charmming.models import ProgramSet, InputScriptTemplate, Task, TaskGroup, TaskField, TaskConditional, TaskFieldOption
from views.mixins import PermissionsMixin
from views.helpers import most_recent_sessions
import json
import pystache
import re

class AdminTaskNewView(PermissionsMixin, generic.TemplateView):
  template_name = 'admin/task/new.html'
  permissions = 2

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    context['program_sets'] = ProgramSet.objects.all()[:5]
    return self.render_to_response(context)

  def post(self, request, *args, **kwargs):
    post = request.POST
    data = json.loads(post.get('data'))

    _template = data['template']
    parsed_template = pystache.parse(_template)
    # Get each node in the parsed template
    nodes = parsed_template._parse_tree
    node_vars = []
    for node in nodes:
      if hasattr(node, 'key'):
        # node_vars_dict[node.key] = True # We don't care about the value, just the key
        node_vars.append(node.key)
 
    # Don't commit anything to the database unless everything is okay
    # Django 1.6+ uses transaction.atomic()
    with transaction.commit_on_success():
      template = InputScriptTemplate()
      template.template = _template
      template.save()

      task = Task()
      task.name = data['name']
      task.slug = slugify(data['name'])
      task.description = data['description']
      task.program_set_id = data['program_set']
      task.input_script_template = template
      task.type_id = 1 # TODO: non-1 type
      task.order = 1 
      task.visible = True
      task.save()

      # Go through each group and ensure it's valid
      for index, _group in enumerate(data['groups']):
        if len(_group['fields']) == 0:
          pass
        
        group = TaskGroup()
        group.task = task
        group.name = _group['name']
        group.order = _group['order']
        group.internal_id = _group['id']
        group.save()

        if _group.get('conditionals'):
          for _conditional in _group['conditionals']:
            print(_conditional)
            conditional = TaskConditional()
            conditional.object_type = 'TaskGroup'
            conditional.object_id = group.pk

            field = _conditional['field'].split('-')
            conditional.group = field[0]
            conditional.field = field[1]

            conditional.condition = _conditional['condition']
            # may not have value
            if _conditional.get('value'):
              conditional.value = _conditional['value']
            conditional.thenDo = _conditional['thenDo']
            conditional.save()

        for _field in _group['fields']:
          field = TaskField()
          field.group = group
          field.internal_id = _field['id']
          field.name = _field['name']
          field.display_name = _field['displayName']
          field.description = _field['description']
          field.order = _field['order']
          field.required = _field['required']
          field.type = _field['type']['name']
          field.save()

          if field.type == 'radio' or field.type == 'dropdown':
            for _option in _field['type']['values']:
              option = TaskFieldOption()
              option.field = field
              option.label = _option['name']
              option.value = _option['value']
              option.save()

          for _conditional in _field['conditionals']:
            conditional = TaskConditional()
            conditional.object_type = 'TaskField'
            conditional.object_id = field.pk

            field = _conditional['field'].split('-')
            conditional.group = field[0]
            conditional.field = field[1]
            
            conditional.condition = _conditional['condition']
            if _conditional.get('value'):
              conditional.value = _conditional['value']
            conditional.thenDo = _conditional['thenDo']
            conditional.save()


#    template_params = {}
#    for var in post:
#      # every template parameter has a name
#      if re.match('^name-[0-9]*$', var):
#        _s = re.split('-', var)
#        name = post.get('name-' + _s[1])
#        if template_params.get(name):
#          return HttpResponse(json.dumps({"error": "You had two parameters with the name " + name + ". Please rename one of the parameters."}), content_type='application/json', status=422)
#        template_params[name] = {
#          "name": name,
#          "description": post.get('description-' + _s[1]),
#          "type": {
#            "type": post.get('type-' + _s[1])
#          },
#          "conditional": post.get('conditional-' + _s[1]),
#          "order": post.get('order-' + _s[1]),
#          "required": post.get('required-' + _s[1]) == True
#        }
#        if template_params[name]["type"]["type"] == 'radio':
#          template_params[name]["type"]["values"] = []
#          # Iterate over each variable in the post and get the type-n-radio-type-n code
#          for v in post:
#            if re.match('^type-' + _s[1] + '-radio-name-[0-9]*$', v):
#              type_split = re.split('-', v)
#              template_params[name]["type"]["values"].append({
#                "name": post.get(v),
#                "value": post.get('type-' + _s[1] + '-radio-value-' + type_split[len(type_split)-1])
#              })
#        elif template_params[name]["type"]["type"] == 'checkbox':
#          template_params[name]["type"]["checked"] = post.get('type-' + _s[1] + '-checkbox-checked')
#    
#    invalid_params = []
#    for var in template_params:
#      if not node_vars_dict.get(var):
#        invalid_params.append(var)
#    if invalid_params:
#      _p = 'parameters' if len(invalid_params) > 1 else 'parameter'
#      _w = 'were' if len(invalid_params) > 1 else 'was'
#      _s = 's' if len(invalid_params) > 1 else ''
#      return HttpResponse(json.dumps({"error": "The following " + _p + ' ' + _w + " defined but are not present in the task template: " + ', '.join(invalid_params) + '. Please double check the spelling or the existance of the variable' + _s + ' in the task template, or remove the defined parameter.'}), content_type='application/json', status=422)
    return HttpResponse(json.dumps({"success": True}), content_type='application/json', status=200)
