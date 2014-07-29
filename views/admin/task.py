from django.views import generic
from django.http import HttpResponse
from charmming.models import ProgramSet, InputScriptTemplate
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
    name = post.get('name')
    description = post.get('description')
    template = post.get('template')
    parsed_template = pystache.parse(template)
    # Get each node in the parsed template
    nodes = parsed_template._parse_tree
    node_vars_dict = {}
    for node in nodes:
      if hasattr(node, 'key'):
        node_vars_dict[node.key] = True # We don't care about the value, just the key

    template_params = {}
    for var in post:
      # Every template parameter has a name
      if re.match('^name-[0-9]*$', var):
        _s = re.split('-', var)
        name = post.get('name-' + _s[1])
        if template_params.get(name):
          return HttpResponse(json.dumps({"error": "You had two parameters with the name " + name + ". Please rename one of the parameters."}), content_type='application/json', status=422)
        template_params[name] = {
          "name": name,
          "description": post.get('description-' + _s[1]),
          "type": {
            "type": post.get('type-' + _s[1])
          },
          "conditional": post.get('conditional-' + _s[1]),
          "order": post.get('order-' + _s[1]),
          "required": post.get('required-' + _s[1]) == True
        }
        if template_params[name]["type"]["type"] == 'radio':
          template_params[name]["type"]["values"] = []
          # Iterate over each variable in the post and get the type-n-radio-type-n code
          for v in post:
            if re.match('^type-' + _s[1] + '-radio-name-[0-9]*$', v):
              type_split = re.split('-', v)
              template_params[name]["type"]["values"].append({
                "name": post.get(v),
                "value": post.get('type-' + _s[1] + '-radio-value-' + type_split[len(type_split)-1])
              })
        elif template_params[name]["type"]["type"] == 'checkbox':
          template_params[name]["type"]["checked"] = post.get('type-' + _s[1] + '-checkbox-checked')
    
    invalid_params = []
    for var in template_params:
      if not node_vars_dict.get(var):
        invalid_params.append(var)
    if invalid_params:
      _p = 'parameters' if len(invalid_params) > 1 else 'parameter'
      _w = 'were' if len(invalid_params) > 1 else 'was'
      _s = 's' if len(invalid_params) > 1 else ''
      return HttpResponse(json.dumps({"error": "The following " + _p + ' ' + _w + " defined but are not present in the task template: " + ', '.join(invalid_params) + '. Please double check the spelling or the existance of the variable' + _s + ' in the task template, or remove the defined parameter.'}), content_type='application/json', status=422)
    return HttpResponse(json.dumps(template_params), content_type='application/json', status=200)
