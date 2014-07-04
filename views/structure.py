from django.views import generic
from django.http import HttpResponse
from charmming.models import Structure
from views.mixins import LoggedInMixin

class StructureIndexView(LoggedInMixin, generic.TemplateView):
  template_name = 'html/structure/index.html'

  def get(self, request, *args, **kwargs):
    # view all structures
    context = self.get_context_data()
    context['structures'] = Structure.objects.filter(user=context['current_user'])
    return self.render_to_response(context)

#  def post(self, request, *args, **kwargs):
#    # submit a structure

class StructureNewView(LoggedInMixin, generic.TemplateView):
  template_name = 'html/structure/new.html'

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    return self.render_to_response(context)

class StructureShowView(LoggedInMixin, generic.TemplateView):
  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    return render_to_response(context)
    # show a single structure with its working structures, etc.

#  def delete(self, request, *args, **kwargs):
#   # delete a structure and all associated working structures
