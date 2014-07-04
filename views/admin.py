from django.views import generic
from django.http import HttpResponse
from charmming.models import Structure
from views.mixins import AdminMixin

class AdminIndexView(AdminMixin, generic.TemplateView):
  template_name = 'html/admin/index.html'

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    return self.render_to_response(context)

#  def post(self, request, *args, **kwargs):
#    # update a server setting
