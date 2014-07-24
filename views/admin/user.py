import bcrypt
import json
from django.views import generic
from django.http import HttpResponse
from charmming.models import User
from views.mixins import PermissionsMixin

class AdminUserIndexView(PermissionsMixin, generic.TemplateView):
  template_name = 'admin/user/index.html'
  permissions = 2

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    context['users'] = User.objects.all()
    return self.render_to_response(context)

class AdminUserShowView(PermissionsMixin, generic.TemplateView):
  template_name = 'admin/user/show.html'
  permissions = 2

class AdminUserNewView(PermissionsMixin, generic.TemplateView):
  template_name = 'admin/user/new.html'
  permissions = 2

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    return self.render_to_response(context)

  def post(self, request, *args, **kwargs):
    username = request.POST.get('username')
    password = request.POST.get('password')
    permissions = request.POST.get('permissions')
    email_address = request.POST.get('email_address')
    dni = []
    if not username:
      dni.append('username')
    if not password:
      dni.append('password')
    if not permissions:
      dni.append('permissions')
    if not email_address:
      dni.append('email address')
    if len(dni) > 0:
      return HttpResponse(json.dumps({'error': 'You did not enter the following fields: ' + ', '.join(dni) + '. Please resubmit with these fields inputted.'}), content_type='application/json', status=422)
    new_user = User()
    new_user.username = username
    new_user.email_address = email_address
    new_user.password_digest = bcrypt.hashpw(password.encode('UTF-8'), bcrypt.gensalt(12))
    new_user.permissions = permissions
    new_user.save()
    return HttpResponse(json.dumps({'success': True}), content_type='application/json', status=200)
