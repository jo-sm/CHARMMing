from django.views import generic
from django.http import HttpResponse, HttpResponseRedirect

import json
import bcrypt
import datetime, random, string
from views.mixins import PermissionsMixin
from charmming.models import User, Session

class DashboardView(PermissionsMixin):
  template_name = "dashboard.html"

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    return self.render_to_response(context)

class LoginView(generic.TemplateView):
  template_name = "login.html"

  def get(self, request, *args, **kwargs):
    return self.render_to_response([])

  def post(self, request, *args, **kwargs):
    # user logs in
    user = User.objects.filter(username = request.POST.get('username'))
    if user.count() == 0:
      # Invalid username
      error = {"error": "Username and/or password incorrect"}
      return HttpResponse(json.dumps(error), content_type="application/json", status=422)

    user = user[0]
    password = request.POST.get('password').encode('utf-8')
    password_digest = user.password_digest.encode('utf-8')
    if not bcrypt.hashpw(password, password_digest) == password_digest:
      # Invalid password
      error = {"error": "Username and/or password incorrect"}
      return HttpResponse(json.dumps(error), content_type="application/json", status=422)

    # User is authenticated, create a new session
    session_key = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(12))
    session = Session(user = user, key = session_key, user_agent = request.META['HTTP_USER_AGENT'], ip_address = request.META['REMOTE_ADDR'], last_accessed = datetime.datetime.now()).save()

    response = HttpResponse(json.dumps({"success": True}), content_type="application/json")
    response.set_cookie('session_key', session_key) #encrypt_cookie_value
    return response

class LogoutView(PermissionsMixin):
  def get(self, request, *args, **kwargs):
    # Delete the current session from the database
    # We know this is set due to using LoggedInMixin
    session_key = request.COOKIES.get('session_key')
    Session.objects.filter(key = session_key).delete()
    response = HttpResponseRedirect('/login')
    response.delete_cookie('session_key')
    return response
