from django.views import generic
from django.http import HttpResponseRedirect
from charmming.models import User, Session
import datetime

class PermissionsMixin(object):
  def __init__(self):
    self.permissions = 0

  def get_current_user_session(self, request):
    if request.COOKIES.get('session_key'):
      session_key = request.COOKIES.get('session_key')
      current_user_session = Session.objects.filter(key=session_key)
      if current_user_session.count() == 0:
        return None
      return current_user_session[0]

  def get_current_user(self, request):
    current_user_session = self.get_current_user_session(request)
    if not current_user_session:
      return None
    self.current_user_session = current_user_session
    current_user = User.objects.filter(id=current_user_session.user.id)[0]
    return current_user

  def get_context_data(self, **kwargs):
    context = super(PermissionsMixin, self).get_context_data(**kwargs)
    current_user = self.get_current_user(self.request)
    if current_user:
      context['current_user'] = current_user
    return context

  def dispatch(self, request, *args, **kwargs):
    current_user = self.get_current_user(request)
    if not current_user:
      return HttpResponseRedirect('/login')
    else:
      # TODO: Move this into Redis as to not update the DB every page load
      self.current_user_session.last_accessed = datetime.datetime.now()
      self.current_user_session.save()

      # TODO: This will be moved into a group-style permissions set after the initial recoding
      # Once that's done, self.permissions will accept a GroupSet name
      # Currently accepts permissions in the form of an int or, if not given, defaults to 0 (logged in)
      if current_user.permissions < self.permissions:
        return HttpResponseRedirect('/')
      else:
        return super(PermissionsMixin, self).dispatch(request, *args, **kwargs)
