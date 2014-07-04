from django.views import generic
from django.http import HttpResponseRedirect
from charmming.models import User, Session

class LoggedInMixin(object):
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
    current_user = User.objects.filter(id=current_user_session.user.id)[0]
    return current_user

  def get_context_data(self, **kwargs):
      context = super(LoggedInMixin, self).get_context_data(**kwargs)
      current_user = self.get_current_user(self.request)
      if current_user:
        context['current_user'] = current_user
      return context

  def dispatch(self, request, *args, **kwargs):
    current_user_session = self.get_current_user_session(request)
    if not current_user_session:
      return HttpResponseRedirect('/login')
    else:
      return super(LoggedInMixin, self).dispatch(request, *args, **kwargs)

class AdminMixin(LoggedInMixin):
  def dispatch(self, request, *args, **kwargs):
    #session_key = request.COOKIES.get('session_key')
    #current_user_session = Session.objects.filter(key=session_key)[0]
    #current_user = User.objects.filter(id=current_user_session.user.id)[0]
    current_user = self.get_current_user(request)
    if not current_user.is_admin:
      return HttpResponseRedirect('/') # add descriptive error message
    else:
      return super(AdminMixin, self).dispatch(request, *args, **kwargs)
