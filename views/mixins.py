from django.views import generic
from django.http import HttpResponseRedirect
from charmming.models import User, Session, ProgramSet
import datetime
import re
import logging

class APIMixin(object):

  def dispatch(self, request, *args, **kwargs):
    self.request = request
    self.args = args
    self.kwargs = kwargs

    # If the method type is an accepted API method, such as JSON
    # look for its defined class within the view
    # If the class exists, call that handler, and if it doesn't,
    # default to the original implementation
    
    # In order to allow the API class, such as class JSON, to inherit
    # nothing rather than the generic.TemplateView class, im_func is
    # used to copy the function rather than the use the reference
    if re.search('application/json', request.META.get('HTTP_ACCEPT')):
      # Look for the JSON class and call the method if it exists
      klass = getattr(self, 'JSON', None)
      if klass:
        handler = getattr(klass, request.method.lower(), None)
        if handler:
          setattr(self, "__json_" + handler.__name__, handler.im_func)
          _handler = getattr(self, "__json_" + request.method.lower())
          return _handler(self, request, *args, **kwargs)

    return super(APIMixin, self).dispatch(request, *args, **kwargs)

class BaseMixin(generic.TemplateView):
  def get_context_data(self, **kwargs):
    context = super(BaseMixin, self).get_context_data(**kwargs)
    program_sets = ProgramSet.objects.all()
    tasks = []
    for program_set in program_sets:
      parent = {}
      parent['name'] = program_set.name
      parent['slug'] = program_set.slug
      parent['tasks'] = []
      for task in program_set.task_set.all():
        parent['tasks'].append(task)
      if len(parent['tasks']) > 0:
        tasks.append(parent)
    context['tasks_with_parents'] = tasks
    return context

class PermissionsMixin(BaseMixin):
  def __init__(self):
    self.permissions = 0

  def get_current_user_session(self, request):
    if request.COOKIES.get('session_key'):
      session_key = request.COOKIES.get('session_key')
      current_user_session = Session.objects.filter(key=session_key)
      if current_user_session.count() == 0:
        return None
      return current_user_session[0]
    return None

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
    self.request = request

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
