# templatetags/navigation.py
# http://blog.scur.pl/2012/09/highlighting-current-active-page-django/
 
from django import template
from django.core import urlresolvers
import logging
 
register = template.Library()
 
@register.simple_tag(takes_context=True)
def is_current_page(context, url_name, **kwargs):
    matches = current_url_equals(context, url_name, **kwargs)
    return 'active' if matches else ''
 
@register.simple_tag(takes_context=True)
def give_parent_menu_class(context, parent_name, **kwargs):
  """
  Determines parent menu class based off of two factors:
  Is the menu expanded?, and, is it active (a sub page or itself is currently being viewed)
  """
  _class = []
  url_name = urlresolvers.resolve(context.get('request').path).url_name
  if url_name.split('.')[0] == parent_name:
    _class.append('active')
    if not url_name == parent_name + ".index":
      _class.append('with-submenu')
  if submenu_expanded(context, parent_name):
    _class.append('expanded') #if submenu_expanded(context, parent_name)
  return ' '.join(_class)

@register.simple_tag(takes_context=True)
def is_submenu_expanded(context, menu_name, **kwargs):
  _expanded = submenu_expanded(context, menu_name)
  return 'ion-minus-round' if _expanded else 'ion-plus-round'

def submenu_expanded(context, menu_name):
  """
  Determines if the submenu is expanded based off of current page and cookie settings
  Current page takes precedence over cookie settings
  """
  if urlresolvers.resolve(context.get('request').path).url_name.split('.')[0] == menu_name:
    return True

  menu_settings = context.get('request').COOKIES.get('menu_settings')
  if menu_settings:
    if menu_settings.get(menu_name) == 'expanded':
      return True
  return False

def current_url_equals(context, url_name, **kwargs):
    resolved = False
    try:
        resolved = urlresolvers.resolve(context.get('request').path)
    except:
        pass
    matches = resolved and resolved.url_name == url_name
    if matches and kwargs:
        for key in kwargs:
            kwarg = kwargs.get(key)
            resolved_kwarg = resolved.kwargs.get(key)
            if not resolved_kwarg or kwarg != resolved_kwarg:
                return False
    return matches
