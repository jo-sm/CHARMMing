#-*- coding: utf-8 -*-
from itertools import chain

from django import forms
from django.utils.html import conditional_escape
from django.utils.translation import ugettext, ugettext_lazy
from django.utils.encoding import force_text, python_2_unicode_compatible
from django.utils.safestring import mark_safe
from django.conf import settings
from assays.rest import get_list_aids, get_url_base

class SearchForm(forms.Form):
  query = forms.CharField(max_length=200,required=True)

# see /usr/lib/python2.7/dist-packages/django/forms/widgets.py
class CheckboxSelectMultipleWithUrl(forms.CheckboxSelectMultiple):
    def render(self, name, value, attrs=None, choices=()):
        if value is None: value = []
        has_id = attrs and 'id' in attrs
        final_attrs = self.build_attrs(attrs, name=name)
        output = ['<ul style="list-style-type:none">']
        # Normalize to strings
        str_values = set([force_text(v) for v in value])
        for i, (option_value, option_label) in enumerate(chain(self.choices, choices)):
            # If an ID attribute was given, add a numeric index as a suffix,
            # so that the checkboxes don't all have the same ID attribute.
            if has_id:
                final_attrs = dict(final_attrs, id='%s_%s' % (attrs['id'], i))
                label_for = u' for="%s"' % final_attrs['id']
            else:
                label_for = ''

            cb = forms.CheckboxInput(final_attrs, check_test=lambda value: value in str_values)
            option_value = force_text(option_value)
            rendered_cb = cb.render(name, option_value)
            option_label = force_text(option_label)
            output.append(u'<li style="margin-bottom:0.5em"><label%s>%s <a href="%s">%s</a> %s</label></li>' %
                                      (label_for, rendered_cb, get_url_base()+option_value,option_value, option_label))
        output.append('</ul>')
        return mark_safe('\n'.join(output))

class AssayForm(forms.Form):
    def __init__(self, *args, **kwargs):
      aids = []
      step = 0
      total = 0
      if 'choices' in kwargs:
        aids = kwargs['choices']
        del kwargs['choices']
      elif 'query' in kwargs:
        query = kwargs['query']
        del kwargs['query']
        (step,total,aids) = get_list_aids(query)
      if 'step' in kwargs:
        step = kwargs['step']
        del kwargs['step']
      if 'total' in kwargs:
        total = kwargs['total']
        del kwargs['total']
      
      current = ""
      for a in aids:
        if current:
          current += "|"
        current += a[0]+":"+a[1]
          
      super(AssayForm, self).__init__(*args, **kwargs)
      self.fields['assays'] = forms.MultipleChoiceField(choices=aids, label="assays",required=False,
                              widget = CheckboxSelectMultipleWithUrl())
      self.fields['remove'] = forms.BooleanField(label="remove",required=False)
      current_widget = forms.HiddenInput(attrs={'value' : current})
      self.fields['current'] = forms.CharField(widget=current_widget, required = True); 
      step_widget = forms.HiddenInput(attrs={'value' : step})
      self.fields['step'] = forms.CharField(widget=step_widget, required = False);
      total_widget = forms.HiddenInput(attrs={'value' : total})
      self.fields['total'] = forms.CharField(widget=total_widget, required = False);
