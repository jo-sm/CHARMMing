# -*- coding: utf-8 -*-
from django import forms
from django.conf import settings
from rdkit import Chem
from rdkit.Chem import AllChem
from qsar.train import get_sd_properties
from qsar.models import model_types


class TrainingUploadForm(forms.Form):
    
    model_choices = [(obj.id, obj.model_type_name) for obj in model_types.objects.all()]
    model_type=forms.ChoiceField(label="Select Model Type",choices=model_choices)

    model_name=forms.CharField(max_length=100,required=True)
    trainfile = forms.FileField(
        label='Select a training file',
        help_text='SDF with structures and activities'
    )

        
class SelectProperty(forms.Form):
  def __init__(self, *args, **kwargs):
    name = kwargs['filename']
    model_id=kwargs['qsar_model_id']
    fullfilename = kwargs['fullfilename']
    del kwargs['filename']
    del kwargs['qsar_model_id']
    del kwargs['fullfilename']
    query = ""
    if 'query' in kwargs:
     query = kwargs['query']
     del kwargs['query']
    num_mols = 0
    if 'num_mols' in kwargs:
     num_mols = kwargs['num_mols']
     del kwargs['num_mols']
    super(SelectProperty, self).__init__(*args, **kwargs)
    filewidget = forms.HiddenInput(attrs={'value' : name})
    modelwidget = forms.HiddenInput(attrs={'value' : model_id})
    self.fields['filename'] = forms.CharField(widget=filewidget, required=True)
    self.fields['qsar_model_id'] = forms.CharField(widget=modelwidget, required=True)
    self.fields['activity_property'] = forms.ChoiceField(choices=get_sd_properties(fullfilename), label='Select property which contains activity')
    query_widget = forms.HiddenInput(attrs={'value' : query})
    self.fields['query'] = forms.CharField(widget=query_widget, required = True);
    num_mols_widget = forms.HiddenInput(attrs={'value' : num_mols})   
    self.fields['num_mols'] = forms.CharField(widget=num_mols_widget, required = True);

class TestUploadForm(forms.Form):
  def __init__(self, *args, **kwargs):
    active = kwargs['active']
    inactive = kwargs['inactive']
    saved_model = kwargs['saved_model']
    threshold = kwargs['threshold']
    activity_property = kwargs['activity_property']
    filename = kwargs['filename']
    del kwargs['active']
    del kwargs['inactive']
    del kwargs['saved_model']
    del kwargs['threshold']
    del kwargs['activity_property']
    del kwargs['filename']
    super(TestUploadForm, self).__init__(*args, **kwargs)
    activewidget = forms.HiddenInput(attrs={'value' : active})
    inactivewidget = forms.HiddenInput(attrs={'value' : inactive})
    saved_modelwidget = forms.HiddenInput(attrs={'value' : saved_model})
    thresholdwidget = forms.HiddenInput(attrs={'value' : threshold})
    activity_propertywidget = forms.HiddenInput(attrs={'value' : activity_property})
    filenamewidget = forms.HiddenInput(attrs={'value' : filename})
    self.fields['active'] = forms.CharField(widget=activewidget, required=True)
    self.fields['inactive'] = forms.CharField(widget=inactivewidget, required=True)
    self.fields['saved_model'] = forms.CharField(widget=saved_modelwidget, required=True)
    self.fields['threshold'] = forms.CharField(widget=thresholdwidget, required=True)
    self.fields['activity_property'] = forms.CharField(widget=activity_propertywidget, required=True)
    self.fields['filename'] = forms.CharField(widget=filenamewidget, required=True)
    self.fields['predictfile'] = forms.FileField(label='Select a file',help_text='SDF with structures')
    
