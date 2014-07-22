# Create your views here.
from django.http import HttpResponseRedirect,HttpResponse
from django.shortcuts import render_to_response
from django.core.urlresolvers import reverse
from django.template import RequestContext
from tempfile import mkstemp
import os, os.path
import subprocess
import mimetypes
import fnmatch
import urllib, urllib2
import unicodedata
import random
import charmming_config
from django.contrib.auth.models import User
from assays.forms import SearchForm, AssayForm, TrainingSubmitForm
from assays.rest import get_sd_file
from qsar.models import qsar_models, model_types

def get_dir(request):
  username = request.user.username
  u = User.objects.get(username=username)
  location = charmming_config.user_home + '/' + request.user.username + '/qsar/'
  return str(location)


def search(request):
  if not request.user.is_authenticated():
    return render_to_response('html/loggedout.html')
  if request.method == 'POST':
    form = SearchForm(request.POST)
    if form.is_valid():
      query = request.POST['query']
      return assays(request,query)
    else:
      return render_to_response('assays/search.html', {'form': form}, context_instance=RequestContext(request)  )
  else:
    form = SearchForm()
  return render_to_response('assays/search.html', {'form': form}, context_instance=RequestContext(request)  )
  
def assays(request,query=None):
  if not request.user.is_authenticated():
    return render_to_response('html/loggedout.html')
  start = 0
  total = 0
  step = 0
  selected = dict()
  if request.method == 'POST':
    form = AssayForm(request.POST)
    aids = []
    if 'start_form' in request.POST: 
        start = int(request.POST['start_form'])
    if 'assays' in request.POST:
        aids = request.POST.getlist('assays')
    if 'selected' in request.POST:
        selected_packed = request.POST['selected']
        for x in selected_packed.split("|"):
            p =  x.split(":",1)
            if len(p) == 2:
                selected[p[0]] = p[1].split(",")
        selected[str(start)] = aids
    if 'step' in request.POST:
        step = int(request.POST['step'])
    if 'query' in request.POST:
        query = request.POST['query']                                               
    if 'total' in request.POST:
        total = request.POST['total']
        
    if 'back' in request.POST:
      return HttpResponseRedirect(reverse('search', args=()))
    elif 'next' in request.POST:
        start += step
    elif 'prev' in request.POST:
        start -= step
    elif 'download' in request.POST or 'continue' in request.POST:
      aids_packed = request.POST['current']
      choices = [x.split(":",1) for x in aids_packed.split("|")]
      form = AssayForm(request.POST,choices=choices,step=step,total=total,selected=selected)
      work_dir = get_dir(request)
      (fh,filename) = mkstemp(dir=work_dir,prefix="pubchem",suffix=".sdf")      
      os.close(fh)
      remove_inconclusive = False
      if 'remove' in request.POST:
        remove_inconclusive = True
      all_aids = []
      for s in selected:
          all_aids += selected[s]
      num_records = get_sd_file(all_aids,filename,remove_inconclusive)
      if not num_records:
        err_message = "Empty file returned"
        return render_to_response('assays/assays.html', {'form': form, 'message' : err_message}, context_instance=RequestContext(request)  )
      if 'download' in request.POST:
        new_name = "pubchem-assays.sdf"
        mimetype,encoding = mimetypes.guess_type(filename)
        response = HttpResponse(mimetype=mimetype)
        response['Content-Disposition'] = 'attachment; filename=%s' % new_name
        response.write(file(filename, "rb").read())
        return response
      else:
        form = TrainingSubmitForm(filename=filename,query=query,num_mols=num_records) 
        return render_to_response('assays/continue.html', {'form': form, 'query': query, 'num_mols': num_records}, context_instance=RequestContext(request)  )
  form = AssayForm(query=query, start = start, selected=selected)
  if form.len_aids > 0 :
    return render_to_response('assays/assays.html', {'form': form}, context_instance=RequestContext(request)  )
  else:
    return HttpResponseRedirect(reverse('search', args=()))

def cont(request,filename="",query="",num_mols=0):
    if not request.user.is_authenticated():
            return render_to_response('html/loggedout.html')
    log=open("/tmp/cont.log", "w")
    err_message=""
    if request.method == 'POST':
        if 'back' in request.POST:
            query = request.POST['query']
            form = AssayForm(request.POST,query=query)
            return render_to_response('assays/assays.html', {'form': form},context_instance=RequestContext(request)  ) 
        form = TrainingSubmitForm(request.POST)
        log.write("model type:%s\n" % (request.POST['model_type']))
        if form.is_valid():
            name = request.POST['filename']
            activity_property_choice_length = 1
            #activity_property_choice_length = len(get_sd_properties(name));
            #if not activity_property_choice_length:
            #    err_message="Not an SD file or no SD properties found. Please upload an SD file with Activity Property"
            #    return render_to_response('qsar/newmodel.html', {'form': form, 'message':err_message}, context_instance=RequestContext(request)  )
            u = User.objects.get(username=request.user.username)
            newmodel=qsar_models()
            newmodel.model_owner=request.user
            newmodel.model_owner_index=newmodel.getNewModelOwnerIndex(u)
            newmodel.model_name=request.POST['model_name']
            newmodel.model_type=model_types.objects.get(id=request.POST['model_type'])
            newmodel.save()
            filename = os.path.basename(name)
            log.write("newmodel:%s" % (newmodel.id))
            query = request.POST['query']
            num_mols = request.POST['num_mols']
            return HttpResponseRedirect(reverse('property', kwargs={'filename':filename,'qsar_model_id':newmodel.id,'query':query,'num_mols' : num_mols}))
    else:
      form = TrainingSubmitForm(filename=filename,query=query,num_mols=num_mols)
      return render_to_response('assays/continue.html', {'form': form, 'query': query, 'num_mols': num_mols}, context_instance=RequestContext(request))
      