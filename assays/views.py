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
from assays.forms import SearchForm, AssayForm
from assays.rest import get_sd_file

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
                  
  if request.method == 'POST':
    form = AssayForm(request.POST)
    if 'back' in request.POST:
      return HttpResponseRedirect(reverse('search', args=()))
    elif 'download' in request.POST:
      aids_packed = request.POST['current']
      choices = [x.split(":",1) for x in aids_packed.split("|")]
      step = request.POST['step']
      total = request.POST['total']
      form = AssayForm(request.POST,choices=choices,step=step,total=total)
      work_dir = get_dir(request)
      (fh,filename) = mkstemp(dir=work_dir,prefix="pubchem",suffix=".sdf")      
      os.close(fh)
      remove_inconclusive = False
      if 'remove' in request.POST:
        remove_inconclusive = True
      aids = request.POST.getlist('assays')
      num_records = get_sd_file(aids,filename,remove_inconclusive)
      if not num_records:
        err_message = "Empty file returned"
        return render_to_response('assays/assays.html', {'form': form, 'message' : err_message}, context_instance=RequestContext(request)  )
      new_name = "pubchem-assays.sdf"
      mimetype,encoding = mimetypes.guess_type(filename)
      response = HttpResponse(mimetype=mimetype)
      response['Content-Disposition'] = 'attachment; filename=%s' % new_name
      response.write(file(filename, "rb").read())
      return response
  form = AssayForm(request.POST,query=query)
  return render_to_response('assays/assays.html', {'form': form}, context_instance=RequestContext(request)  )

# paging
# insert into charmming/qsar

