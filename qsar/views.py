# -*- coding: utf-8 -*-
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect,HttpResponse
from django.core.urlresolvers import reverse
from qsar.forms import TrainingUploadForm,SelectProperty,TestUploadForm
from tempfile import mkstemp
import os
import subprocess
import mimetypes
import charmming_config
from django.contrib.auth.models import User
from rdkit import Chem
from qsar.train import active_inactive,get_sd_properties
from qsar.models import jobs, job_types
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from account.views import checkPermissions
import MySQLdb
import MySQLdb.cursors
import settings
from qsar import common

def get_dir(request):
    username = request.user.username
    u = User.objects.get(username=username)
    location = charmming_config.user_home + '/' + request.user.username + '/'
    return location

def handle_uploaded_file(f,filename):
  path = filename 
  with open(path, 'wb+') as destination:
    for chunk in f.chunks():
      destination.write(chunk)
                        
def upload(request):
    if not request.user.is_authenticated():
            return render_to_response('html/loggedout.html')
            
    if request.method == 'POST':
        form = TrainingUploadForm(request.POST, request.FILES)
        if form.is_valid():
            (fh,name) = mkstemp(dir=get_dir(request))
            os.close(fh);        
            handle_uploaded_file(request.FILES['trainfile'],name)
            #return HttpResponseRedirect(reverse('property', args=[name]))
            return property(request,name)
    else:
        form = TrainingUploadForm() # A empty, unbound form

    return render_to_response('qsar/upload.html', {'form': form}, context_instance=RequestContext(request)  )

def newModel(request):
    if not request.user.is_authenticated():
            return render_to_response('html/loggedout.html')

    if request.method == 'POST':
        form = TrainingUploadForm(request.POST, request.FILES)
        if form.is_valid():
            (fh,name) = mkstemp(dir=get_dir(request))
            os.close(fh);
            handle_uploaded_file(request.FILES['trainfile'],name)
            #return HttpResponseRedirect(reverse('property', args=[name]))
            u = User.objects.get(username=username)
            newmodel=qsar_model()
            newmodel.model_owner=request.user
            newmodel.model_owner_index=common.getNewModelOwnerIndex(u)
            newmodel.model_name=request.POST['model_name']
            newmodel.model_type=request.POST['model_type']
            newmodel.save()
            return property(request,name,qsar_model)
        else:
            return render_to_response('qsar/newmodel.html', {'form': form}, context_instance=RequestContext(request)  )
    else:
        form = TrainingUploadForm() # A empty, unbound form

    return render_to_response('qsar/newmodel.html', {'form': form}, context_instance=RequestContext(request)  )


def property(request,filename=None,qsar_model=None):
  if not request.user.is_authenticated():
    return render_to_response('html/loggedout.html')
        
  if request.method == 'POST':
    if 'filename' in request.POST:
      filename = request.POST['filename']
    form = SelectProperty(request.POST,filename=filename)
    if 'back' in request.POST:
      return HttpResponseRedirect(reverse('upload', args=()))
    elif 'next' in request.POST:
      if form.is_valid():
        activity_property = request.POST["activity_property"]
        common.AssignObjectAttribute(request,qsar_model.id,"qsar_qsar_models","Activity Property",activity_property)
        model_type=model_types.objects.get(id=qsar_model.model_type_id)
        return train(request,filename,activity_property,qsar_model)
  else:
    form = SelectProperty(filename=filename)
  activity_property_choice_length = 0;
  if filename is not None:
    activity_property_choice_length = len(get_sd_properties(filename));
  return render_to_response('qsar/property.html', {'form': form,'filename':filename,'activity_property_choice_length':activity_property_choice_length, 
                            'qsar_model':qsar_model, 'model_type':model_type}, context_instance=RequestContext(request))

def train(request,filename=None,activity_property=None,qsar_model=None):
  if not request.user.is_authenticated():
    return render_to_response('html/loggedout.html')
  
        
  if filename is not None and activity_property is not None:
    count = 0
    ms = []
    for x in Chem.SDMolSupplier(str(filename)):
        if x is not None:
          ms.append(x)
          count += 1
          if count > 10000:
            break
                    
    if count < 100:
      return HttpResponse('<h4> The training file should contain at least 100 structures<h4>')
    if count > 10000:
      return HttpResponse('<h4> Training file is too big <h4>')
    
    active,inactive,not_two = active_inactive(ms,activity_property)
    if not_two:
      return HttpResponse('<h4> Property %s should take exactly 2 values</h4>' % activity_property)
    

    ####Iwona
    #work_dir = get_dir(request)
    #(fh,model_name) = mkstemp(dir=work_dir)
    #os.close(fh);  
    #(fh,output_txt) = mkstemp(dir=work_dir)
    #os.close(fh);

    #subprocess.call(["python","/var/www/charmming/qsar/create_model.py", filename, model_name, activity_property, active, output_txt])
    #f = open(output_txt,"r")    
    #self_auc = float(f.readline())
    #recall = float(f.readline())
    #precision = float(f.readline())
    #threshold = float(f.readline())
    #auc_rand = float(f.readline())
    #cross_auc = float(f.readline())
    #f.close()
    ####end Iwona

    ###YP

    os.system("chmod g+rw %s" % (filename))
    os.system("chmod g+rw %s" % (model_name))
    os.system("chmod g+rw %s" % (output_txt))
    username=request.user.username
    u = User.objects.get(username=request.user.username)
    #job_owner_id=jobs.getNewJobOwnerIndex(u)
    #job_folder = charmming_config.user_home + '/' + username
    job_folder = charmming_config.user_qsar_jobs_home + "/qsar_job_" + common.getNewJobOwnerIndex(u)
    model_folder = charmming_config.user_qsar_models_home + "/qsar_model_" + common.getNewModelOwnerIndex(u)
    train_submitscript = open(job_folder + "/train_submitscript.inp", 'w')
    train_submitscript.write("#! /bin/bash\n")
    train_submitscript.write("cd %s\n" % (job_folder))
    #train_submitscript.write("export PYTHONPATH=$PYTHONPATH:%s/\n" % (job_folder))
    train_submitscript.write("python /var/www/charmming/qsar/create_model.py %s %s %s %s %s\n" % (filename, model_name, activity_property, active, output_txt))
    train_submitscript.write("echo 'NORMAL TERMINATION'\n")
    train_submitscript.close()
    
    scriptlist=[]
    scriptlist.append(job_folder + "/train_submitscript.inp")
    exedict={job_folder + "/train_submitscript.inp": 'sh'}
    si = schedInterface()
    newSchedulerJobID = si.submitJob(request.user.id,job_folder,scriptlist,exedict)    

    NewJob=jobs()
    NewJob.job_owner=request.user
    NewJob.job_scheduler_id=newSchedulerJobID
    NewJob.job_owner_index=NewJob.getNewJobOwnerIndex(u)
    NewJob.description = "QSAR Training Job"
    jobtype=job_types.objects.get(job_type_name='Train')
    NewJob.job_type=jobtype
    NewJob.save()

    self_auc = 0
    recall = 0
    precision = 0
    threshold = 0
    auc_rand = 0
    cross_auc = 0 
    ###end YP

    form = TestUploadForm(active=active,inactive=inactive,saved_model=model_name,threshold=threshold,activity_property=activity_property,filename=filename)
  return render_to_response('qsar/train.html', {'job_owner_index':NewJob.job_owner_index, 'self_auc': self_auc, 'auc_rand': auc_rand, 'cross_auc':cross_auc,
                                                 'threshold':threshold,'recall':recall,'precision':precision, 'form':form}, context_instance=RequestContext(request))
def predict(request):
    if not request.user.is_authenticated():
      return render_to_response('html/loggedout.html')
        
    if request.method == 'POST':
      active = request.POST['active']
      inactive = request.POST['inactive']
      saved_model = request.POST['saved_model']
      threshold = float(request.POST['threshold'])
      activity_property = request.POST['activity_property']
      filename = request.POST['filename']
      if 'back' in request.POST:
        return HttpResponseRedirect(reverse('property', args=[filename]))
      if 'predict' in request.POST:
        form = TestUploadForm(request.POST, request.FILES,active=active,inactive=inactive,saved_model=saved_model,threshold=threshold,activity_property=activity_property,filename=filename)
        if form.is_valid():
          orig_name = str(request.FILES['predictfile'].name)
          work_dir = get_dir(request)
          (fh,name) = mkstemp(dir=work_dir)
          os.close(fh);        
          handle_uploaded_file(request.FILES['predictfile'],name)
          (fh,out) = mkstemp(dir=work_dir) 
          os.close(fh);
          (fh,output_txt) = mkstemp(dir=work_dir)
          os.close(fh);
              
          subprocess.call(["python","/var/www/charmming/qsar/run_model.py",saved_model,str(name),str(threshold),active,inactive,activity_property,str(out),output_txt])
          f = open(output_txt,"r")       
          recall = float(f.readline())   
          precision = float(f.readline())
          f.close()
                                      
          output_exists = os.path.isfile(out)
          return render_to_response('qsar/predict.html', {'orig_name': orig_name, 'output_exists': output_exists, 'output': out,'recall':recall,'precision':precision,
                                                           'activity_property':activity_property,'active':active}, context_instance=RequestContext(request))
          

def download(request):
  if not request.user.is_authenticated():
    return render_to_response('html/loggedout.html')
      
  if request.method == 'POST':
    if 'download' in request.POST:
      filename = request.POST['output']
      orig_name = request.POST['orig_name']
      name,ext = fileName, fileExtension = os.path.splitext(orig_name)
      new_name = name + "-predicted"+ext
      mimetype,encoding = mimetypes.guess_type(filename)
      response = HttpResponse(mimetype=mimetype)
      response['Content-Disposition'] = 'attachment; filename=%s' % new_name
      response.write(file(filename, "rb").read())
      return response
                       
def viewJobs(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    #lesson_ok, dd_ok = checkPermissions(request)
    #if not dd_ok:
    #    return render_to_response('html/unauthorized.html')


    #try:
    dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)
    #except:
    #    pass

    jobs=[]

    cursor = dba.cursor(MySQLdb.cursors.DictCursor)
    jobssql = "SELECT js.id, dj.job_scheduler_id, jt.job_type_name, dj.job_owner_index, dj.description, js.queued, " \
              "if (state=1,timediff(now(),js.queued),timediff(started,js.queued)) 'Waited', " \
              "if (state=2,timediff(now(),js.started),timediff(ended,started)) 'Runtime', " \
              "if (state=2,timediff(now(),js.queued),timediff(ended,js.queued)) 'Total Time', " \
              "(select case state when 1 then 'Queued' when 2 then 'Running' when 3 then 'Done' when 4 then 'Failed' end) as 'Status'\
, " \
              "(select case state when 1 then 'Yellow' when 2 then 'Blue' when 3 then 'Green' when 4 then 'Red' end) as 'Status Color\
' " \
              "FROM job_scheduler js, qsar_jobs dj, qsar_job_types jt  WHERE dj.job_scheduler_id=js.id and dj.job_type_id=jt.id and " \
              " js.userid=%d order by dj.job_owner_index desc" % (request.user.id)
    cursor.execute(jobssql)
    rows = cursor.fetchall()
    cursor.close()
    dba.close()

    for row in rows:
        job=common.objQSARJob()
        job.id=row['id']
        job.scheduler_id=row['job_scheduler_id']
        job.owner_index=row['job_owner_index']
        job.description=row['description']
        job.type=row['job_type_name']
        job.status="<font color='%s'>%s</font>" % (row['Status Color'],row['Status'])
        job.queued=row['queued']
        job.waited=row['Waited']
        job.runtime=row['Runtime']
        job.totaltime=row['Total Time']
        jobs.append(job)


    lesson_ok, dd_ok = checkPermissions(request)

    return render_to_response('qsar/viewjobs.html', {'jobs': jobs,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

