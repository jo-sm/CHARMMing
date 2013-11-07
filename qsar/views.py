# -*- coding: utf-8 -*-
from django.shortcuts import render_to_response
from django.template import RequestContext
from django.http import HttpResponseRedirect,HttpResponse
from django.core.urlresolvers import reverse
from qsar.forms import TrainingUploadForm,SelectProperty,TestUploadForm
from tempfile import mkstemp
import os, os.path
import subprocess
import mimetypes
import fnmatch
import urllib, urllib2
import unicodedata
import charmming_config
from django.contrib.auth.models import User
from rdkit import Chem
from qsar.train import active_inactive,get_sd_properties,check_regression_property, is_valid_sd
from qsar.models import jobs, job_types, qsar_models, model_types, jobs_models
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from account.views import checkPermissions
import MySQLdb
import MySQLdb.cursors
import settings
import common

def get_dir(request):
    username = request.user.username
    u = User.objects.get(username=username)
    location = charmming_config.user_home + '/' + request.user.username + '/qsar/'
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
    log=open("/tmp/newmodel.log", "w")
    err_message=""
    if request.method == 'POST':
        form = TrainingUploadForm(request.POST, request.FILES)
        log.write("model type:%s\n" % (request.POST['model_type']))
        if form.is_valid():
            (fh,name) = mkstemp(dir=get_dir(request))
            os.close(fh);
            handle_uploaded_file(request.FILES['trainfile'],name)
            activity_property_choice_length = len(get_sd_properties(name));
            if not activity_property_choice_length:
                err_message="Not an SD file or no SD properties found. Please upload an SD file with Activity Property"
                return render_to_response('qsar/newmodel.html', {'form': form, 'message':err_message}, context_instance=RequestContext(request)  )
            #if not is_valid_sd(name):
            #    err_message="No structures found. Please upload an SD file with at least one structure"
            #    return render_to_response('qsar/newmodel.html', {'form': form, 'message':err_message}, context_instance=RequestContext(request)  )
            #return HttpResponseRedirect(reverse('property', args=[name]))
            u = User.objects.get(username=request.user.username)
            newmodel=qsar_models()
            newmodel.model_owner=request.user
            newmodel.model_owner_index=newmodel.getNewModelOwnerIndex(u)
            newmodel.model_name=request.POST['model_name']
            newmodel.model_type=model_types.objects.get(id=request.POST['model_type'])
            newmodel.save()
            log.write("newmodel:%s" % (newmodel.id))
            return property(request,name,newmodel.id)
        else:
            return render_to_response('qsar/newmodel.html', {'form': form}, context_instance=RequestContext(request)  )
    else:
        form = TrainingUploadForm() # A empty, unbound form

    return render_to_response('qsar/newmodel.html', {'form': form}, context_instance=RequestContext(request)  )


def property(request,filename=None,qsar_model_id=None,message=""):
  if not request.user.is_authenticated():
    return render_to_response('html/loggedout.html')

  log=open("/tmp/qsar.log","w")
  log.write("model:%s" % (qsar_model_id))
  if request.method == 'POST':
    if 'filename' in request.POST:
      filename = request.POST['filename']
      #model_type=model_types.objects.get(id=qsar_model.model_type_id)
      qsar_model_id=request.POST['qsar_model_id']
      log.write("model set\n")
    form = SelectProperty(request.POST,filename=filename,qsar_model_id=qsar_model_id)
    if 'back' in request.POST:
      return HttpResponseRedirect(reverse('newModel', args=()))
      log.write("back\n")
    elif 'next' in request.POST:
      if form.is_valid():
        log.write("next\n")
        activity_property = request.POST["activity_property"]
        #qsar_model_id=request.POST['qsar_model_id']
        #common.AssignObjectAttribute(request.user.id,qsar_model_id,"qsar_qsar_models","Activity Property",activity_property)
        #model_type=model_types.objects.get(id=qsar_model.model_type_id)
        return train(request,filename,activity_property,qsar_model_id)
  else:
    model_type=model_types.objects.get(id=qsar_model.model_type_id)
    form = SelectProperty(filename=filename,qsar_model_id=qsar_model_id)
  activity_property_choice_length = 0;
  if filename is not None:
    qsar_model=qsar_models.objects.get(id=qsar_model_id)
    model_type=model_types.objects.get(id=qsar_model.model_type_id)
    activity_property_choice_length = len(get_sd_properties(filename));
    #if notactivity_property_choice_length
  return render_to_response('qsar/property.html', {'form': form,'filename':filename,'activity_property_choice_length':activity_property_choice_length, 'qsar_model':qsar_model, 'model_type':model_type}, context_instance=RequestContext(request))

def train(request,filename=None,activity_property=None,qsar_model_id=None):
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
    
    qsar_model=qsar_models.objects.get(id=qsar_model_id)
    err_message=""

    if qsar_model.model_type.model_type_name=="Random Forest (SAR) Categorization":
        active,inactive,not_two = active_inactive(ms,activity_property)
        if not_two:
            err_message="Property %s should take exactly 2 values" % activity_property

        
    if qsar_model.model_type.model_type_name=="Random Forest (QSAR) Regression":
        if not check_regression_property(ms,activity_property):
            err_message="Property %s is not numerical or takes fewer than 10 values" % activity_property


    if err_message !="":
        model_type=model_types.objects.get(id=qsar_model.model_type_id)
        form = SelectProperty(filename=filename,qsar_model_id=qsar_model.id)
        activity_property_choice_length = 0;
        activity_property_choice_length = len(get_sd_properties(filename))
        return render_to_response('qsar/property.html', {'form': form,'filename':filename,'activity_property_choice_length':activity_property_choice_length, 'qsar_model':qsar_model, 'model_type':model_type, 'message':err_message}, context_instance=RequestContext(request))


    common.AssignObjectAttribute(request.user.id,qsar_model_id,"qsar_qsar_models","Activity Property",activity_property)
    

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
    log=open("/tmp/qjob.log","w")
    
    
    username=request.user.username
    u = User.objects.get(username=request.user.username)
    #job_owner_id=jobs.getNewJobOwnerIndex(u)
    #job_folder = charmming_config.user_home + '/' + username
    NewJob=jobs()
    qsar_folder = charmming_config.user_home + "/" + username + "/qsar"
    qsar_jobs_folder = qsar_folder + "/jobs"
    job_folder = qsar_folder + "/jobs/qsar_job_" + str(NewJob.getNewJobOwnerIndex(u))
    qsar_models_folder = qsar_folder + "/models"
    model_folder = qsar_folder + "/models/qsar_model_" + str(qsar_model.model_owner_index) #common.getNewModelOwnerIndex(u)
    os.system("mkdir %s" % (qsar_folder))
    log.write("mkdir %s\n" % (qsar_folder))
    os.system("mkdir %s" % (qsar_jobs_folder))
    os.system("mkdir %s" % (job_folder))
    os.system("chmod -R g+w %s" % (job_folder))
    log.write("mkdir %s\n" % (job_folder))
    os.system("mkdir %s" % (qsar_models_folder))
    os.system("mkdir %s" % (model_folder))
    os.system("chmod -R g+w %s" % (model_folder))
    log.write("mkdir %s\n" % (model_folder))
    model_file=model_folder + "/model"
    training_output=model_folder + "/training_output"
    os.system("chmod g+rw %s" % (filename))
    os.system("chmod g+rw %s" % (model_file))
    os.system("chmod g+rw %s" % (training_output))
    train_submitscript = open(job_folder + "/train_submitscript.inp", 'w')
    train_submitscript.write("#! /bin/bash\n")
    train_submitscript.write("cd %s\n" % (job_folder))
    train_submitscript.write("export PYTHONPATH=$PYTHONPATH:%s/\n" % ("/var/www/charmming")) #job_folder))
    train_submitscript.write("export DJANGO_SETTINGS_MODULE=settings\n")
    log.write("model type: %s\n" % (qsar_model.model_type.model_type_name))
    if qsar_model.model_type.model_type_name=="Random Forest (SAR) Categorization":
        train_submitscript.write("python /var/www/charmming/qsar/create_model.py %s %s %s %s %s %s %s\n" % (filename, model_file, activity_property, active, training_output, str(qsar_model.id), str(u.id)))
    elif qsar_model.model_type.model_type_name=="Random Forest (QSAR) Regression":
        train_submitscript.write("python /var/www/charmming/qsar/create_model_regression.py %s %s %s %s %s %s\n" % (filename, model_file, activity_property, training_output, str(qsar_model.id), str(u.id)))
    
    train_submitscript.write("echo 'NORMAL TERMINATION'\n")
    train_submitscript.close()
    
    scriptlist=[]
    scriptlist.append(job_folder + "/train_submitscript.inp")
    exedict={job_folder + "/train_submitscript.inp": 'sh'}
    si = schedInterface()
    newSchedulerJobID = si.submitJob(request.user.id,job_folder,scriptlist,exedict)    


    NewJob.job_owner=request.user
    NewJob.job_scheduler_id=newSchedulerJobID
    NewJob.job_owner_index=NewJob.getNewJobOwnerIndex(u)
    NewJob.description = "Building a model " + qsar_model.model_name #"QSAR Training Job"
    try:
        jobtype=job_types.objects.get(job_type_name='Train')
    except:
        jobtype=job_types()
        jobtype.job_type_name="Train"
        jobtype.save()
  
    NewJob.job_type=jobtype
    NewJob.save()

    
    qsar_model.model_file=model_file
    if qsar_model.model_type.model_type_name=="Random Forest (SAR) Categorization":
        common.AssignObjectAttribute(request.user.id,qsar_model.id,"qsar_qsar_models","Active",active)
        common.AssignObjectAttribute(request.user.id,qsar_model.id,"qsar_qsar_models","Inactive",inactive)
    
    qsar_model.save()
    
    NewJobModel=jobs_models()
    NewJobModel.owner = request.user
    NewJobModel.job_id = NewJob.id
    NewJobModel.qsar_model_id=qsar_model.id
    NewJobModel.save()

    self_auc = 0
    recall = 0
    precision = 0
    threshold = 0
    auc_rand = 0
    cross_auc = 0 
    ###end YP

    #form = TestUploadForm(active=active,inactive=inactive,saved_model=model_file,threshold=threshold,activity_property=activity_property,filename=filename)
  #return render_to_response('qsar/train.html', {'job_owner_index':NewJob.job_owner_index, 'self_auc': self_auc, 'auc_rand': auc_rand, 'cross_auc':cross_auc,
  #                                               'threshold':threshold,'recall':recall,'precision':precision, 'form':form}, context_instance=RequestContext(request))
  return viewJobs(request)

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


def newPredict(request,model_id):
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

    return render_to_response('qsar/viewjobs.html')
    
def viewJobsDiv(request):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    #lesson_ok, dd_ok = checkPermissions(request)
    #if not dd_ok:
    #    return render_to_response('html/unauthorized.html')


    #try:
    dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db="charmming2", compress=1)
    #except:
    #    pass

    jobs=[]
    log=open("/tmp/viewjobs.log","w")
    cursor = dba.cursor(MySQLdb.cursors.DictCursor)
    jobssql = "SELECT dj.id, dj.job_scheduler_id, jt.job_type_name, dj.job_owner_index, dj.description, js.queued, " \
              "if (state=1,timediff(now(),js.queued),timediff(started,js.queued)) 'Waited', " \
              "if (state=2,timediff(now(),js.started),timediff(ended,started)) 'Runtime', " \
              "if (state=2,timediff(now(),js.queued),timediff(ended,js.queued)) 'Total Time', " \
              "(select case state when 1 then 'Queued' when 2 then 'Running' when 3 then 'Done' when 4 then 'Failed' end) as 'Status', " \
              "(select case state when 1 then 'Yellow' when 2 then 'Blue' when 3 then 'Green' when 4 then 'Red' end) as 'Status Color' " \
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
        #try:
        output_data=""
        log.write("looking for output data\n")
        try:
        #if 1==1:
            if job.type=="Predict":
                log.write("predict job found\n")
                qsar_folder = charmming_config.user_home + "/" + request.user.username + "/qsar"
                qsar_jobs_folder = qsar_folder + "/jobs"
                job_folder = qsar_folder + "/jobs/qsar_job_" + str(job.owner_index)
                log.write("opening predict output: %s\n" % (job_folder + "/predict_output"))
                #log.write("does file exist:" + str(os.path.isfile(job_folder + "/predict_output"))+"\n")
                #if os.path.isfile(job_folder + "/predict_output"):
                predict_output=open(job_folder + "/predict_output","r")
                log.write("predict output opened\n")
                for line in predict_output:
                    log.write("reading predit output value: %s\n" % (line))
                    log.write("converted value: %s\n" % str(round(float(line),3)))
                    output_data=output_data+str(round(float(line.strip("\n")),3))+"/"
                    
                output_data=output_data[:-1]
                predict_output.close()

        except:
            output_data=""
        
        job.output_data=output_data
        log.write("job output data: %s\n" % (job.output_data))
        try:
            log.write("owner_id: %s, job_id: %s\n" % (request.user.id, job.id))
            qsar_job_model=jobs_models.objects.get(owner_id=request.user.id, job_id = job.id)
            job.model_id=qsar_job_model.qsar_model_id
        except:
            job.model_id=0
        
        jobs.append(job)


    lesson_ok, dd_ok = checkPermissions(request)

    return render_to_response('qsar/viewjobsdiv.html', {'jobs': jobs,'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

def viewModels(request,message=""):
    if not request.user.is_authenticated():
        return render_to_response('html/loggedout.html')

    dba = MySQLdb.connect("localhost", user=settings.DATABASE_USER, passwd=settings.DATABASE_PASSWORD, \
                             db=settings.DATABASE_NAME, compress=1)

    models=[]

    #cursor = dba.cursor(MySQLdb.cursors.DictCursor)
    
    #modelssql "select qm.*, qmt.model_type_name from qsar_qsar_models qm left join qsar_model_types qmt on " \
    #          " qm.model_type_id=qmt.id where model_owner_id=%s" % (request.user.id)

    #cursor.execute(sql)
    #rows = cursor.fetchall()
    #cursor.close()
    #dba.close()
    orig_name=""
    username=request.user.username
    log=open("/tmp/predictjob.log","w")
    if request.method == 'POST' and message=="":
        models_list = qsar_models.objects.select_related().filter(model_owner_id=request.user.id)
        for qsar_model in models_list:
            try:  
                if request.FILES['predict_file_' + str(qsar_model.id)]:
                    log.write("file: predict_file_" + str(qsar_model.id))
                    orig_name=str(request.FILES['predict_file_' + str(qsar_model.id)].name)
            except:
                continue      
            
            name,ext=fileName, fileExtension = os.path.splitext(orig_name)
            new_name=name + "-predicted"+ext
            u = User.objects.get(username=request.user.username)
            #job_owner_id=jobs.getNewJobOwnerIndex(u) 
            #job_folder = charmming_config.user_home + '/' + username
            NewJob=jobs()
            qsar_folder = charmming_config.user_home + "/" + username + "/qsar"
            qsar_jobs_folder = qsar_folder + "/jobs"
            job_folder = qsar_folder + "/jobs/qsar_job_" + str(NewJob.getNewJobOwnerIndex(u))
            os.system("mkdir %s" % (qsar_folder))
            log.write("mkdir %s\n" % (qsar_folder))
            os.system("mkdir %s" % (qsar_jobs_folder))
            os.system("mkdir %s" % (job_folder))
            os.system("chmod -R g+w %s" % (job_folder))
            log.write("mkdir %s\n" % (job_folder))
            predict_output_file=job_folder + "/predict_output"
            predict_results_file=job_folder +"/" +  new_name #"/predict_results.sdf"
            #predict_input_file=request.POST['predict_file_' + str(qsar_model.id)]
            predict_input_file = open("/%s/predict_input.sdf" % (job_folder), 'w')
            for fchunk in request.FILES['predict_file_' + str(qsar_model.id)].chunks():
                predict_input_file.write(fchunk)
            predict_input_file.close()
            predict_input_file="/%s/predict_input.sdf" % (job_folder)
            if not is_valid_sd(predict_input_file):
                err_message="No structures found in prediction file. Please upload an SD file with at least one structure"
                return viewModels(request,err_message)
            
            log.write("wrote predict input\n")
            #os.system("chmod g+rw %s" % (predict_file))
            os.system("chmod g+rw %s" % (predict_results_file))
            os.system("chmod g+rw %s" % (predict_output_file))
            predict_submitscript = open(job_folder + "/predict_submitscript.inp", 'w')
            predict_submitscript.write("#! /bin/bash\n")
            predict_submitscript.write("cd %s\n" % (job_folder))
            predict_submitscript.write("export PYTHONPATH=$PYTHONPATH:%s/\n" % ("/var/www/charmming")) #job_folder))
            predict_submitscript.write("export DJANGO_SETTINGS_MODULE=settings\n")
            #threshold=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Recommended threshold")
            #activity_property=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Activity Property")
            #active=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Active")
            #inactive=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Inactive")
            #subprocess.call(["python","/var/www/charmming/qsar/run_model.py",saved_model,str(name),str(threshold),active,inactive,activity_property,str(out),output_txt])
            log.write("model type: %s" % (qsar_model.model_type.model_type_name))
            if qsar_model.model_type.model_type_name=="Random Forest (SAR) Categorization":
                threshold=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Recommended threshold")
                activity_property=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Activity Property")
                active=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Active")
                inactive=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Inactive")

                predict_submitscript.write("python /var/www/charmming/qsar/run_model.py %s %s %s %s %s %s %s %s\n" %  
                                      (qsar_model.model_file, predict_input_file, threshold, active, inactive, activity_property, predict_results_file, predict_output_file))
            elif qsar_model.model_type.model_type_name=="Random Forest (QSAR) Regression":
                log.write("regression model executing\n")
                activity_property=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Activity Property")
                #common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "Self R2", str(self_r2))
                #common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "Y-randomization", str(r2_rand))
                #common.AssignObjectAttribute(user_id, model_id, "qsar_qsar_models", "5-fold cross-validation", str(cross_r2))

                predict_submitscript.write("python /var/www/charmming/qsar/run_model_regression.py %s %s %s %s %s\n" %
                                      (qsar_model.model_file, predict_input_file, activity_property, predict_results_file, predict_output_file))

            predict_submitscript.write("echo 'NORMAL TERMINATION'\n")
            predict_submitscript.close()

            scriptlist=[]
            log.write("submitting predict script\n")
            scriptlist.append(job_folder + "/predict_submitscript.inp")
            exedict={job_folder + "/predict_submitscript.inp": 'sh'}
            si = schedInterface()
            newSchedulerJobID = si.submitJob(request.user.id,job_folder,scriptlist,exedict)

            NewJob.job_owner=request.user
            NewJob.job_scheduler_id=newSchedulerJobID
            NewJob.job_owner_index=NewJob.getNewJobOwnerIndex(u)
            NewJob.description = "Based on " + qsar_model.model_name#"QSAR Prediction Job"
            log.write("creating a new job object\n")
            try:
                jobtype=job_types.objects.get(job_type_name='Predict')
            except:
                jobtype=job_types()
                jobtype.job_type_name="Predict"
                jobtype.save()

            NewJob.job_type=jobtype
            NewJob.save()
            log.write("redirecting to submit_predict\n")
            return viewJobs(request)
            #return render_to_response ('qsar/submit_predict.html', {'job_owner_index' : NewJob.job_owner_index})
            #return HttpResponse("Predict Job %s was successfully submitted." % (NewJob.job_owner_index))
          #except:
          #  continue
        if orig_name=="":
                err_message="No prediction file specified. Please upload an SD file with at least one structure"
                return viewModels(request,err_message)
    else:
        models_list = qsar_models.objects.select_related().filter(model_owner_id=request.user.id).order_by("-model_owner_index")
        models_with_attributes_list=[]
        for qsar_model in models_list:
        
            model_with_attributes=common.ObjectWithAttributes()
            model_with_attributes.object_id=qsar_model.id
            model_with_attributes.object_owner_id=qsar_model.model_owner_id
            model_with_attributes.object_owner_index=qsar_model.model_owner_index
            model_with_attributes.object_name=qsar_model.model_name
            model_with_attributes.object_type=qsar_model.model_type.model_type_name
            model_with_attributes.FillObjectAttributeList("qsar_qsar_models")
        
            qsar_folder = charmming_config.user_home + "/" + username + "/qsar"
            qsar_models_folder = qsar_folder + "/models"
            model_folder = qsar_folder + "/models/qsar_model_" + str(qsar_model.model_owner_index)

            training_output = model_folder + "/training_output"
            model_with_attributes.object_status="<font color='Red'>Incomplete</font>"
            if os.path.exists(training_output):
                if os.path.getsize(training_output)>0:
                    model_with_attributes.object_status="<font color='Green'>Ready</font>"

            models_with_attributes_list.append(model_with_attributes)

        lesson_ok, dd_ok = checkPermissions(request)

        return render_to_response('qsar/viewmodels.html', {'models': models_with_attributes_list, 'message':message, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})

def viewModelDetails(request,qsar_model_id,message=""):
    if not request.user.is_authenticated():
        return render_to_reponse('html/loggedout.html')

    username = request.user.username

    u = User.objects.get(username=username)

    qsar_model = qsar_models.objects.get(id=qsar_model_id,model_owner_id=request.user.id)
    log=open("/tmp/modeldetails.log","w")
    if request.method == 'POST' and message=="":
        
        log.write("file: predict_file_" + str(qsar_model_id))

        u = User.objects.get(username=request.user.username)
        #job_owner_id=jobs.getNewJobOwnerIndex(u)                                                                                                                                                    
        #job_folder = charmming_config.user_home + '/' + username                                                                                                                                    
        NewJob=jobs()
        qsar_folder = charmming_config.user_home + "/" + username + "/qsar"
        qsar_jobs_folder = qsar_folder + "/jobs"
        job_folder = qsar_folder + "/jobs/qsar_job_" + str(NewJob.getNewJobOwnerIndex(u))
        os.system("mkdir %s" % (qsar_folder))
        log.write("mkdir %s\n" % (qsar_folder))
        os.system("mkdir %s" % (qsar_jobs_folder))
        os.system("mkdir %s" % (job_folder))
        os.system("chmod -R g+w %s" % (job_folder))
        log.write("mkdir %s\n" % (job_folder))
        predict_output_file=job_folder + "/predict_output"
        try:
            orig_name=str(request.FILES['predict_file'].name)
        except:
            err_message="No prediction file specified. Please upload an SD file with at least one structure"
            return viewModelDetails(request,qsar_model_id,err_message)
        name,ext=fileName, fileExtension = os.path.splitext(orig_name)
        new_name=name + "-predicted"+ext
        predict_results_file=job_folder + "/" + new_name #"/predict_results.sdf"
        #predict_input_file=request.POST['predict_file_' + str(qsar_model.id)]                                                                                                                       
        predict_input_file = open("/%s/predict_input.sdf" % (job_folder), 'w')
        #for fchunk in request.FILES['predict_file_' + str(qsar_model.id)].chunks():
        for fchunk in request.FILES['predict_file'].chunks():
            predict_input_file.write(fchunk)
        predict_input_file.close()
        predict_input_file="/%s/predict_input.sdf" % (job_folder)
        log.write("wrote predict input\n")
        #os.system("chmod g+rw %s" % (predict_file))
        if not is_valid_sd(predict_input_file):
                err_message="No structures found in prediction file. Please upload an SD file with at least one structure"
                return viewModelDetails(request,qsar_model_id,err_message)
                                      
        os.system("chmod g+rw %s" % (predict_results_file))
        os.system("chmod g+rw %s" % (predict_output_file))
        predict_submitscript = open(job_folder + "/predict_submitscript.inp", 'w')
        predict_submitscript.write("#! /bin/bash\n")
        predict_submitscript.write("cd %s\n" % (job_folder))
        predict_submitscript.write("export PYTHONPATH=$PYTHONPATH:%s/\n" % ("/var/www/charmming")) #job_folder))                                                                                     
        predict_submitscript.write("export DJANGO_SETTINGS_MODULE=settings\n")
        #threshold=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Recommended threshold")
        #activity_property=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Activity Property")
        #active=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Active")
        #inactive=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Inactive")
        #subprocess.call(["python","/var/www/charmming/qsar/run_model.py",saved_model,str(name),str(threshold),active,inactive,activity_property,str(out),output_txt])                               
        #predict_submitscript.write("python /var/www/charmming/qsar/run_model.py %s %s %s %s %s %s %s %s\n" %
        #                          (qsar_model.model_file, predict_input_file, threshold, active, inactive, activity_property, predict_results_file, predict_output_file))
        if qsar_model.model_type.model_type_name=="Random Forest (SAR) Categorization":
            
            threshold=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Recommended threshold")
            activity_property=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Activity Property")
            active=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Active")
            inactive=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Inactive")
            predict_submitscript.write("python /var/www/charmming/qsar/run_model.py %s %s %s %s %s %s %s %s\n" %
                                      (qsar_model.model_file, predict_input_file, threshold, active, inactive, activity_property, predict_results_file, predict_output_file))
        elif qsar_model.model_type.model_type_name=="Random Forest (QSAR) Regression":
            activity_property=common.GetObjectAttributeValue(request.user.id, qsar_model.id, "qsar_qsar_models", "Activity Property")
            predict_submitscript.write("python /var/www/charmming/qsar/run_model_regression.py %s %s %s %s %s\n" %
                                  (qsar_model.model_file, predict_input_file, activity_property, predict_results_file, predict_output_file))

        predict_submitscript.write("echo 'NORMAL TERMINATION'\n")
        predict_submitscript.close()

        scriptlist=[]
        log.write("submitting predict script\n")
        scriptlist.append(job_folder + "/predict_submitscript.inp")
        exedict={job_folder + "/predict_submitscript.inp": 'sh'}
        si = schedInterface()
        newSchedulerJobID = si.submitJob(request.user.id,job_folder,scriptlist,exedict)

        NewJob.job_owner=request.user
        NewJob.job_scheduler_id=newSchedulerJobID
        NewJob.job_owner_index=NewJob.getNewJobOwnerIndex(u)
        NewJob.description = "Based on " + qsar_model.model_name#"QSAR Prediction Job"
        log.write("creating a new job object\n")
        try:
            jobtype=job_types.objects.get(job_type_name='Predict')
        except:
            jobtype=job_types()
            jobtype.job_type_name="Predict"
            jobtype.save()

        NewJob.job_type=jobtype
        NewJob.save()
        log.write("redirecting to submit_predict\n")
        return viewJobs(request)
        #return 
        #return render_to_response ('qsar/submit_predict.html', {'job_owner_index' : NewJob.job_owner_index})

    else:
    
        model_with_attributes=common.ObjectWithAttributes()
        model_with_attributes.object_id=qsar_model.id
        model_with_attributes.object_owner_id=qsar_model.model_owner_id
        model_with_attributes.object_owner_index=qsar_model.model_owner_index
        model_with_attributes.object_name=qsar_model.model_name
        model_with_attributes.object_type=qsar_model.model_type.model_type_name
        model_with_attributes.FillObjectAttributeList("qsar_qsar_models")

        qsar_folder = charmming_config.user_home + "/" + username + "/qsar"
        qsar_models_folder = qsar_folder + "/models"
        model_folder = qsar_folder + "/models/qsar_model_" + str(qsar_model.model_owner_index)

        training_output = model_folder + "/training_output"
        model_with_attributes.object_status="<font color='Red'>Incomplete</font>"
        if os.path.exists(training_output):
            if os.path.getsize(training_output)>0:
                model_with_attributes.object_status="<font color='Green'>Ready</font>"

        

        lesson_ok, dd_ok = checkPermissions(request)

        return render_to_response('qsar/viewmodeldetails.html', {'model': model_with_attributes, 'message':message, 'lesson_ok': lesson_ok, 'dd_ok': dd_ok})


def downloadQSARResults(request,qsar_job_id,mimetype=None):
    if not request.user.is_authenticated():
        return render_to_reponse('html/loggedout.html')

    username = request.user.username

    u = User.objects.get(username=username)
    log=open("/tmp/qsardownload.log","w")

    qsar_folder = charmming_config.user_home + "/" + username + "/qsar"
    #qsar_folder ="/schedd/" + username + "/qsar"
    qsar_jobs_folder = qsar_folder + "/jobs"
    job_folder = qsar_jobs_folder + "/qsar_job_" + qsar_job_id
    
   
    log.write("calling: " + "ls -1 " + job_folder+"/*-predicted.*\n")
    p = subprocess.Popen("ls -1 " + job_folder+"/*-predicted.*", shell=True, stdout=subprocess.PIPE)
    pout, perr = p.communicate()
    log.write ("pout is: %s\n" % (perr))
    filename=pout.rstrip()
    
    log.write("filename is: %s\n" % (filename))

    #if mimetype is None:
    mimetype, encoding = mimetypes.guess_type("%s" % (filename))
    try:
        os.stat("%s" % (filename))
    except:
        return HttpResponse('Oops ... that file no longer exists.')
    
    
    response = HttpResponse(mimetype=mimetype)
    response['Content-Disposition'] = 'attachment; filename=%s' % (os.path.basename(filename))
    response.write(file(filename,"rb").read())
    return response


