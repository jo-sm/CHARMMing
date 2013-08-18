# -*- coding: utf-8 -*-
from django.db import models
from django.contrib.auth.models import User

class job_types(models.Model):

    job_type_name = models.CharField(max_length=200)
    description = models.CharField(max_length=250,null=True, blank=True)
    #application = models.ForeignKey(applications,null=True, blank=True)

class jobs(models.Model):

    job_owner = models.ForeignKey(User,related_name='job_owner')
    job_scheduler_id = models.PositiveIntegerField(default=0)
    job_owner_index = models.PositiveIntegerField(default=0)
    job_name = models.CharField(max_length=200,null=True, blank=True)
    description = models.CharField(max_length=250,null=True, blank=True)
    job_start_time = models.DateTimeField(null=True, blank=True)
    job_end_time = models.DateTimeField(null=True, blank=True)
    job_type = models.ForeignKey(job_types,null=True, blank=True)

    def getNewJobOwnerIndex(self,user):

        try:
            last_job = jobs.objects.filter(job_owner=user).order_by('-id')[0]
            next_job_owner_index=int(last_job.job_owner_index)+1

        except:
            next_job_owner_index=1

        return next_job_owner_index

class model_types(models.Model):

    model_type_name = models.CharField(max_length=200)
    description = models.CharField(max_length=250,null=True, blank=True)

class qsar_models(models.Model):

    model_owner = models.ForeignKey(User)
    model_owner_index = models.PositiveIntegerField(default=0)
    model_name = models.CharField(max_length=200,null=True, blank=True)
    model_description = models.CharField(max_length=250,null=True, blank=True)
    model_type = models.ForeignKey(model_types,null=True, blank=True)
    model_file = models.CharField(max_length=250,null=True, blank=True)
    
    def getNewModelOwnerIndex(self,user):

        try:
            last_model = qsar_models.objects.filter(model_owner=user.id).order_by('-id')[0]
            next_model_owner_index=int(last_model.model_owner_index)+1

        except:
            next_model_owner_index=1

        return next_model_owner_index

class data_types(models.Model):

    data_type = models.CharField(max_length=100)


class units(models.Model):

    unit_short_name = models.CharField(max_length=10)
    unit_name = models.CharField(max_length=100)

class attributes(models.Model):

    attribute_owner = models.ForeignKey(User,related_name='attribute_owner')
    attribute_short_name = models.CharField(max_length=50,null=True, blank=True)
    attribute_name = models.CharField(max_length=100,null=True, blank=True)
    description = models.CharField(max_length=250,null=True, blank=True)
    units = models.ForeignKey(units,null=True, blank=True)
    data_type = models.ForeignKey(data_types,null=True, blank=True)


class object_attributes(models.Model):

    object_attribute_owner = models.ForeignKey(User,related_name='object_attribute_owner')
    attribute = models.ForeignKey(attributes)
    object_table_name = models.CharField(max_length=100)
    object_id = models.PositiveIntegerField(default=0)
    value = models.CharField(max_length=250,null=True, blank=True)
