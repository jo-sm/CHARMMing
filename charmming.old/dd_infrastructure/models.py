#
#                            PUBLIC DOMAIN NOTICE
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the authors' official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software is freely available
#  to the public for use.  There is no restriction on its use or
#  reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, NIH and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. NIH, NHLBI, and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any
#  particular purpose.
from django.template import *
from django.template.loader import get_template
from django.db import models
from django import forms
from django.contrib.auth.models import User
from django.core.mail import mail_admins
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
from structure.models import Task
import charmming_config
import commands, datetime, sys, re, os, glob, smtplib
import lesson1, normalmodes, dynamics, minimization
import string, output, charmming_config


class projects(models.Model):

    owner = models.ForeignKey(User)
    project_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True)
    date_created = models.DateTimeField(null=True, blank=True)
    selected = models.CharField(max_length=1)


class applications(models.Model):

    application_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True)


class job_types(models.Model):

    job_type_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True)
    application = models.ForeignKey(applications,null=True, blank=True)


class jobs(models.Model):

    owner = models.ForeignKey(User)
    job_scheduler_id = models.PositiveIntegerField(default=0)
    job_owner_index = models.PositiveIntegerField(default=0) 
    job_name = models.CharField(max_length=200,null=True, blank=True) 
    description = models.CharField(max_length=250,null=True, blank=True)
    job_start_time = models.DateTimeField(null=True, blank=True)
    job_end_time = models.DateTimeField(null=True, blank=True)
    job_type = models.ForeignKey(job_types,null=True, blank=True)


class files(models.Model):

    owner = models.ForeignKey(User)
    file_name = models.CharField(max_length=100)
    file_location = models.CharField(max_length=300,null=True, blank=True)
    description = models.CharField(max_length=250,null=True, blank=True)


class files_objects(models.Model):

    owner = models.ForeignKey(User)
    file = models.ForeignKey(files)
    object_table_name = models.CharField(max_length=100)
    object_id = models.PositiveIntegerField(default=0)
 

class sources(models.Model):

    source_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True)
    source_object_table_name = models.CharField(max_length=100,null=True, blank=True)
    source_object_id = models.PositiveIntegerField(default=0,null=True, blank=True)

class file_types(models.Model):

    file_type_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True)


class ddtask(Task):
    
    def start(self,job_folder,scriptlist,exedict):
        st = self.workstruct.structure

        #logfp = open('/tmp/startjob.txt', 'w')
        #logfp.write('In job start routine.\n')

        si = schedInterface()
        #if kwargs.has_key('altexe'):
        #    logfp.write('Got alt exe.\n')
        #    exedict = {}
        #    for inpscript in self.scriptList:
        #        exedict[inpscript] = kwargs['altexe']
        #    logfp.write('exedict = %s\n' % exedict)
        self.jobID = si.submitJob(st.owner.id,job_folder,scriptlist,exedict)
        #else:
        #    logfp.write('No alt exe.\n')
        #    self.jobID = si.submitJob(st.owner.id,st.location,self.scriptList)
        #logfp.close()
        if self.jobID > 0:
            self.save()
            self.query()
        else:
            raise AssertionError('Job submission fails')    

class dockingtask(ddtask):
    pass

class dsfdockingtask(dockingtask):
    pass




"""
class projects(models.Model):

    author = models.CharField(max_length=250) 
    journal = models.CharField(max_length=250)
    pub_date = models.DateTimeField(default=datetime.datetime.now)
    rtf = models.CharField(max_length=250)  
    rtf_append_replace = models.CharField(max_length=250,null=True)  
    prm = models.CharField(max_length=250)  
    prm_append_replace = models.CharField(max_length=250,null=True)  
    selected = models.CharField(max_length=1)  
    good_het = models.CharField(max_length=250)  
    nongood_het = models.CharField(max_length=250)  
    segids = models.CharField(max_length=250)
    pid = models.CharField(max_length=6)
    append_status = models.CharField(max_length=250) 
    solvation_structure = models.CharField(max_length=50) 
    crystal_x = models.DecimalField(max_digits=8,decimal_places=5,null=True)
    crystal_y = models.DecimalField(max_digits=8,decimal_places=5,null=True)
    crystal_z = models.DecimalField(max_digits=8,decimal_places=5,null=True)
    domains = models.CharField(max_length=250,default='')
    fes4list = models.CharField(max_length=250,default='')
    doblncharge = models.CharField(max_length=1,default='f')
    
    minimization_jobID = models.PositiveIntegerField(default=0)
    solvation_jobID = models.PositiveIntegerField(default=0)
    nma_jobID = models.PositiveIntegerField(default=0)
    ld_jobID = models.PositiveIntegerField(default=0)
    md_jobID = models.PositiveIntegerField(default=0)
    sgld_jobID = models.PositiveIntegerField(default=0)
    redox_jobID = models.PositiveIntegerField(default=0)
"""
    
