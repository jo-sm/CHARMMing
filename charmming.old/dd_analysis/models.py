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
#from scheduler.schedInterface import schedInterface
#from scheduler.statsDisplay import statsDisplay
#import charmming_config
import commands, datetime, sys, re, os, glob, smtplib
#import lesson1, normalmodes, dynamics, minimization
import string, output, charmming_config


class scoring_criteria(models.Model):

    owner = models.ForeignKey(User)
    scoring_criteria_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True)


class units(models.Model):

    unit_short_name = models.CharField(max_length=10) 
    unit_name = models.CharField(max_length=100)

class data_types(models.Model):

    data_type = models.CharField(max_length=100)


class scoring(models.Model):

    owner = models.ForeignKey(User)
    object_table_name = models.CharField(max_length=100)
    object_id = models.PositiveIntegerField(default=0)
    scoring_criteria = models.ForeignKey(scoring_criteria)
    value = models.DecimalField(max_digits=8,decimal_places=5,null=True)
    units = models.ForeignKey(units)


class attributes(models.Model):

    owner = models.ForeignKey(User)
    attribute_short_name = models.CharField(max_length=10,null=True, blank=True) 
    attribute_name = models.CharField(max_length=100,null=True, blank=True)
    description = models.CharField(max_length=250,null=True, blank=True)
    units = models.ForeignKey(units)
    data_type = models.ForeignKey(data_types)
   

class object_attributes(models.Model):

    owner = models.ForeignKey(User)
    attribute = models.ForeignKey(attributes)
    object_table_name = models.CharField(max_length=100)
    object_id = models.PositiveIntegerField(default=0)
    value = models.CharField(max_length=250,null=True, blank=True)


class interaction_types(models.Model):

    interaction_type_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250)


class object_attributes(models.Model):

    owner = models.ForeignKey(User)
    object_table_name1 = models.CharField(max_length=100)
    object_id1 = models.PositiveIntegerField(default=0)
    object_table_name2 = models.CharField(max_length=100)
    object_id2 = models.PositiveIntegerField(default=0)
    interaction_type = models.ForeignKey(interaction_types)
