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
from django.db.models import get_model
from django.template import *
from django.template.loader import get_template
from django.db import models
from django import forms
from django.contrib.auth.models import User
from django.core.mail import mail_admins
from scheduler.schedInterface import schedInterface
from scheduler.statsDisplay import statsDisplay
#from dd_target import models
#from dd_infrastructure import models
import charmming_config
import commands, datetime, sys, re, os, glob, smtplib
import lesson1, normalmodes, dynamics, minimization
import string, output, charmming_config
import dd_infrastructure, dd_target


class fragments(models.Model):

    owner = models.ForeignKey(User)
    fragment_owner_index = models.PositiveIntegerField(default=0)
    fragment_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True) 
    source = models.ForeignKey(dd_infrastructure.models.sources)

class fragment_sets(models.Model):

    owner = models.ForeignKey(User)
    fragment_set_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True) 


class fragment_sets_fragments(models.Model):

    owner = models.ForeignKey(User)
    fragment_set = models.ForeignKey(fragment_sets)
    fragment = models.ForeignKey(fragments)


class ligands(models.Model):

    owner = models.ForeignKey(User)
    ligand_owner_index = models.PositiveIntegerField(default=0)
    ligand_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True) 
    source = models.ForeignKey(dd_infrastructure.models.sources)


class ligand_sets(models.Model):

    owner = models.ForeignKey(User)
    ligand_set_name = models.CharField(max_length=200) 
    description = models.CharField(max_length=250,null=True, blank=True)
    public = models.CharField(max_length=1)
     


class ligand_sets_ligands(models.Model):

    owner = models.ForeignKey(User)
    ligands_set = models.ForeignKey(ligand_sets)
    ligands = models.ForeignKey(ligands) 


class poses(models.Model):
    
    pose_object_table_name = models.CharField(max_length=100)
    pose_object_id = models.PositiveIntegerField(default=0)
    pose_name = models.CharField(max_length=200)
    description = models.CharField(max_length=250,null=True, blank=True)
    source = models.ForeignKey(dd_infrastructure.models.sources)
   
class poses_binding_sites(models.Model):
    
    owner = models.ForeignKey(User)
    pose = models.ForeignKey(poses)
    binding_site = models.ForeignKey(get_model('dd_target', 'binding_sites'))

