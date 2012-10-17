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
#from dd_substrate import models
import charmming_config
import commands, datetime, sys, re, os, glob, smtplib
import lesson1, normalmodes, dynamics, minimization
import string, output, charmming_config
import dd_infrastructure, dd_substrate




class proteins(models.Model):

    owner = models.ForeignKey(User)
    protein_owner_index = models.PositiveIntegerField(default=0)
    pdb_code = models.CharField(max_length=5,null=True, blank=True)
    #structure = models.ForeignKey(structure)
    #pdbinfo_pdbfile = models.ForeignKey(pdbinfo_pdbfile)
    protein_name = models.CharField(max_length=200,null=True, blank=True) 
    description = models.CharField(max_length=250,null=True, blank=True)
    source = models.ForeignKey(dd_infrastructure.models.sources,null=True, blank=True)


class protein_conformations(models.Model):

    owner = models.ForeignKey(User)
    protein = models.ForeignKey(proteins)
    conformation_protein_index = models.PositiveIntegerField(default=0)
    conformation_name = models.CharField(max_length=200,null=True, blank=True)
    description = models.CharField(max_length=250,null=True, blank=True)


class projects_protein_conformations(models.Model):

    owner = models.ForeignKey(User)
    project = models.ForeignKey(dd_infrastructure.models.projects)
    protein_conformation = models.ForeignKey(protein_conformations)


class amino_acids(models.Model):
    
    amino_acid_name =  models.CharField(max_length=200)
    three_letter_code = models.CharField(max_length=3,null=True, blank=True)
    one_letter_code = models.CharField(max_length=1,null=True, blank=True)
    description = models.CharField(max_length=250,null=True, blank=True)

class protein_residues(models.Model):

    amino_acid = models.ForeignKey(amino_acids)
    protein = models.ForeignKey(proteins)
    chain = models.CharField(max_length=50,null=True, blank=True)
    residue_number = models.PositiveIntegerField(default=0)

class binding_sites(models.Model):

    protein_conformation = models.ForeignKey(protein_conformations)
    binding_site_name =  models.CharField(max_length=200)
    description = models.CharField(max_length=250,null=True, blank=True)
    
class binding_sites_residues(models.Model):

    binding_site = models.ForeignKey(binding_sites)
    residue = models.ForeignKey(protein_residues)
    
class protein_conformation_native_ligands(models.Model):

    #protein_conformation_file = models.ForeignKey(dd_infrastructure.models.files,related_name="protein_conformation_file",null=True, blank=True)
    protein_conformation = models.ForeignKey(protein_conformations)
    #ligand_file = models.ForeignKey(dd_infrastructure.models.files,related_name="ligand_file",null=True, blank=True)
    ligand = models.ForeignKey(get_model('dd_substrate', 'ligands'))
    #project = models.ForeignKey(dd_infrastructure.models.projects)
    description = models.CharField(max_length=250,null=True, blank=True)
