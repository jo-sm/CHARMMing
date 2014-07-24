from django.views import generic
from django.http import HttpResponse, HttpResponseRedirect
from charmming.models import Structure
from views.mixins import PermissionsMixin
from views.helpers import id_generator
import datetime
import logging
import json
import os
from pychm.io.pdb import PDBFile

class StructureIndexView(PermissionsMixin, generic.TemplateView):
  template_name = 'structure/index.html'

  def get(self, request, *args, **kwargs):
    # view all structures
    context = self.get_context_data()
    context['structures'] = Structure.objects.filter(user=context['current_user'])
    return self.render_to_response(context)

class StructureShowView(PermissionsMixin, generic.TemplateView):
  template_name = 'structure/show.html'

  def get(self, request, *args, **kwargs):
    return self.render_to_response([])

  def delete(self, request, *args, **kwargs):
    id = kwargs.get('id')
    Structure.objects.filter(user=self.get_current_user(request), id=id).delete()
    return HttpResponse(json.dumps({'success': True}), content_type="application/json", status=200)

class StructureShowSegmentsView(PermissionsMixin, generic.TemplateView):
  # Currently evaluates the requested structure and, on the first request, writes the segments to 
  # the filesystem for easy retrieval later. This is a workaround to using a background processor to do this 
  # when the structure is first uploaded, and will be implemented in a background processing
  # queue e.g. Celery.
  def get(self, request, *args, **kwargs):
    id = kwargs.get('id')
    current_user = self.get_current_user(request)
    structure = Structure.objects.filter(user=current_user, id=id)
    if structure:
      structure = structure[0]
      
      segments = {}

      filename = '/home/schedd/' + current_user.username + '/' + structure.directory + '/' + structure.file_name
      pdb = PDBFile(filename)
      for model in pdb.iter_models():
        dir = '/home/schedd/' + current_user.username + '/' + structure.directory + '/segments/' + model.name
        if not os.path.isdir(dir):
          os.mkdir(dir)
        segments[model.name] = {}
    #    model.parse() # Save time
        for segment in model.iter_seg():
          if not os.path.exists(dir + '/' + segment.segid + '.pdb'):
            segment.write(dir + '/' + segment.segid + '.pdb', outformat='charmm')
          segments[model.name][segment.segid] = {}
          segments[model.name][segment.segid]['segmentType'] = segment.segType

          default_first_patches = ['NTER','GLYP','PROP','ACE','ACP']
          first_patches = []
          last_patches = []
          first_residue = segment.iter_res().next()
          # Determine the list of patches
          if segment.segType == 'pro':
            # Sets don't preserve order, but we want to in this scenario
            # We mimic the functionality of a set for the protein residue patch list
            last_patches = ['CTER','CT1','CT2','CT3','ACP']
            if first_residue == 'gly':
              first_patches.append('GLYP')
            elif first_residue == 'pro':
              first_patches.append('PROP')
            else:
              first_patches.append('NTER')
            default_first_patches.remove(first_patches[0])
            first_patches.extend(default_first_patches)
          elif segment.segType == 'dna' or segment.segType == 'rna':
            first.patches.update(('5TER','5MET','5PHO','5POM','CY35'))
            last_patches = ['3TER','3PHO','3POM','3CO3','CY35']

          if len(first_patches):
            first_patches.append('NONE')
          if len(last_patches):
            last_patches.append('NONE')
          segments[model.name][segment.segid]['firstPatches'] = first_patches
          segments[model.name][segment.segid]['lastPatches'] = last_patches
          if segment.segType == 'bad':
            segments[model.name][segment.segid]['residue'] = segment.iter_res().next().resName
            segments[model.name][segment.segid]['additionalToppar'] = {}
            if current_user.is_admin:
              segments[model.name][segment.segid]['additionalToppar']['dogmans'] = 'Force Dogmans'
              segments[model.name][segment.segid]['additionalToppar']['match'] = 'Force MATCH'
              segments[model.name][segment.segid]['additionalToppar']['ante'] = 'Force Antechamber'
              segments[model.name][segment.segid]['additionalToppar']['genrtf'] = 'Force GENRTF'
      return HttpResponse(json.dumps(segments), content_type='application/json', status=200)
    else:
      return HttpResponse(status=404)

class StructureNewView(PermissionsMixin, generic.TemplateView):
  template_name = 'structure/new.html'

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    return self.render_to_response(context)

  def post(self, request, *args, **kwargs):
    # User submits a new structure
    inputs = request.POST
    _type = inputs.get('type')
    current_user = self.get_current_user(request)
    _id = id_generator()

    location = '/home/schedd'
    structure_dir = location + '/' + current_user.username + '/' + _id
    os.mkdir(structure_dir, 0775)
    os.mkdir(structure_dir + '/segments')
    os.mkdir(structure_dir + '/tasks')
    if _type == 'upload_pdb':
      # User is uploading a pdb
      _name = inputs.get('upload_pdb[name]')
      _desc = inputs.get('upload_pdb[description]')
      _pdb = request.FILES['pdb']
      
      # TODO: Check if the input is correct (no empty input)
      # Save file on disk
      tempfile = open(structure_dir + '/' + _pdb.name, 'w')
      for chunk in _pdb.chunks():
        tempfile.write(chunk)
      tempfile.close();

      structure = Structure.objects.create(user = current_user, directory = _id, file_name = _pdb.name, name = _name, description = _desc, uploaded_at = datetime.datetime.now())
    return HttpResponseRedirect('/structure')

