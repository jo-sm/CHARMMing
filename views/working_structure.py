from django.views import generic
from django.http import HttpResponse, HttpResponseRedirect
from charmming.models import Structure
from views.mixins import PermissionsMixin
from views.helpers import id_generator
import datetime
import json
import os

class WorkingStructureNewView(PermissionsMixin, generic.TemplateView):
  template_name = 'working_structure/new.html'

  def get(self, request, *args, **kwargs):
    context = self.get_context_data()
    context['structures'] = Structure.objects.filter(user=context['current_user'])[:5]
    return self.render_to_response(context)

  # A user has submitted working structure selections to be generated into a working structure
  def post(self, request, *args, **kwargs):
    name = request.POST.get('name')
    desc = request.POST.get('description')
    # Currently this only supports CHARMM all-atom force fields type working struture creation
    # Will be updated to support GO and BLN models at a later time

    # First, check to see what kind of segments they have (a-good, a-bad)
    # If the segment type is bad, and the topology will be automatically generated
    # Do those processes now
    # We get the SDF from RCSB for each ligand and then send to Dogmans to get the RTF and PRM files
    # We'll use the user's supplied topology files if they've given them

    # Else, we have the topology files and can build the structure in CHARMM
    return self.render_to_response([])
