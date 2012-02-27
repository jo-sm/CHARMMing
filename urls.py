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
from django.template import Context, loader
from django.contrib.auth.views import login, logout
from django.conf.urls.defaults import *
from django.contrib import admin
import charmming_config

admin.autodiscover()

urlpatterns = patterns('',
     (r'^charmming/$', 'account.views.loadFrontPage'),
     (r'^charmming/about/$', 'account.views.loadAboutPage'),
     (r'^charmming/loginuser/$', 'account.views.redirectBasedOnLogin'),
     (r'^charmming/minimize/', include('minimization.urls')),
     (r'^charmming/analysis/', include('analysis.urls')),
     (r'^charmming/solvate/', include('solvation.urls')),
     (r'^charmming/dynamics/', include('dynamics.urls')),
     (r'^charmming/wiki/', include('wiki.urls')),
     (r'^charmming/lessons/', include('lessons.urls')),
     (r'^charmming/normalmodes/', include('normalmodes.urls')),
     (r'^charmming/images/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/mytemplates/images/' % charmming_config.charmming_root}),
     (r'^charmming/tabpage/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/mytemplates/tabs/' % charmming_config.charmming_root}),
     (r'^charmming/js/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/mytemplates/js/' % charmming_config.charmming_root}),
     (r'^charmming/css/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/mytemplates/css/' % charmming_config.charmming_root}),
     (r'^charmming/login/$',  login),
     (r'^charmming/loginfrontpage/$',  'account.views.login'),
     (r'^charmming/accounts/resetpassword/$',  'account.views.resetUserPassword'),
     (r'^charmming/accounts/changepassword/$',  'account.views.changeUserPassword'),
     (r'^charmming/accounts/profile/$',  'account.views.validate'),
     (r'^charmming/accounts/register/teacher/$',  'account.views.registerTeacher'),
     (r'^charmming/accounts/register/student/refreshclasslist/(?P<teacher_name>.*)$',  'account.views.refreshClassDiv'),
     (r'^charmming/accounts/register/student/$',  'account.views.registerStudent'),
     (r'^charmming/accounts/register/$',  'account.views.register'),
     (r'^charmming/accounts/$', 'django.views.generic.simple.direct_to_template', {'template': 'html/frontpage.html'}),
     (r'^charmming/html/skeleton/$',  'account.views.skeleton'),
     (r'^charmming/html/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/mytemplates/html/' % charmming_config.charmming_root}),
     (r'^charmming/accounts/logout/$', logout, {'template_name': 'html/frontpage.html'}),
     #(r'^charmming/contact/$', 'contactform.views.contact'),
     (r'^charmming/fileupload/$', 'structure.views.newupload'),
     (r'^charmming/buildstruct/$', 'structure.views.buildstruct'),
     (r'^charmming/modstruct/$', 'structure.views.modstruct'),
     (r'^charmming/swapstruct/$', 'structure.views.swap'),
     (r'^charmming/energy/$', 'structure.views.energyform'),
     (r'^charmming/fileuploadform/$', 'structure.views.fileuploadform'),
     (r'^charmming/editpdbinfo/(?P<edit_pdb>.*)$', 'structure.views.editPDB'),
     (r'^charmming/calcjobtime/$', 'structure.views.getJobTime'),
     (r'^charmming/viewpdbs/$', 'structure.views.viewpdbs'),
     (r'^charmming/status/$', 'structure.views.viewstatus'),
     (r'^charmming/jmolview/(?P<filename>.*)$', 'structure.views.jmol'),
     (r'^charmming/chemdoodleview/(?P<filename>.*)$', 'structure.views.chemdoodle'),
     (r'^charmming/jmolviewhl/(?P<filename>.*)/(?P<segid>.*)/(?P<resid>.*)/$', 'structure.views.jmolHL'),
     (r'^charmming/chemaxonview/(?P<filename>.*)$', 'structure.views.chemaxon'),
     (r'^charmming/visualize/(?P<filename>.*)$', 'structure.views.visualize'),
     (r'^charmming/jmol/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/jmol/' % charmming_config.charmming_root}),
     # Uncomment the next two lines to enable ChemAxom's Marvin software
     #(r'^charmming/marvin/pdbuploads/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '/var/www/html/charmming/pdbuploads/'}),#way around Marvin's path system
     #(r'^charmming/marvin/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '/var/www/html/charmming/marvin/'}),
     (r'^charmming/pdbuploads/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/pdbuploads/' % charmming_config.charmming_root}),
     (r'^charmming/select/(?P<switch_id>.*)$', 'structure.views.switchpdbs'),
     (r'^charmming/editstructure/(?P<edit_pdb>.*)$', 'structure.views.editPDB'),
     (r'^charmming/viewdynamicoutputcontainer/(?P<jobtype>.*)/$', 'structure.views.viewDynamicOutputContainer'),
     (r'^charmming/viewdynamicoutput/$', 'structure.views.viewDynamicOutput'),
     (r'^charmming/editprotonation/$', 'structure.views.editProtonation'),
     (r'^charmming/greybox/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/mytemplates/GreyBox_v5_5/' % charmming_config.charmming_root}),
     (r'^charmming/viewprocessfiles/$', 'structure.views.viewProcessFiles'),
     (r'^charmming/viewprocessfiles/download/(?P<filename>.*)$', 'structure.views.downloadProcessFiles'),
     (r'^charmming/viewprocessfiles/viewfile/(?P<filename>.*)$', 'structure.views.viewFiles'),
     (r'^charmming/viewprocessfilescontainer/(?P<filename>.*)$', 'structure.views.viewProcessContainer'),
     (r'^charmming/deletefile/$', 'structure.views.deleteFile'),
     (r'^charmming/downloadfiles/downloadtar/$', 'structure.views.downloadTarFile'),
     (r'^charmming/downloadfiles/$', 'structure.views.downloadFilesPage'),
     (r'^charmming/downloadfilescontainer/$', 'structure.views.downloadFilesContainer'),
     (r'^charmming/viewpdbscontainer/$', 'structure.views.viewPDBsContainer'),
     (r'^charmming/reporterror/$', 'structure.views.reportError'),
     (r'^charmming/upload_parse_error/$', 'structure.views.uploadError', {'problem': 'parse_error'}),
     (r'^charmming/upload_overflow/$', 'structure.views.uploadError', {'problem': 'overflow'}),
     (r'^charmming/viewprotores/(?P<segid>.*)/(?P<resid>.*)$', 'structure.views.viewprotores'),
     (r'^charmming/analysis/trajanal/$', 'trajanal.views.trajanalformdisplay'),
     (r'^charmming/analysis/updatetrajanal/$', 'trajanal.views.updatetrajanal'),
     (r'^charmming/analysis/redox/$', 'apbs.views.redoxformdisplay'),

     # edit structure container
     (r'^charmming/editstructurecontainer/$', 'structure.views.editstructurecontainer'),

     # job killing
     (r'^charmming/killjob/(?P<taskid>.*)/$', 'structure.views.killTask'),

     # Uncomment this for admin:
     (r'^charmming/admin/unapproved/', 'account.admin_views.showUnapproved'),
     (r'^charmming/admin/(.*)', admin.site.root),
     (r'^charmming/admin-media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '/usr/share/pyshared/django/contrib/admin/media'}),
)
