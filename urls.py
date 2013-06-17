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
     #(r'^charmming/fileuploadform/$', 'structure.views.fileuploadform'),
     (r'^charmming/editpdbinfo/(?P<edit_pdb>.*)$', 'structure.views.editPDB'),
     (r'^charmming/calcjobtime/$', 'structure.views.getJobTime'),
     (r'^charmming/viewpdbs/$', 'structure.views.viewpdbs'),
     (r'^charmming/status/$', 'structure.views.viewstatus'),
     (r'^charmming/ligand_design/', 'ligdes.views.build_ligand'),
     (r'^charmming/mutation/$', 'mutation.views.selectstructure'),
     (r'^charmming/mutestruct/$', 'mutation.views.mutestructure'),
     (r'^charmming/jmolview/(?P<filename>.*)$', 'structure.views.jmol'),
#     (r'^charmming/chemdoodleview/(?P<filename>.*)$', 'structure.views.chemdoodle'),
     (r'^charmming/jmolviewhl/(?P<filename>.*)/(?P<segid>.*)/(?P<resid>.*)/$', 'structure.views.jmolHL'),
     (r'^charmming/visualize/(?P<filename>.*)$', 'structure.views.visualize'),
     (r'^charmming/glmolview/(?P<filename>.*)$', 'structure.views.glmol'),
     #(r'^charmming/jmol/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '%s/jmol/' % charmming_config.charmming_root}),
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
     (r'^charmming/deletefile/$', 'structure.views.deleteFile'),
     (r'^charmming/downloadfiles/downloadtar/$', 'structure.views.downloadTarFile'),
     (r'^charmming/downloadfiles/$', 'structure.views.downloadFilesPage'),
     (r'^charmming/reporterror/$', 'structure.views.reportError'),
     (r'^charmming/upload_parse_error/$', 'structure.views.uploadError', {'problem': 'parse_error'}),
     (r'^charmming/upload_overflow/$', 'structure.views.uploadError', {'problem': 'overflow'}),
     (r'^charmming/viewprotores/(?P<segid>.*)/(?P<resid>.*)$', 'structure.views.viewprotores'),
     (r'^charmming/analysis/trajanal/$', 'trajanal.views.trajanalformdisplay'),
     (r'^charmming/analysis/updatetrajanal/$', 'trajanal.views.updatetrajanal'),
     (r'^charmming/analysis/redox/$', 'apbs.views.redoxformdisplay'),

     # job killing
     (r'^charmming/killjob/(?P<taskid>.*)/$', 'structure.views.killTask'),

     # Uncomment this for admin:
     (r'^charmming/admin/unapproved/', 'account.admin_views.showUnapproved'),
     (r'^charmming/admin/(.*)', include(admin.site.urls)),
     (r'^charmming/admin-media/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '/usr/share/pyshared/django/contrib/admin/media'}),


     # drug design
     (r'^charmming/dd/(?P<path>.*)$', 'django.views.static.serve', {'document_root': '/var/www/html/charmming/dd/'}),
     (r'^charmming/dd_infrastructure/projectscontainer/$', 'dd_infrastructure.views.viewProjectsContainer'),
     (r'^charmming/dd_infrastructure/viewprojects/$', 'dd_infrastructure.views.viewProjects'),
     (r'^charmming/dd_infrastructure/updateprojectinfo/(?P<projectid>.*)/(?P<action>.*)/$', 'dd_infrastructure.views.viewUpdateProjectInfo'),
     (r'^charmming/dd_infrastructure/projectsaddnew/$', 'dd_infrastructure.views.projectsaddnew'),
     (r'^charmming/dd_infrastructure/select_project/(?P<switch_id>.*)$', 'dd_infrastructure.views.switchProjects'),
     (r'^charmming/dd_infrastructure/projectdetails/(?P<project_id>.*)$', 'dd_infrastructure.views.projectDetails'),
     (r'^charmming/dd_infrastructure/projectavailableconformations/(?P<protein_id>.*)/$', 'dd_infrastructure.views.projectAvailableConformations'),
     (r'^charmming/dd_infrastructure/projectconformations/(?P<project_id>.*)/(?P<addedids>.*)/(?P<removedids>.*)/$', 'dd_infrastructure.views.projectConformations'),
     (r'^charmming/dd_infrastructure/dsfcontainer/', 'dd_infrastructure.views.DSFContainer'),
     (r'^charmming/dd_infrastructure/dsfform/', 'dd_infrastructure.views.DSFFormDisplay'),
     (r'^charmming/dd_infrastructure/updatedsfligands/(?P<selected_set_id>.*)$', 'dd_infrastructure.views.updateDSFLigands'),
     (r'^charmming/dd_infrastructure/updatedsfligandschemspider/(?P<selected_set_id>.*)$', 'dd_infrastructure.views.updateDSFLigandsChemSpider'),
     (r'^charmming/dd_infrastructure/viewjobscontainer/$', 'dd_infrastructure.views.viewJobsContainer'),
     (r'^charmming/dd_infrastructure/viewjobs/$', 'dd_infrastructure.views.viewJobs'),
     (r'^charmming/dd_infrastructure/viewjobdetails/(?P<job_id>.*)$', 'dd_infrastructure.views.viewJobDetails'),
     (r'^charmming/dd_infrastructure/viewjobinfo/(?P<job_id>.*)$', 'dd_infrastructure.views.viewJobInfo'),
     (r'^charmming/dd_infrastructure/viewjobresults/(?P<job_id>.*)$', 'dd_infrastructure.views.viewJobResults'),
     (r'^charmming/dd_infrastructure/viewresultpose/(?P<result_file_id>.*)/(?P<job_id>.*)/$', 'dd_infrastructure.views.viewResultPose'),
     (r'^charmming/dd_infrastructure/downloadfile/(?P<file_id>.*)$', 'dd_infrastructure.views.downloadDDFile'),
     (r'^charmming/dd_infrastructure/viewfile/(?P<file_id>.*)$', 'dd_infrastructure.views.viewDDFile'),
     (r'^charmming/dd_target/targetfileupload/$', 'dd_target.views.newUpload'),
     (r'^charmming/dd_target/targetfileuploadform/$', 'dd_target.views.targetFileUploadForm'),
     (r'^charmming/dd_target/viewtargetdropdown/$', 'dd_target.views.targetDropdown'),
     (r'^charmming/dd_infrastructure/projectd/(?P<projectid>.*)$', 'dd_infrastructure.views.projectD'),
     (r'^charmming/dd_infrastructure/projectdet/$', 'dd_infrastructure.views.projectDet'),
     (r'^charmming/dd_substrate/ligandfileupload/$', 'dd_substrate.views.newLigandUpload'),
     (r'^charmming/dd_substrate/ligandfileuploadform/$', 'dd_substrate.views.ligandFileUploadForm'),
     (r'^charmming/dd_substrate/viewliganddropdown/$', 'dd_substrate.views.ligandDropdown'),
     (r'^charmming/dd_substrate/refreshligands/$', 'dd_substrate.views.RefreshLigands'),
     (r'^charmming/dd_substrate/deleteligandset/(?P<set_id>.*)$', 'dd_substrate.views.deleteLigandSet'),
     (r'^charmming/dd_substrate/viewligandsetscontainer/$', 'dd_substrate.views.viewLigandSetsContainer'),
     (r'^charmming/dd_substrate/viewligandsets/$', 'dd_substrate.views.viewLigandSets'),
     (r'^charmming/dd_substrate/updateligandsetinfo/(?P<ligandsetid>.*)/(?P<action>.*)/$', 'dd_substrate.views.viewUpdateLigandSetInfo'),
     (r'^charmming/dd_substrate/ligandsetsaddnew/$', 'dd_substrate.views.ligandSetsAddNew'),
     (r'^charmming/dd_substrate/ligandsetdetails/(?P<ligandset_id>.*)$', 'dd_substrate.views.ligandSetDetails'),
     (r'^charmming/dd_substrate/setavailableligands/(?P<selected_set_id>.*)$', 'dd_substrate.views.updateAvailableSetLigands'),
     (r'^charmming/dd_substrate/setligands/(?P<ligandset_id>.*)/(?P<addedids>.*)/(?P<removedids>.*)/$', 'dd_substrate.views.setLigands'),
     (r'^charmming/dd_substrate/viewligands/$', 'dd_substrate.views.viewligands'),
     (r'^charmming/dd_substrate/deleteligand/$', 'dd_substrate.views.deleteligand'),
     (r'^charmming/dd_substrate/viewligandjmol/(?P<ligand_file_id>.*)/$', 'dd_substrate.views.viewLigandJmol'),
     (r'^charmming/dd_target/viewtargetscontainer/$', 'dd_target.views.viewTargetsContainer'),
     (r'^charmming/dd_target/viewtargets/$', 'dd_target.views.viewTargets'),

)
