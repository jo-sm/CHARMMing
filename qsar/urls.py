# -*- coding: utf-8 -*-
from django.conf.urls.defaults import patterns, url

urlpatterns = patterns('qsar.views',
    url(r'^$', 'upload', name='upload'),
    url(r'^upload/$', 'upload', name='upload'),
    url(r'^(?P<filename>.+)/property/$', 'property', name='property'),
    url(r'^property/$', 'property', name='property'),
    url(r'^(?P<activity_property>.+)/train/$', 'train', name='train'),
    url(r'^train/$', 'train', name='train'),
    url(r'^predict/$', 'predict', name='predict'),
    url(r'^download/$', 'download', name='download'),
    url(r'^viewjobs/$', 'viewJobs', name='viewJobs'),
    url(r'^newmodel/$', 'newModel', name='newModel'),
    url(r'^viewmodels/$', 'viewModels', name='viewModels'),
    url(r'^newpredict/(?P<model_id>.*)$', 'newPredict', name='newPredict'),
    url(r'^downloadqsarresults/(?P<qsar_job_id>.*)$', 'downloadQSARResults', name='downloadQSARResults'),
)

