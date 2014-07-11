from django.conf.urls import patterns, url

from assays import views

urlpatterns = patterns('',
    url(r'^$', views.search, name='search'),
    url(r'^search/$', views.search, name='search'),
    url(r'^assays/$', views.assays, name='assays'),
    url(r'^cont/$', views.cont, name='cont')
    )
    
