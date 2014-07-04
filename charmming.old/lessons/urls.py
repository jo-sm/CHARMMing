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
# the lesson_num_lis from lesson_config is needed
from lesson_config import *

urlpatterns = patterns('',
    (r'^download/(?P<lessonfolder>.*)/(?P<filename>.*)$', 'lessons.aux.downloadLessonFiles'),
    (r'^lessonstatus/$', 'lessons.views.displayLessonStatus'),
    # we need a line with admin/ for the lessons generator 
)
# add all lessons' urls
for num in lesson_num_lis:
    urlpatterns += patterns('',(r'^lesson'+num+'/', include('lesson'+num+'.urls')),)

#urlpatterns += patterns('',
#    (r'^lesson1/', include('lesson1.urls')),
#    (r'^lesson2/', include('lesson2.urls')),
#    (r'^lesson3/', include('lesson3.urls')),
#    (r'^lesson4/', include('lesson4.urls')),
#    (r'^lesson96/', include('lesson96.urls')),
#    (r'^lesson95/', include('lesson95.urls')),
#    (r'^lesson94/', include('lesson94.urls')),
#)    
