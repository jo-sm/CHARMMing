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

urlpatterns = patterns('normalmodes.views',
    #(r'^calcatomnumber/$', 'calcAtomNumber'),
    (r'^$', 'normalmodesformdisplay'),

    # Uncomment this for admin:
    (r'^admin/', include(admin.site.urls)),
)
