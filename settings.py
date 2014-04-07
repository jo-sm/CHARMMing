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
# Django settings for charmming.
# lesson_config has a list called 'lesson_num_lis'
# which contains all lessons created

import os
from lesson_config import * 

DEBUG = True
TEMPLATE_DEBUG = True

ADMINS = (
    ('Tim Miller', 'btmiller@helix.nih.gov'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
                  'ENGINE': 'django.db.backends.mysql',
                  'HOST': '/var/run/mysqld/mysqld.sock',
                  'NAME': 'charmming',
                  'USER': 'charmming',
                  'PASSWORD': 'qwerty',
               }
}

#DATABASE_ENGINE = 'mysql' # 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'ado_mssql'.
DATABASE_NAME = 'charmming'             # Or path to database file if using sqlite3.
DATABASE_USER = 'charmming'         # Not used with sqlite3.
DATABASE_PASSWORD = 'qwerty'         # Not used with sqlite3.
DATABASE_HOST = ''             # Set to empty string for localhost. Not used with sqlite3.
DATABASE_PORT = ''             # Set to empty string for default. Not used with sqlite3.
TIME_ZONE = 'EST5EDT'

# Language code for this installation. All choices can be found here:
# http://www.w3.org/TR/REC-html40/struct/dirlang.html#langcodes
# http://blogs.law.harvard.edu/tech/stories/storyReader$15
LANGUAGE_CODE = 'en-us'

SITE_ID = 1

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = False

# Absolute path to the directory that holds media.
# Example: "/home/media/media.lawrence.com/"
MEDIA_ROOT = '/home/schedd'

# URL that handles the media served from MEDIA_ROOT.
# Example: "http://media.lawrence.com"
MEDIA_URL = '/charmming/schedd/'

# static URL config.
STATIC_ROOT = os.path.join(os.path.dirname(__file__), 'static/')
STATIC_URL = '/charmming/static/'

# URL prefix for admin media -- CSS, JavaScript and images. Make sure to use a
# trailing slash.
# Examples: "http://foo.com/media/", "/media/".
ADMIN_MEDIA_PREFIX = '/charmming/static/admin/'

# Make this unique, and don't share it with anybody.
SECRET_KEY = 'f0!9z)dca4v=y#c0-ojs@gzu@hf&x53drd5txbv!rdx=lol#-r'

AUTH_PROFILE_MODULE = "UserProfile"

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    # Uncomment the next line for simple clickjacking protection:
    # 'django.middleware.clickjacking.XFrameOptionsMiddleware',
    # ToDo: add CSRF middleware.
    #'django.middleware.csrf.CsrfViewMiddleware',
    # Uncomment the next line for simple clickjacking protection:
    # 'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

TEMPLATE_CONTEXT_PROCESSORS = (
    "django.contrib.auth.context_processors.auth",
    "django.core.context_processors.debug",
    "django.core.context_processors.i18n",
    "django.core.context_processors.media",
    "django.core.context_processors.static",
    "django.core.context_processors.tz",
    "django.contrib.messages.context_processors.messages"
)

# List of callables that know how to import templates from various sources.
#TEMPLATE_LOADERS = (
#     'django.template.loaders.filesystem.Loader',
#     'django.template.loaders.app_directories.Loader',
#    'django.template.loaders.filesystem.load_template_source',
#    'django.template.loaders.app_directories.load_template_source',
#     'django.template.loaders.eggs.load_template_source',
#)

#MIDDLEWARE_CLASSES = (
#    'django.middleware.common.CommonMiddleware',
#    'django.contrib.sessions.middleware.SessionMiddleware',
#    'django.contrib.auth.middleware.AuthenticationMiddleware',
#    'django.middleware.doc.XViewMiddleware',
#)

ROOT_URLCONF = 'urls'

TEMPLATE_DIRS = (
    "/var/www/charmming/mytemplates",
    # Put strings here, like "/home/html/django_templates" or "C:/www/django/templates".
    # Always use forward slashes, even on Windows.
    # Don't forget to use absolute paths, not relative paths.
)

INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.messages',
    'django.contrib.sessions',
    'django.contrib.admin',
    'django.contrib.sites',
    'django.contrib.staticfiles',
    'account',
    'analysis',
    'apbs',
    'dd_analysis',
    'dd_infrastructure',
    'dd_substrate',
    'dd_target',
    'dynamics',
    'lessons',
    'ligdes',
    'minimization',
    'mutation',
    'normalmodes',
    'selection',
    'solvation',
    'statistics',
    'structure',
    'solvation',
    'toppar',
    'trajanal',
    'qsar',
    'api',
)

# all lessons in lesson_num_lis(lesson_config.py) will be put in extra_apps and then added into INSTALLED_APPS
extra_apps = ()
for num in lesson_num_lis:
    extra_apps = extra_apps + ('lesson'+num,)

INSTALLED_APPS = INSTALLED_APPS + extra_apps 

#print INSTALLED_APPS
