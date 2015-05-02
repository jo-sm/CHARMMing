#!/bin/bash

# Kill any processes that weren't killed previously
# Kill livereload
# Kill django

# Build JSX to JS
grunt &

# Watch JS for changes
grunt watch:js &
# Watch css for livereload
grunt watch:css & 

# django local server
python manage.py runserver
