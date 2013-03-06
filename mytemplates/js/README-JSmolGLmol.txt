Important! This file contains the files that need to be replaced.
While a simple diff may suffice, this is probably for the best in case 
something goes wrong.

Add the following folders:
/charmming/mytemplates/js/jsmol - JSmol libraries, ~10 MB. Anything ending in .htm or .php or .xml in that folder can probably be removed, they're just test pages for JSmol. I can't do it now because it got chown'd to schedd.
/charmming/mytemplates/js/glmol - GLmol libraries, ~500 KB.

Add the following files:
/charmming/mytemplates/html/glmol.html - GLmol page.

Replace the following existing files:
/charmming/mytemplates/html/jmol.html - Has been replaced with JSmol visualizer
/charmming/mytemplates/html/visualize.html - Has JSmol and GLmol options, also opens a new window for each of them. 
/charmming/mytemplates/html/ddviewresultposejmol.html - Replaces drug design viewer with JSmol viewer. Since it opens in a new window there should be no issues; JSmol and Jmol scripting is identical.
/charmming/structure/views.py - Adds function "glmol" for rendering glmol page properly, as well as flags for rendering proteins/good/bad as cartoons/ball-and-stick respectively. Please check file open code - Python can only access files on local disk; it does not go through Apache's server structure. If path on actual server is NOT /var/www/charmming, the code will break.
/charmming/urls.py - Adds glmol.html

Original files are backed up as:
/charmming/mytemplates/html/jmoloriginal.html
/charmming/mytemplates/html/originalvisualize.html
/charmming/mytemplates/html/ddviewresultposeoriginal.html - Please verify this one against Yuri's, he may have a different setup in the SVN at this point.
/charmming/urlsoriginal.py
/charmming/structure/viewsoriginal.py

Replacements in default libraries:
JSmol - Nothing was changed in the base libraries; however, if you wish to remove the "File" options, the following files will see a change:
/charmming/mytemplates/js/jsmol/j2s/org/jmol/popup/MainPopupResourceBundle.js
/charmming/mytemplates/js/jsmol/j2s/core/coremenu.z.js

The operation is detailed in the file /charmming/mytemplates/js/jsmol/j2s/core/README.txt. Please read that to understand what is being done to the menu in question.

GLmol - Changes are more substantial, but are limited to a single file, /charmming/mytemplates/js/glmol/js/GLmol.js

This file is the basic GLmol library file. The changes are to the functions removeSolvents, getNonbonded, and parsePDB2, all of which are required for GLmol to properly parse our PDB files.
