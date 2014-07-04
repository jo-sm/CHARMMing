#!/bin/bash -x 

mv schedd.log.4.gz schedd.log.5.gz >& /dev/null
mv schedd.log.3.gz schedd.log.4.gz >& /dev/null
mv schedd.log.2.gz schedd.log.3.gz >& /dev/null
mv schedd.log.1.gz schedd.log.2.gz >& /dev/null
mv schedd.log schedd.log.1
killall -HUP schedd.py
gzip schedd.log.1
