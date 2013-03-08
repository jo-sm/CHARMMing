#!/bin/sh

# This script is invoked via [py_setup.sh BINDIR PYTHON], so $1 is the bin
# directory and $2 is the python interpreter that we're using.

if [ $# -ne 2 ]; then
	# We must not have a PYTHON interpreter, so we can't do any python programs
	exit 0
fi

# Launch with the appropriate python

sed -e "s@PYTHONEXE@$2@g" < cpinutil.py > $1/cpinutil.py
/bin/chmod +x $1/cpinutil.py

$2 -c "from cpinutils import cpin_data, cpin_utilities" > /dev/null 2>&1

if [ $? -gt 0 ]; then
   echo "Error importing cpinutil.py python modules! cpinutil will not work."
   exit 1
fi

/bin/cp -LR cpinutils $1/

# put a copy of softcore_setup.py in the BINDIR:
sed -e "s@PYTHONEXE@$2@g" < softcore_setup.py > $1/softcore_setup.py
/bin/chmod +x $1/softcore_setup.py
