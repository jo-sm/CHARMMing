#!/bin/csh

# T bath values must match those in replica_exchange.c, rex.py

# LINK CHECK
if ( ! -e ./charmm ) then
 echo "exiting, ./charmm link not found"
 exit(1)
endif

if ( $1 == '' ) then
 echo -n 'Config file name ==> '
 set cfg = $<
else
 set cfg = $1
endif
if ( ! -e $cfg ) then
 echo "exiting, $cfg not found"
 exit(1)
endif

# CREATE T SUBDIRS
echo 'An * indicates a pre-existing subdir'
@ n = 0
foreach t ( `sed '1,5d' $cfg` )
 echo -n "$t"
 if ( ! -d $t ) then
  mkdir $t
  echo -n " "
 else
  echo -n "* "
 endif
 @ n += 1
 if ( 0 == $n % 10 ) echo ''
end
echo ''

