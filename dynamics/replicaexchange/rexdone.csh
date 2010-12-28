#!/bin/csh
# final .out, .res cleanup

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

foreach t ( `sed '1,5d' $cfg` )
 pushd $t >& /dev/null
 ls -1tr rex.res.* > r.t
 set ff = `head -1 r.t | cut -d. -f3`
 set lf = `tail -1 r.t | cut -d. -f3`
 @ i = $ff
# compress output, move to subdir; wipe restart
 while ( $i < $lf )
  gzip rex$i.out
  rm rex.res.$i
  mv rex$i.out.gz Out
  @ i += 1
 end
 gzip rex.res.$lf merge.out
 popd >& /dev/null
end

