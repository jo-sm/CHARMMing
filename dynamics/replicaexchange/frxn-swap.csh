#!/bin/csh

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

cat `ls -1tr rexswap*.log` > tmp.log
# get last swap step from swap.log
@ ns = `tail -10 tmp.log | grep switched | tail -1 | awk '{print $8}'`
@ hs = $ns / 2

@ nt = 0
echo " T  a(swap), N = $ns"
set b = ( `sed '1,5d' $cfg` )
set nb = ` sed '1,5d' $cfg | wc -l`
foreach t ( $b )
 echo -n "$t "
 set k = `grep " $t " tmp.log | wc -l`
 if ( $t == $b[1] || $t == $b[$nb] ) then
  echo "3 k $k $hs / p q" | dc
 else
  echo "3 k $k $ns / p q" | dc
 endif
 @ nt += 1
end

@ ht = $nt / 2
@ m1 = $ht - 1
# COMPUTE OVERALL AND BATH ACCEPTANCE RATES
set s = `grep switched tmp.log | wc -l`
echo -n "Overall a(swap) = "
echo "3 k $s $hs $ht * $hs $m1 * + / p q" | dc
rm tmp.log

