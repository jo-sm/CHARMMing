#!/bin/csh

# query for input data
echo -n "Expt name? "
set enam = $<
echo -n "First swap step? [1] "
set fstp = $<
if ( $fstp == '' ) set fstp = 1
echo -n "Cleanup frequency [50] "
set clen = $<
if ( $clen == '' ) set clen = 50
echo -n "Last swap step? "
set lstp = $<
echo -n "Reference P? [1.0] "
set refp = $<
if ( $refp == '' ) set refp = 1.0
echo -n "No. of baths? "
set ntba = $<
set cfg = $enam$ntba.cfg
echo -n "Lowest T? "
set tlow = $<
echo -n "Initial T incr? "
set iinc = $<
echo -n "Final T incr? "
set finc = $<

# delta delta
set dd = `echo "2 k $finc $iinc - $ntba / p q" | dc`

# initial output
set dt = $iinc
set tb = $tlow
echo "FIRST $fstp" > $cfg
echo "LAST  $lstp" >> $cfg
echo "CLEAN $clen" >> $cfg
echo "PREF  $refp" >> $cfg
echo "NBATH $ntba" >> $cfg
echo $tlow >> $cfg

# computed bath T loop
@ ib = 2
while ( $ib <= $ntba )
 set t = `echo "2 k $tb $dt + p q" | dc`
 echo $t | awk '{print int($1)}' >> $cfg
 set d = $dt
 set tb = $t
 set dt = `echo "2 k $d $dd + p q" | dc`
 @ ib += 1
end

echo "Config file $cfg created"

