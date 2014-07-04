#!/bin/csh

if ( $1 == '' ) exit(1)
set p1 = "plot '$1/trkswp.dat' t '$1'"
set p2 = "plot '$1/trkswp.dat' u 1:3 t '$1'"

if ( $2 != '' ) then
 foreach t ( $argv[2-] )
  set p1 = "$p1, '$t/trkswp.dat' t '$t'"
  set p2 = "$p2, '$t/trkswp.dat' u 1:3 t '$t'"
 end
endif

set c = `echo $cwd | cut -d/ -f5-`

cat > g.t << FIN
se da st li
se tit "$c"
#set nokey
se xl 'Swap step'
se yl 'Bath History'
#se yr [280:560]
$p1
pause -1
set te post color solid "Palatino" 18
set out "trk1.ps"
replot
se out
se te x11
se yl 'Bath Identity'
$p2
pause -1
set te post color solid "Palatino" 18
set out "trk2.ps"
replot
se out

FIN
gnuplot g.t
#mpage -2 -c -m0 -M0 trk?.ps > trkswp.ps
#ps2pdf trkswp.ps
#rm trk?.ps

