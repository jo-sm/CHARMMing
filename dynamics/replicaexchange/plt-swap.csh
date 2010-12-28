#!/bin/csh
# use gnuplot to show swap data
# build plot file in-line

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

echo -n 'Gen data? [n] '
set a = $<
if ( $a == y ) then
 frxn-swap.csh $cfg > swap.txt 
 sed -e '1d' -e '$d' swap.txt > swap.dat
 sed '1d' swap.dat > up.t
 sed '$d' swap.dat | paste up.t - | awk '{print $1-$3}' > dff.dat
 rm up.t
endif

set ns = `head -1 swap.txt | awk '{print $5}'`
set as = `tail -1 swap.txt | awk '{print $4}'`

cat > g.t << FIN
set styl data linesp
set tit "DMPC gel, model III NPTrex, $ns steps; <a(swap)>=$as"
set xl 'Bath No.'
set yl 'Temp (deg K)'
set y2l 'Frxn Swap Acceptance'
set yti nomirror
set y2ti
set y2r [0.0:0.5]
plot [] [275:500] 'swap.dat' u 1 t 'bath T', \
     'swap.dat' u 2 axes x1y2 t 'a(swap)'
pause -1
set te post color solid 'Palatino' 16
set out 'btswp.ps'
replot
set out
se te x11

plot [] [2:7] 'dff.dat' u 1 t 'delta T', \
     'swap.dat' u 2 axes x1y2 t 'a(swap)'
pause -1
set te post color solid 'Palatino' 16
set out 'dtswp.ps'
replot
set out
se te x11
f(x) = a + b*x
fit f(x) 'swap.dat' via a, b
set styl data points
set tit "DMPC gel, model III NPTrex, Linear Fit, a(swap) vs. T"
set xl 'Temp (deg K)'
set yl 'Frxn Swap Acceptance'
se y2l
set yti mirror
set noy2ti
plot [275:] 'swap.dat', f(x) t 'fit'
pause -1
set te post color solid 'Palatino' 16
set out 'fitswp.ps'
replot
set out

FIN
gnuplot g.t
mpage -c -4 -m0 -M0 btswp.ps dtswp.ps fitswp.ps > bathswapII.ps
ps2pdf bathswapII.ps
rm g.t btswp.ps dtswp.ps fitswp.ps

