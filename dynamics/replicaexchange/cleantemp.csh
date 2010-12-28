#!/bin/csh
# utility script, called from NPTrex prog; cleanup files

# arg1 is T value, arg2 is swap step counter, arg3 is interval
# arg2 must be an exact multiple of arg3
# arg3 must be less than 90; Fortran limit on the no. of open files
# arg4 is the pid; allows special processing for pid=0

# synch nfs mounted disks
sync

# preserve random state for possible restart
if ( $4 == 0 ) then
 if ( ! -e SavRand ) mkdir SavRand
 cp randstate.bdt SavRand/randstate$2.bdt
endif

pushd $1 >& /dev/null
if ( ! -e Out ) mkdir Out
if ( ! -e Res ) mkdir Res
# list of unprocessed files
ls -1tr rex.res.* > r.t
set ff = `head -1 r.t | cut -d. -f3`
set lf = `tail -1 r.t | cut -d. -f3`
# preserve last restart file and rex.str of each merge group
gzip < rex.res.$lf > Res/rex.res.$lf.gz
cp rex.str Res/rex$lf.str

@ lr = $lf - 1
@ i = $ff
# wipe restart
while ( $i <= $lr )
 rm rex.res.$i
 @ i += 1
end

@ i = $ff + 1
date > Out/mrg$lf.out
# merge and remove output
while ( $i <= $lf )
 cat rex$i.out >> Out/mrg$lf.out
 rm rex$i.out
 @ i += 1
end

popd >& /dev/null

# merge .trj files; orig removed via CLOSE w. DISP DELETE
./charmm < merge.inp TACT:$1 END:$2 I:$3 PID:$4 >> $1/Out/mrg$lf.out

# swap step output merged, .trj merge output appended; compressed
gzip -f $1/Out/mrg$lf.out

