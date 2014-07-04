	     NPTREX ad hoc interface to CHARMM for REMD
                       venabler@nhlbi.nih.gov
                          November 2007

	An interface to CHARMM for temperature based replica exchange is
described, using a custom Fortran program, and taking advantage of the
dynamics restart capabilities.  Some editing of the scripts is needed
as part of the setup; however, the interface is fairly flexible for this
same reason.  The programs and scripts are presented in the form used to
perform some test REMD simulations on gel phase bilayers of the lipid
DMPC (48 T baths, 500 MD steps, 2000 swaps; 1 ns).  This is an extension
to the AdHocTrex package, modified to handle NPT simulations for the case
where the same reference pressure is used for each T bath replica.  The
implementation is based on the equations given in 

Molecular Simulation of Sucrose Solutions near the Glass
Transition Temperature
Nancy C. Ekdawi-Sever, Paul B. Conrad, and Juan J. de Pablo
J. Phys. Chem. A 105, 734-742 (2001)

	The Fortran program controls the synchronization of the MD
simulations and swaps, and calls various other programs (including CHARMM).
After each MD run, the final potential energy and volume are used to 
evaluate the Metropolis criterion; for a swap, the CHARMM restart files 
are exchanged, and the velocities are scaled based on the T change.  From 
dynamc.doc, the feature used to scale the velocities is:

   ----------------------------------------------------------------
ISCALE      0     This option is to allow the user to scale the velocities
                  by a factor SCALE at the beginning of a restart run.
                  This may be useful in changing the desired temperature.
                  .eq. 0  no scaling done (usual input value)
                  .ne. 0  scale velocities by SCALE.
        WARNING:  Please use this option only when you are changing the
                  temperature of the run.

SCALE      1.     Scale factor for the previous option.
   ----------------------------------------------------------------

The value of SCALE is 1.0 unless a swap occurs; for a swap, the value
computed from the square root of the T ratio is used.  The values are
communicated to CHARMM by a formatted file, which is read via the STREAM
command.  The final energy (?ENER) and volume (?VOLU) from each MD run is 
written to a file using the WRITE TITLE feature.

	REMD is embarassingly parallel; it is ideally suited to Linux
clusters, as there is no communication between processes except via
files between the short MD runs, and each process is self contained. 
Parallel MD, on the other hand, must exchange a fair amount of data
between all processes on every integration step, esp. if particle-mesh
Ewald (PME) summation is being used.  It's simplest to use one processor
for each T bath, although it is certainly possible and reasonable (but
more complicated) to use both processors in a dual-processor compute
node.  However, for a fixed number of nodes, it may often be better to
have twice as many T baths than to double the speed.

	In order to help manage the output, a cleanup script is called
every N swap steps from the Trex program, performing: [1] deletion of
unneeded restart files, [2] gzip compression and relocation of output
files, and [3] merge of the trajectories into larger files and deletion
of the originals (via merge.inp).  Without these steps, one can easily
become overwhelmed by the number of files produced and the disk space used.

	The following listing outlines the steps required to use this
collection of scripts and programs to perform REMD with CHARMM:

 * compile NPTrex.f via e.g.

% set path = ( /usr/local/mpich-gnu/bin $path )
% mpif77 NPTrex.f -o NPTrex.ix

	On a given cluster, the program should only need to be compiled
once in a while, for new machine types or possibly software updates.  The
remaining steps apply to each REMD simulation.

 * Choose the number of T baths and the spacing between baths; some
iterative tests may be needed to get the optimum spacing and to maximize
the T range for a fixed number of baths; integer T values are assumed
 * Decide how many MD steps to perform between swaps, and how many swap
steps for a particular run; again, some initial tests may be needed. 
 * Based on these choices, create the config file (see below); the 
script Tbathconfig.csh can be used to create a suitable config file
 * Create a link named "charmm" in the working directory that points to
a single threaded (non-parallel) version of CHARMM; on NIH Biowulf use
something like

% ln -s /usr/local/charmm/c31b2/pgi-med.one charmm

 * Edit rexstart.inp, rex.inp, and init.str based on the number of MD
steps and protocol, and the actual files (RTF, PARAM, PSF, COOR, etc.)
that will be used; note that at least one coord set must be written to
the .trj file.
 * Edit the merge.inp script, to set the SKIP interval for MERGE
 * Edit the run_rex.csh script if needed, e.g. for some other queuing
system besides PBS, or for the location of the MPICH installation.
 * Run the rex_init.csh script; it sets up the subdirs for each T
bath, and requires the config filename as an argument
 * Submit (PBS) the run_rex.csh script to start the replica exchange;
the example passes 2 args, the number of baths and the config filename

		Monitoring the swap acceptance

	Also included are a couple means of monitoring the swaps and the
acceptance rate:

 * frxn-swap.csh; computes the overall acceptance, and that for each bath
from data in the rexswapN.log files; config file arg, prints a table
 * plt-swap.csh; sample script using gnuplot for frxn swap data
 * trackswap.f; a Fortran program which reads the config file, and makes
bath swap time series files from log data read via stdin (UNIT 5)
 * plt-trk.csh; sample script using gnuplot for bath swap time series

 * smove-temp.py; a Python program, which must be edited (T bath list)
prior to use; this produces data suitable for the 'xmgrace' program via

% smove-temp.py > aswap.dat
% xmgrace aswap.dat


		Restart; extending a run, or after a failure

	Continuation after a successful completion is straightforward; 
simply change the first and last swap step in the config file, as
indicated in the detailed description of the config file below.  Then
submit the run_rex.csh script to the PBS queue as before.
	If a run stops early because of a node failure, network
disruption, or some other problem, it can also be continued, but a
couple of additional steps are needed.  The steps are:

[1] determine the last two file numbers, which should represent the
swap step MD in progress, and the last completed swap step; I use

% ls -str 293

[2] verify the status for all baths; if the last 2 file numbers are
e.g. 694 and 695, the command

% ls -s ???/rex.trj.69[45]

will produce a listing of the last two .trj files for all baths, and
their size (in blocks); for the running step which failed, one or more
of the rex.trj.695 files will have a zero size

[3] in order to preserve the log data, change the name of the log file
to reflect the last completed step, e.g.

% mv rexswap4000.log rexswap694.log

[4] edit the config file, and change the first step number to that of
the step that failed in progress (695 for this example)

[5] resubmit run_rex.csh to a PBS batch queue


		Additional notes on CHARMM input scripts

	The filenames rexstart.inp and rex.inp should probably be
preserved in order to avoid breaking the package; however, there is
a great deal of latitude for the contents of the setup (init.str) and 
the short dynamics run (rex*inp) scripts.  Either NVT (PREF 0.0) or NPT
can be performed; for constant P, either LEAP CPT (pressure.doc) or VV2
(tpcntrl.doc) methods for NPT dynamics should work.  The only part that
should not be changed is the very last section which creates the rex.ene
file for the current swap step. 
	The merge.inp script currently prunes the trajectory files to
1 ps intervals for the lowest T bath, and to 5 ps intervals for all other
baths.  Other than changing the SKIP value for saving coord sets, this
file should not require any further editing.  The example1 subdir has
another script which further merges the files to 1 ns sizes, in order to
faciliate analysis of the results. 

		Configuration file format

	My convention has been to use .cfg as the extension for the
config file, but any extension may be used.  The first 4 lines of the file
specify key parameters; the remaining lines are the integer values for the
T baths.  For the first 4 lines, the first 5 chars of each is a label;
the integer values should start in col 6 or 7.  The data is read based on
order; the first line should be the number of the first swap step for this
particular run, etc.  A simple example:

FIRST 1
LAST  2000
CLEAN 50
PREF  1.0
NBATH 40
300
304
308
312
316
:

The labels indicate:

FIRST	the first swap step; 1 to start, e.g. 4001 to continue
LAST	the last swap step; 1+(LAST-FIRST) CHARMM runs of rex*.inp
CLEAN	frequency for running cleantemp.csh; must be < 90
PREF    external reference pressure for NPT; 0.0 for NVT or vacuum
NBATH	the number of T bath values on the following lines

Note that rextart.inp is only run for the very first swap step, i.e.
step 1 when FIRST is 1.  To continue the above run, the config file
would be changed as indicated below, and run_rex.csh re-submitted to
the PBS queue.

FIRST 2001
LAST  4000
CLEAN 50
PREF  1.0
NBATH 40
300
304
308
312
316
:


		File manifest

NPTrex.f    		master Fortran prog

Tbathconfig.csh		utility to help assign T values

rex_init.csh		set up the T bath subdirs
run_rex.csh		script to submit to PBS queue for REX
cleantemp.csh		periodically called to cleanup
rexdone.csh		optional cleanup after completion 
frxn-swap.csh		monitor swap acceptance
plt-swap.csh		gnuplot example for swap acceptance
trackswap.f		extract bath swap history from logs
plt-trk.csh             gnuplot example for bath swap history
smove-temp.py		build data for moves plot (xmgrace)

example1		subdir illustrating vacuum phase sampling
 dppc10.cfg		 sample config file; 10 T baths
 init.str		 CHARMM initialization; read some files
 rexstart.inp		 the first MD run; write .res and .ene files
 rex.inp		 all other MD runs; reads .res
 merge.inp		 merges .trj files (cleantemp.csh)
 postmerg.inp		 further merge of .trj files, ex post facto

example2		subdir illustrating condensed phase sampling
 dmpc48.cfg		 sample config file; 48 T baths
 init.str		 CHARMM initialization; read some files
 rexstart.inp		 the first MD run; write .res and .ene files
 rex.inp		 all other MD runs; reads .res
 merge.inp		 merges .trj files (cleantemp.csh)

README.txt		this description

