      program nptrex
C temperature replica exchange; venabler AT nhlbi*nih*gov
C constant P version in the AdHocTrex series; fixed PREF
C compile with MPICH mpif77 (or mpif90)

      implicit none
      include 'mpif.h'
      integer nbath, mxba, istrt,istop, ierr, nclen, i,j,k, nran
      parameter (mxba=256,nran=24)
      integer ibath(0:mxba-1), pid, nprocs, is, kt,ks, lastnb
      real bt1,bt2, delta,tscal, pex,ra, random, rt(nran+1), pfac
      real emd(0:mxba-1), vmd(0:mxba-1), tb,tn, pref, dbdu, dbdv
C gasj     gas constant in joules
C kcalj    1 kcal = 4184. joules
C atmosp   1 ATMOSP = 1.4584007 10^-5 Kcal/Mol/Angs.^3; consta.fcm
      real gasj, kcalj, atmosp
      parameter (atmosp=1.4584007D-05)
      parameter (gasj=8.31451, kcalj=4184.)
      character*80  str
      character*12  tmp,stb,sst,dbg
      character*512 cmd
      external random, lastnb

C MPI initialization and setup
      call MPI_INIT(ierr)
      call MPI_COMM_RANK( MPI_COMM_WORLD, pid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr )

C read the Trex config file from process 0
      if (pid.eq.0) then
C cannot check for j.eq.1 as there are (apparently) other args
        j = iargc()
        if (j.eq.0) then
          write(6,*) 'Missing arg'
          stop
        else
          call getarg(1,str)
          if (j.gt.1) call getarg(2,dbg)
        endif
        open(1,file=str,form='formatted',status='old')
        read(1,'(A5,I10)') tmp,istrt
        read(1,'(A5,I10)') tmp,istop
        read(1,'(A5,I10)') tmp,nclen
        read(1,'(A5,F10.2)') tmp,pref
        read(1,'(A5,I10)') tmp,nbath
        if (nbath.gt.mxba) stop 'MXBA limit exceeded'
        do i=0,nbath-1
          read(1,'(I10)',end=750,err=750) ibath(i)
        enddo
        close(1)
        write(6,701) istrt,istop,nclen,pref,nbath,nprocs
        if (nbath.ne.nprocs) stop 'no. procs, no. baths mismatch'
C check for random number state file; print warning if not found
        if (istrt.ne.1) then
          open(17,file='randstate.bdt',form='unformatted',
     *       err=120,status='old')
          read(17) (rt(j),j=1,nran+1)
          close(17)
          ra=random(rt,-nran)
          write(6,721)
          goto 150
 120      write(6,722)
        endif
 150    j = istop
        call leftint(j,tmp)
        str= 'rexswap'//tmp(1:j)//'.log'
        open(12,file=str,form='formatted',status='new')
C DEBUG LOG FILE
        if (dbg.EQ.'debug') then
          str= 'swapstep'//tmp(1:j)//'.log'
          open(22,file=str,form='formatted',status='new')
        endif
      endif
 701  format('START=',I6,3X,'STOP=',I6,3X,'CLEAN=',I3,3X,
     *  'PREF=',F7.2,3X, 'NBATH=',I3,3X,'NPROCS=',I3)
 721  format('Random number state read from file')
 722  format('Problem reading random number state file; resetting')

C broadcast the data read by pid 0; must exec on each pid
      call MPI_BCAST(istrt,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(istop,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nclen,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nbath,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibath,nbath,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

C reporting in
      write(6,702) pid,ibath(pid)
 702  format('Process ',I3,5X,'Temperature = ',I4)
      call flush()

C main loop over swap steps
      do is=istrt,istop

C build command line, exec via system
      ks = is
      call leftint(ks,sst)
      kt = ibath(pid)
      call leftint(kt,stb)
      cmd='./charmm tact:'
      j = lastnb(cmd)
      cmd=cmd(1:j)//stb(1:kt)//' i:'
      j = lastnb(cmd)
      if (is.eq.1) then
        cmd=cmd(1:j)//sst(1:ks)//' < rexstart.inp > '
      else
        cmd=cmd(1:j)//sst(1:ks)//' < rex.inp > '
      endif
      j = lastnb(cmd)
      cmd=cmd(1:j)//stb(1:kt)//'/rex'
      j = lastnb(cmd)
      cmd=cmd(1:j)//sst(1:ks)//'.out; sync'
      call system(cmd,ierr)
C wait for all processes to finish CHARMM run
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C the next bits are for swap evaluation
      if (pid.eq.0) then
        tscal=1.0
        do i=0,nbath-1
          kt= ibath(i)
          call leftint(kt,stb)
C write stream files for restart pointer, assuming no swaps
          str=stb(1:kt)//'/rex.str'
          open(2,file=str,form='formatted')
          write(2,'(A)') '* restart pointer stream'
          write(2,'(A)') '*'
          str='set oldrest '//stb(1:kt)//'/rex.res.@J'
          write(2,'(A)') str
          write(2,800) tscal
          write(2,'(A)') 'return'
          close(2)
C fill energy array for eval
          str=stb(1:kt)//'/rex.ene'
          open(3,file=str,form='formatted',status='old')
          if (pref.gt.0.0) then
            read(3,801) emd(i),vmd(i)
          else
            read(3,811) emd(i)
          endif
          close(3)
        enddo
 800    format('set sfactor ',F12.6)
 801    format(1X,F20.6,1X,F20.6)
 811    format(1X,F20.6)

C randomly pick the starting bath, 0 or 1
        ra=random(rt,nran)
        k=0
        if (ra.gt.0.5) k=1
 
C pressure * factor; factor converts PV to kcal
        pfac = pref * atmosp

        do i=k,nbath-2,2
          tb = real(ibath(i))
          tn = real(ibath(i+1))
C compute Upot term
          bt1 = 1./(gasj*tb)
          bt2 = 1./(gasj*tn)
          dbdu = kcalj*(bt1-bt2)*(emd(i+1)-emd(i))
C compute PV term
          if (pref.gt.0.0) then
            dbdv = kcalj*(bt1-bt2)*pfac*(vmd(i+1)-vmd(i))
          else
            dbdv = 0.
          endif
          delta = dbdu + dbdv
          pex = exp(-delta)
          ra=random(rt,nran)
          if (pex.gt.ra) then
            kt=ibath(i)
            call leftint(kt,stb)
            ks=ibath(i+1)
            call leftint(ks,sst)
C downward swap
            tscal=sqrt(tb/tn)
            str=stb(1:kt)//'/rex.str'
            open(2,file=str,form='formatted')
            write(2,'(A)') '* restart pointer stream'
            write(2,'(A)') '*'
            str='set oldrest '//sst(1:ks)//'/rex.res.@J'
            write(2,'(A)') str
            write(2,800) tscal
            write(2,'(A)') 'return'
            close(2)
C upward swap
            tscal=sqrt(tn/tb)
            str=sst(1:ks)//'/rex.str'
            open(2,file=str,form='formatted')
            write(2,'(A)') '* restart pointer stream'
            write(2,'(A)') '*'
            str='set oldrest '//stb(1:kt)//'/rex.res.@J'
            write(2,'(A)') str
            write(2,800) tscal
            write(2,'(A)') 'return'
            close(2)
            if (pref.gt.0.0) then
              write(12,802) (ibath(j),j=i,i+1), is,
     *          emd(i),emd(i+1),pex,
     *          vmd(i),vmd(i+1),dbdu,dbdv
            else
              write(12,812) (ibath(j),j=i,i+1), is,
     *          emd(i),emd(i+1),pex
            endif
            call flush(12)
          endif
C FOR DEBUGGING; PRINT ALL SWAP DATA
          if (dbg.EQ.'debug') then
            if (pref.gt.0.0) then
              write(22,822) is, (ibath(j),j=i,i+1), 
     *          emd(i),emd(i+1),pex,
     *          vmd(i),vmd(i+1),dbdu,dbdv
            else
              write(22,832) is, (ibath(j),j=i,i+1), 
     *          emd(i),emd(i+1),pex
            endif
            call flush(22)
          endif
        enddo
C save random number generator state on each pass
        open(17,file='randstate.bdt',form='unformatted')
        write(17) (rt(j),j=1,nran+1)
        close(17)
C end pid 0 swap processing
      endif
 802  format('the following baths were switched: ',2I5,I9,/,
     *  'e1,e2,prob: ',3F12.3,/,'v1,v2,dBdU,PdBdV: ',2F12.1,2E13.5)
 812  format('the following baths were switched: ',2I5,I9,/,
     *  'e1,e2,prob: ',3F12.3)
 822  format(I9,2I4,3F12.3,2F12.1,2E13.5)
 832  format(I9,2I4,3F12.3)

C cleanup; merge 
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      if (mod(is,nclen).eq.0) then
        write(cmd,803) ibath(pid), is, nclen, pid
        call system(cmd,ierr)
      endif
 803  format('./cleantemp.csh ',I6,1X,I8,1X,I4,1X,I4,' > /dev/null')
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)

C end main loop
      enddo

      call MPI_FINALIZE(ierr)

      stop
 750  write(6,*) 'Error reading T bath values'
      stop
      end

C string utility routines

      integer function lastnb(s)
C last non-blank char in a string
      character*(*) s
      integer n,i
      n=len(s)
      lastnb=0
      do i=n,1,-1
        if (s(i:i).ne.' ') then
          lastnb=i
          return
        endif
      enddo
      return
      end

      subroutine leftint(n,s)
C left justify an integer in a string, report no. of digits
C n inout, s out; both var args are modified
      integer i,j,n
      character*(*) s
      if (n.lt.10) then
        write(s,'(I1)') n
        n = 1
      else if (n.lt.100) then
        write(s,'(I2)') n
        n = 2
      else if (n.lt.1000) then
        write(s,'(I3)') n
        n = 3
      else if (n.lt.10000) then
        write(s,'(I4)') n
        n = 4
      else if (n.lt.100000) then
        write(s,'(I5)') n
        n = 5
      else if (n.lt.1000000) then
        write(s,'(I6)') n
        n = 6
      else if (n.lt.10000000) then
        write(s,'(I7)') n
        n = 7
      else if (n.lt.100000000) then
        write(s,'(I8)') n
        n = 8
      endif
      return
      end

C routines from NETLIB for uniform random numbers

      real function random (t, n)
c
c this random number generator is portable amoung a wide variety of
c computers.  it generates a random number between 0.0 and 1.0 accord-
c ing to the algorithm presented by bays and durham (toms, 2, 59,
c 1976).  the motivation for using this scheme, which resembles the
c maclaren-marsaglia method, is to greatly increase the period of the
c random sequence.  if the period of the basic generator (rand) is p,
c then the expected mean period of the sequence generated by random is
c given by   new mean p = sqrt (pi*factorial(n)/(8*p)),
c where factorial(n) must be much greater than p in this asymptotic
c formula.  generally, n should be 16 to maybe 32.
c
c             input argument --
c n      iabs(n) is the number of random numbers in an auxiliary table.
c        note though that iabs(n)+1 is the number of items in array t.
c        if n is positive and differs from its value in the previous
c        invocation, then the table is initialized for the new value of
c        n.  if n is negative, iabs(n) is the number of items in an
c        auxiliary table, but the tables are now assumed already to
c        be initialized.  this option enables the user to save the
c        table t at the end of a long computer run and to restart with
c        the same sequence.  normally, random would be called at most
c        once with negative n.  subsequent invocations would have n
c        positive and of the correct magnitude.
c
c             input and output argument  --
c t      an array of iabs(n)+1 random numbers from a previous invocation
c        of random.  whenever n is positive and differs from the old
c        n, the table is initialized.  the first iabs(n) numbers are the
c        table discussed in the reference, and the n+1 -st value is y.
c        this array may be saved in order to restart a sequence.
c
c             output value --
c random a random number between 0.0 and 1.0.
c
      dimension t(n)
      external rand
      data nold, floatn / -1, -1.0 /
c
      if (n.eq.nold) go to 20
c
      nold = iabs(n)
      floatn = nold
      if (n.lt.0) dummy = rand (t(nold+1))
      if (n.lt.0) go to 20
c
      do 10 i=1,nold
        t(i) = rand (0.)
 10   continue
      t(nold+1) = rand (0.)
c
 20   j = t(nold+1)*floatn + 1.
      t(nold+1) = t(j)
      random = t(j)
      t(j) = rand (0.)
c
      return
      end

      function rand (r)
c
c      this pseudo-random number generator is portable amoung a wide
c variety of computers.  rand(r) undoubtedly is not as good as many
c readily available installation dependent versions, and so this
c routine is not recommended for widespread usage.  its redeeming
c feature is that the exact same random numbers (to within final round-
c off error) can be generated from machine to machine.  thus, programs
c that make use of random numbers can be easily transported to and
c checked in a new environment.
c      the random numbers are generated by the linear congruential
c method described, e.g., by knuth in seminumerical methods (p.9),
c addison-wesley, 1969.  given the i-th number of a pseudo-random
c sequence, the i+1 -st number is generated from
c             x(i+1) = (a*x(i) + c) mod m,
c where here m = 2**22 = 4194304, c = 1731 and several suitable values
c of the multiplier a are discussed below.  both the multiplier a and
c random number x are represented in double precision as two 11-bit
c words.  the constants are chosen so that the period is the maximum
c possible, 4194304.
c      in order that the same numbers be generated from machine to
c machine, it is necessary that 23-bit integers be reducible modulo
c 2**11 exactly, that 23-bit integers be added exactly, and that 11-bit
c integers be multiplied exactly.  furthermore, if the restart option
c is used (where r is between 0 and 1), then the product r*2**22 =
c r*4194304 must be correct to the nearest integer.
c      the first four random numbers should be .0004127026,
c .6750836372, .1614754200, and .9086198807.  the tenth random number
c is .5527787209, and the hundredth is .3600893021 .  the thousandth
c number should be .2176990509 .
c      in order to generate several effectively independent sequences
c with the same generator, it is necessary to know the random number
c for several widely spaced calls.  the i-th random number times 2**22,
c where i=k*p/8 and p is the period of the sequence (p = 2**22), is
c still of the form l*p/8.  in particular we find the i-th random
c number multiplied by 2**22 is given by
c i   =  0  1*p/8  2*p/8  3*p/8  4*p/8  5*p/8  6*p/8  7*p/8  8*p/8
c rand=  0  5*p/8  2*p/8  7*p/8  4*p/8  1*p/8  6*p/8  3*p/8  0
c thus the 4*p/8 = 2097152 random number is 2097152/2**22.
c      several multipliers have been subjected to the spectral test
c (see knuth, p. 82).  four suitable multipliers roughly in order of
c goodness according to the spectral test are
c    3146757 = 1536*2048 + 1029 = 2**21 + 2**20 + 2**10 + 5
c    2098181 = 1024*2048 + 1029 = 2**21 + 2**10 + 5
c    3146245 = 1536*2048 +  517 = 2**21 + 2**20 + 2**9 + 5
c    2776669 = 1355*2048 + 1629 = 5**9 + 7**7 + 1
c
c      in the table below log10(nu(i)) gives roughly the number of
c random decimal digits in the random numbers considered i at a time.
c c is the primary measure of goodness.  in both cases bigger is better.
c
c                   log10 nu(i)              c(i)
c       a       i=2  i=3  i=4  i=5    i=2  i=3  i=4  i=5
c
c    3146757    3.3  2.0  1.6  1.3    3.1  1.3  4.6  2.6
c    2098181    3.3  2.0  1.6  1.2    3.2  1.3  4.6  1.7
c    3146245    3.3  2.2  1.5  1.1    3.2  4.2  1.1  0.4
c    2776669    3.3  2.1  1.6  1.3    2.5  2.0  1.9  2.6
c   best
c    possible   3.3  2.3  1.7  1.4    3.6  5.9  9.7  14.9
c
c             input argument --
c r      if r=0., the next random number of the sequence is generated.
c        if r.lt.0., the last generated number will be returned for
c          possible use in a restart procedure.
c        if r.gt.0., the sequence of random numbers will start with the
c          seed r mod 1.  this seed is also returned as the value of
c          rand provided the arithmetic is done exactly.
c
c             output value --
c rand   a pseudo-random number between 0. and 1.
c
c ia1 and ia0 are the hi and lo parts of a.  ia1ma0 = ia1 - ia0.
      data ia1, ia0, ia1ma0 /1536, 1029, 507/
      data ic /1731/
      data ix1, ix0 /0, 0/
c
      if (r.lt.0.) go to 10
      if (r.gt.0.) go to 20
c
c           a*x = 2**22*ia1*ix1 + 2**11*(ia1*ix1 + (ia1-ia0)*(ix0-ix1)
c                   + ia0*ix0) + ia0*ix0
c
      iy0 = ia0*ix0
      iy1 = ia1*ix1 + ia1ma0*(ix0-ix1) + iy0
      iy0 = iy0 + ic
      ix0 = mod (iy0, 2048)
      iy1 = iy1 + (iy0-ix0)/2048
      ix1 = mod (iy1, 2048)
c
 10   rand = ix1*2048 + ix0
      rand = rand / 4194304.
      return
c
 20   ix1 = amod(r,1.)*4194304. + 0.5
      ix0 = mod (ix1, 2048)
      ix1 = (ix1-ix0)/2048
      go to 10
c
      end
