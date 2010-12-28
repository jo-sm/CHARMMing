      program trkswp
C temperature replica exchange; venabler AT nhlbi*nih*gov
C process log file to produce a swap track for each bath

      implicit none
      integer nbath, mxba, mxsw, istrt,istop, nclen, i,j,k
      parameter (mxba=128,mxsw=50000)
      integer ibath(mxba), itrk(mxsw,mxba), ndxbath
      integer kb1,kb2, j1,j2, isw,lsw, ksv
      real  pref
      character*80  str,tmp
      external ndxbath

C named common for large array
      common /trak/ itrk

C read the Trex config file
      j = iargc()
      if (j.ne.1) then
        write(6,*) 'Missing arg'
        write(6,*) 'trackswap.ix CONFIG-FILE'
        stop
      else
        call getarg(1,str)
      endif
      open(1,file=str,form='formatted',status='old')
      read(1,'(A5,I10)') tmp,istrt
      read(1,'(A5,I10)') tmp,istop
      read(1,'(A5,I10)') tmp,nclen
      read(1,'(A5,F10.2)') tmp,pref
      read(1,'(A5,I10)') tmp,nbath
      if (nbath.gt.mxba) stop 'MXBA limit exceeded'
      do i=1,nbath
        read(1,'(I10)',end=750,err=750) ibath(i)
      enddo
      close(1)
      write(6,701) istrt,istop,nclen,pref,nbath

 701  format('START=',I6,3X,'STOP=',I6,3X,'CLEAN=',I3,3X,
     *  'PREF=',F7.2,3X, 'NBATH=',I3)

C init array
      do j=1,nbath
        itrk(1,j)=ibath(j)
        do i=2,mxsw
          itrk(i,j)=0
        enddo
      enddo

      write(6,'(10I6)') (ibath(j),j=1,nbath)

      lsw=1
C read from stdin until EOF
  10  read(5,'(A80)',end=100,err=760) tmp
C check first 8 chars; get index of : and read 3 integers
      if (tmp(1:8).eq.'the foll') then
        k=2+index(tmp,':')
        read(tmp(k:k+19),'(2I5,I9)') kb1,kb2,isw
C update bath identity for a new swap step
        if (isw.gt.mxsw) goto 770
        if (isw.gt.lsw) then
          do j=1,nbath
            do i=lsw+1,isw
              itrk(i,j)=itrk(i-1,j)
            enddo
          enddo
          lsw = isw
        endif
C get index from current bath T
        j1 = ndxbath(kb1,nbath,ibath)
        j2 = ndxbath(kb2,nbath,ibath)
        if ((j1.eq.0).or.(j2.eq.0)) goto 780
C swap initial indenty
        ksv = itrk(isw,j1)
        itrk(isw,j1) = itrk(isw,j2)
        itrk(isw,j2) = ksv
      endif
      goto 10
 100  continue

C output as individual data files for each T
      do j=1,nbath
        kb1 = ibath(j)
        if (kb1.gt.999) then
          write (str,804) kb1 
        else
          write (str,803) kb1 
        endif
        open(2,file=str,form='formatted')
        do i=1,lsw
          j1 = 0
C output identity of bath with config of initial T, config in current bath
          do k=1,nbath
            if (kb1.eq.itrk(i,k)) j1 = k
          enddo
          write(2,810) real(i),real(ibath(j1)),real(itrk(i,j))
        enddo
      enddo
 803  format(I3,'/trkswp.dat')
 804  format(I4,'/trkswp.dat')
 810  format(F10.0,2F8.0)

      stop
 750  write(6,*) 'Error reading T bath values'
      stop
 760  write(6,*) 'Error reading log data'
      stop
 770  write(6,*) 'MAX SWAP STEP (mxsw) EXCEEDED'
      stop
 780  write(6,*) 'Bath index problem, ',kb1,kb2,j1,j2
      stop
      end


      integer function ndxbath(k,n,ib)
C index of a bath value
C k  in; bath value
C n  in; max bath index
C ib in; array of length n, bath values
      integer n,i,k,ib(n)
      ndxbath=0
      do i=1,n
        if (ib(i).eq.k) then
          ndxbath=i
          return
        endif
      enddo
      return
      end

