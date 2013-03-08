!<compile=optimized>

!The 3D-RISM-KH software found here is copyright (c) 2010-2012 by 
!Andriy Kovalenko, Tyler Luchko, Takeshi Yamazaki and David A. Case.
!
!This program is free software: you can redistribute it and/or modify it
!under the terms of the GNU General Public License as published by the Free
!Software Foundation, either version 3 of the License, or (at your option)
!any later version.
!
!This program is distributed in the hope that it will be useful, but
!WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
!or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
!for more details.
!
!You should have received a copy of the GNU General Public License in the
!../../LICENSE file.  If not, see <http://www.gnu.org/licenses/>.
!
!Users of the 3D-RISM capability found here are requested to acknowledge
!use of the software in reports and publications.  Such acknowledgement
!should include the following citations:
!
!1) A. Kovalenko and F. Hirata. J. Chem. Phys., 110:10095-10112  (1999); 
!ibid. 112:10391-10417 (2000).   
!
!2) A. Kovalenko,  in:  Molecular  Theory  of  Solvation,  edited  by  
!F. Hirata  (Kluwer Academic Publishers, Dordrecht, 2003), pp.169-275.  
!
!3) T. Luchko, S. Gusarov, D.R. Roe, C. Simmerling, D.A. Case, J. Tuszynski,
!and  A. Kovalenko, J. Chem. Theory Comput., 6:607-624 (2010). 
#include "../include/dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Closure super class for 1D-RISM.  Closure sub-classes, (i.e., actual closure 
!!!implementations) are registered here.  Subroutine calls then call the 
!!!appropriate subroutine of the subclass. interface.  This is an explicit 
!!!implementation of class inheritance. See V. K. Decyk, C. D. Norton, 
!!!B. K. Szymanski.  How to express C++ concepts in Fortran 90. Scientific 
!!!Programming. 6, 363-390 (1997).
!!!
!!!Some closure independent properties are calculated within this class. Uvv is
!!!the site-site potential.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module rism1d_closure_c
  !add new closure modules here
  use rism1d_kh_c
  use rism1d_hnc_c
  use rism1d_py_c
  use rism1d_mv0_c
  use rism1d_psen_c
  use rism1d_potential_c
  use safemem
  implicit none

  integer,private ,parameter :: maxep0=4
  type rism1d_closure
     _REAL_, pointer :: xvv_dT(:,:,:) => NULL()
     !pressk : k-space pressure contribution
     !fek    : k-space free energy contribution
     _REAL_ :: pressk=HUGE(1d0), fek=HUGE(1d0)
     !xvv   : site-site succeseptibility.
     _REAL_, pointer :: xvv(:,:,:) => NULL()
     type(rism1d_potential),pointer :: pot => NULL()
     type(rism1d_kh),pointer :: kh => NULL()
     type(rism1d_hnc),pointer :: hnc => NULL()
     type(rism1d_py),pointer :: py => NULL()
     type(rism1d_mv0),pointer :: mv0 => NULL()
     type(rism1d_psen),pointer :: psen => NULL()
     character(len=8) :: type

  end type rism1d_closure

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Creates a new closure object of the requested type.
!!!IN:
!!!   this  : the closure object
!!!   type : one of 'KH', 'HNC', 'MV0', 'PSEn', 'V*', where 'n' is the
!!!          order of the  PSE-n closure
!!!   pot   : rism1d_potential object.  Must be initialized.
!!!   coeff : (optional) coefficients for the selected closure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_closure_new(this,type,pot,coeff)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    type(rism1d_potential), target, intent(in)  :: pot
    character(len=*), intent(in) :: type
    integer :: order, iostat
    _REAL_, optional, intent(in) :: coeff(:)

    !............... reset pressk and fek so they will be recalculated
    this%pressk = HUGE(1d0)
    this%fek = HUGE(1d0)

    this%pot => pot

    this%type=trim(type)
    if(type .eq. "KH") then
       allocate(this%kh)
       call rism1d_kh_new(this%kh)
    elseif(index(type,"PSE") ==1) then
       read(type(4:),*, iostat=iostat) order
       if(iostat/=0)&
           call rism_report_error(trim(type)//" not a valid closure")
       allocate(this%psen)
       call rism1d_psen_new(this%psen,order)
    elseif(type .eq. "HNC") then
       allocate(this%hnc)
       call rism1d_hnc_new(this%hnc)
    elseif(type .eq. "PY") then
       allocate(this%py)
       call rism1d_py_new(this%py)
    elseif(type .eq. "MV0") then
       allocate(this%mv0)
       call rism1d_mv0_new(this%mv0)
    else
       call rism_report_error(trim(type)//" not a valid closure")
    end if
  end subroutine rism1d_closure_new

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a identifier string for the closure type
!!!IN:
!!!   this : the closure object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_type(this) result(type)
    implicit none
    type(rism1d_closure), intent(in) :: this
    character(len=4) :: type
    type=this%type
  end function rism1d_closure_type

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns true if the full potential is being used.  This is always true except
!!!for MV0
!!!IN:
!!!   this : the closure object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_isCharged(this) result(charged)
    implicit none
    type(rism1d_closure), intent(in) :: this
    logical :: charged
    if(associated(this%MV0))then
       charged = rism1d_mv0_charged(this%mv0)
       return
    endif
    charged = .true.
  end function rism1d_closure_isCharged

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Set to true to use the full potential, false otherwise.  This only applies to
!!!MV0
!!!IN:
!!!   this : the closure object
!!!   charged : true to use full potential
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_closure_useCharged(this,charged)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    logical,intent(in) :: charged
    if(associated(this%MV0))then
       call rism1d_mv0_usecharged(this%mv0, charged)
    endif
  end subroutine rism1d_closure_useCharged

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Gvv from Uvv, Hvv, and Cvv using the associated closure
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site pair correlation function
!!!   hvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_closure_gvv(this,gvv, hvv, cvv)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(out) :: gvv(:,:)
    _REAL_, intent(in) :: hvv(:,:),cvv(:,:)
    if(associated(this%KH))then
       call rism1d_kh_gvv(this%kh,gvv,this%pot%uvv,hvv,cvv)
    elseif(associated(this%PSEN))then
       call rism1d_psen_gvv(this%psen,gvv,this%pot%uvv,hvv,cvv)!
    elseif(associated(this%HNC))then
       call rism1d_hnc_gvv(this%hnc,gvv,this%pot%uvv,hvv,cvv)
    elseif(associated(this%PY))then
       call rism1d_py_gvv(this%py,gvv,this%pot%uvv,hvv,cvv)
    elseif(associated(this%MV0))then
       call rism1d_mv0_gvv(this%mv0,gvv,this%pot%uvv,hvv,cvv)
    end if
  end subroutine rism1d_closure_gvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates temperature derivative Gvv from Gvv, Uvv, Hvv_dT, and Cvv_dT using
!!!the associated closure
!!!IN:
!!!   this   : the closure object
!!!   gvv_dT  : site-site temperature derivative pair correlation function
!!!   gvv    : site-site pair correlation function
!!!   cvv  : site-site direct correlation function
!!!   hvv_dT  : site-site temperature derivative total correlation function
!!!   cvv_dT  : site-site temperature derivative direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_closure_gvv_dT(this,gvv_dT, gvv, cvv, hvv_dT, cvv_dT)
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(out) :: gvv_dT(:,:)
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:), hvv_dT(:,:),cvv_dT(:,:)
    if(associated(this%KH))then
       call rism1d_kh_gvv_dT(this%kh,gvv_dT,this%pot%uvv,gvv,hvv_dT,cvv_dT)
    elseif(associated(this%PSEN))then
       call rism1d_psen_gvv_dT(this%psen,gvv_dT,this%pot%uvv,gvv,cvv,hvv_dT,cvv_dT)
    elseif(associated(this%HNC))then
       call rism1d_hnc_gvv_dT(this%hnc,gvv_dT,this%pot%uvv,gvv,hvv_dT,cvv_dT)
    else
       call rism_report_error("Temperature derivative not supported for "//this%type)
    end if
  end subroutine rism1d_closure_gvv_dT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates Bvv (bridge function) from Uvv, Gvv, and Cvv using the associated closure
!!!IN:
!!!   this : the closure object
!!!   gvv  : site-site total correlation function
!!!   cvv  : site-site direct correlation function
!!!OUT:
!!!   bvv  : site-site birdge function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_bvv(this, gvv, cvv) result(bvv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, pointer :: bvv(:,:)
    _REAL_, intent(in) :: gvv(:,:),cvv(:,:)
    nullify(bvv)
    if(associated(this%KH))then
       bvv=> rism1d_kh_bvv(this%kh,this%pot%uvv,gvv,cvv)
    elseif(associated(this%PSEN))then
       bvv=> rism1d_psen_bvv(this%psen,this%pot%uvv,gvv,cvv)
    elseif(associated(this%HNC))then
       bvv=> rism1d_hnc_bvv(this%hnc,this%pot%uvv,gvv,cvv)
    elseif(associated(this%PY))then
       bvv=> rism1d_py_bvv(this%py,this%pot%uvv,gvv,cvv)
    elseif(associated(this%MV0))then
       bvv=> rism1d_mv0_bvv(this%mv0,this%pot%uvv,gvv,cvv)
    end if
  end function rism1d_closure_bvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Frees memory and resets object state
!!!IN:
!!!   this : the closure object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_closure_destroy(this)
    use safemem
    implicit none
    type(rism1d_closure), intent(inout) :: this
    !............... reset pressk and fek so they will be recalculated
    this%pressk = HUGE(1d0)
    this%fek = HUGE(1d0)
    if(associated(this%KH))then
       call rism1d_kh_destroy(this%kh)
       deallocate(this%kh)
    end if
    if(associated(this%PSEN))then
       call rism1d_psen_destroy(this%psen)
       deallocate(this%psen)
    end if
    if(associated(this%HNC))then
       call rism1d_hnc_destroy(this%hnc)
       deallocate(this%hnc)
    end if
    if(associated(this%PY))then
       call rism1d_py_destroy(this%py)
       deallocate(this%py)
    end if
    if(associated(this%MV0))then
       call rism1d_mv0_destroy(this%mv0)
       deallocate(this%mv0)
    end if
    nullify(this%pot)
    if(safemem_dealloc(this%xvv) /=0)then
       call rism_report_error("deallocating Xvv failed in closure")
    end if
    if(safemem_dealloc(this%xvv_dT) /=0)then
       call rism_report_error("deallocating Xvv_dT failed in closure")
    end if
  end subroutine rism1d_closure_destroy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the compressibility of the solvent [A^3]
!!!IN:
!!!   this : rism1d closure object
!!!   cvv  : direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getCompressibility(this,cvv) result(xikt)
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: cvv(:,:)
    _REAL_ :: xikt
    _REAL_ :: r, ck0, ck0r
    integer :: ir, ivv, iv2, iv1, msym
    ck0 = 0.d0
    do ir=2,this%pot%nr
       r = (ir-1)*this%pot%dr
       ck0r = 0.d0
       ivv = 0
       do iv2=1,this%pot%nv
          do iv1=1,iv2
             ivv = ivv + 1
             if (iv1 == iv2)  then
                msym = 1
             else
                msym = 2
             endif
             ck0r = ck0r + msym*this%pot%rhov(iv1)*this%pot%rhov(iv2)*cvv(ir,ivv)
          enddo
       enddo
       ck0 = ck0 + r**2*ck0r
    enddo
    ck0 = ck0 * 4.d0*pi*this%pot%dr
    xikt = 1.d0 / (sum(this%pot%rhosp)-ck0)
  end function rism1d_closure_getCompressibility

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the the extrapolated value for
!!!DelHv=-Lim_k->0 ( Sum_v1 Qv1*Xv1v2(k)4pi/k^2 - hlkv0 )
!!!This is used by 3D-RISM long range asymptotics
!!!IN:
!!!   this : rism1d closure object
!!!   hvv  : total correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getDelHvLimit(this,hvv) result(delhv0)
    use constants, only : pi
    use rism_util, only : poly_interp
    implicit none
#include "../xblas/f77/blas_namedconstants.fh"    
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: hvv(:,:)
    _REAL_ :: delhv0(this%pot%nv)
     !hlkv0 : molecular version of hlkvv for 3D-RISM
    _REAL_ :: hlkv0(maxep0+1,this%pot%nv)
    integer :: ivv, iv1, iv2, ir
    _REAL_ :: k, ep0(maxep0), cep0(0:maxep0), err0

    !.....getting the long-range asymptotic function of H(k) for 3D-RISM
    do iv2=1,this%pot%nv
       do ir=2,maxep0+1
          k = (ir-1)*this%pot%dk
          hlkv0(ir,iv2) =-this%pot%qspv(iv2) /this%pot%dielconst &
               * 4.d0*pi*exp(-(0.5d0*this%pot%smear*k)**2) /(k**2+this%pot%kappa**2)
       enddo
    enddo
    call rism1d_closure_calcXvv(this,hvv)
    do ir=2,maxep0+1
       ep0(ir-1) = (ir-1)*this%pot%dk
    end do
    do iv2=1,this%pot%nv
       call DGEMV("N",maxep0+1,this%pot%nv,-4d0*pi,this%xvv(:maxep0+1,:,iv2),maxep0+1,&
            this%pot%qv,1,0d0,cep0,1)
!!$       call BLAS_DGEMV_X(BLAS_NO_TRANS,maxep0+1,this%pot%nv,-4d0*pi,this%xvv(:maxep0+1,:,iv2),maxep0+1,&
!!$            this%pot%qv,1,0d0,cep0,1,BLAS_PREC_EXTRA)
       cep0(1:maxep0) = cep0(1:maxep0)/ep0**2
       call DAXPY(maxep0,-1d0,hlkv0(2:maxep0+1,iv2),1,&
            cep0(1:maxep0),1)
!!$       call BLAS_DAXPBY_X(maxep0,-1d0,hlkv0(2:maxep0+1,iv2),1,&
!!$            1d0,cep0(1:maxep0),1,BLAS_PREC_EXTRA)
       call  poly_interp (ep0,cep0(1:maxep0),maxep0, 0.d0,delhv0(iv2), err0)
    enddo
  end function rism1d_closure_getDelHvLimit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns the the extrapolated value for the temperature derivative of 
!!!DelHv=-Lim_k->0 ( Sum_v1 Qv1*Xv1v2_dT(k)4pi/k^2 - hlkv0 )
!!!This is used by 3D-RISM long range asymptotics
!!!IN:
!!!   this   : rism1d closure object
!!!   hvv    : total correlation function
!!!   hvv_dT : temperature derivative total correlation function
!!!OUT:
!!!   temperature derivative delhv0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getDelHvLimit_DT(this,hvv,hvv_dT) result(delhv0_dT)
    use constants, only : pi
    use rism_util, only : poly_interp
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: hvv(:,:), hvv_dT(:,:)
    _REAL_ :: delhv0_dT(this%pot%nv)
     !hlkv0 : molecular version of hlkvv for 3D-RISM
    _REAL_ :: hlkv0(maxep0+1,this%pot%nv)
    integer :: ivv, iv1, iv2, ir
    _REAL_ :: k, ep0(maxep0), cep0(maxep0), err0

    !.....getting the long-range asymptotic function of H(k) for 3D-RISM
!!    do iv2=1,this%pot%nv
!!       do ir=2,maxep0+1
!!          k = (ir-1)*this%pot%dk
!!          hlkv0(ir,iv2) =-this%pot%qspv(iv2) /this%pot%dielconst &
!!               * 4.d0*pi*exp(-(0.5d0*this%pot%smear*k)**2) /(k**2+this%pot%kappa**2)
!!       enddo
!!    enddo
!!!    call rism1d_closure_calcXvv(this,hvv)
    call rism1d_closure_calcXvv_DT(this,hvv_dT)
    do iv2=1,this%pot%nv
       do ir=2,maxep0+1
          k = (ir-1)*this%pot%dk
          ep0(ir-1) = k
          cep0(ir-1) = 0.d0
          do iv1=1,this%pot%nv
!!             cep0(ir-1) = cep0(ir-1) + this%pot%qv(iv1)*(this%xvv_dT(ir,iv1,iv2) - this%xvv(ir,iv1,iv2))
             cep0(ir-1) = cep0(ir-1) + this%pot%qv(iv1)*this%xvv_dT(ir,iv1,iv2)
          enddo

          cep0(ir-1) = - 4.d0*pi*cep0(ir-1)/k**2
!!          cep0(ir-1) = cep0(ir-1) + hlkv0(ir,iv2)
          cep0(ir-1) = cep0(ir-1)
       enddo
       call  poly_interp (ep0,cep0,4, 0.d0,delhv0_dT(iv2), err0)
    enddo
    delhv0_dT = delhv0_dT - rism1d_closure_getDelHvLimit(this,hvv)
  end function rism1d_closure_getDelHvLimit_DT

!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!If Xvv has not been calculated since our last solution, allocate memory and
!!!calculate it
!!!IN:
!!!   this : rism1d closure object
!!!   hvv  : total correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_closure_calcXvv(this,hvv)
    implicit none
#include "../xblas/f77/blas_namedconstants.fh"    
    type(rism1d_closure),intent(inout) :: this
    _REAL_, intent(in) :: hvv(:,:)
    integer :: iv1, iv2, ivv, ir
    if(associated(this%xvv)) return
    this%xvv => safemem_realloc(this%xvv,this%pot%nr,this%pot%nv,this%pot%nv)
    this%xvv=0
    !.................. getting Xvv(k)=Wvv(k)+RhoV*Hvv(k) ..................
    ivv = 0
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          ivv = ivv + 1
!!$          call BLAS_DWAXPBY_X(this%pot%nr,1d0,this%pot%wvv(:,iv1,iv2),1,&
!!$               this%pot%rhov(iv1),hvv(:,ivv),1,&
!!$               this%xvv(:,iv1,iv2),1,BLAS_PREC_EXTRA)
!!$          if (iv1 /= iv2) &
!!$               call BLAS_DWAXPBY_X(this%pot%nr,1d0,this%pot%wvv(:,iv2,iv1),1,&
!!$               this%pot%rhov(iv2),hvv(:,ivv),1,&
!!$               this%xvv(:,iv2,iv1),1,BLAS_PREC_EXTRA)
          do ir=1,this%pot%nr
             this%xvv(ir,iv1,iv2) = this%pot%wvv(ir,iv1,iv2) + this%pot%rhov(iv1)*hvv(ir,ivv)
             if (iv1 /= iv2) &
                  this%xvv(ir,iv2,iv1) = this%pot%wvv(ir,iv2,iv1) + this%pot%rhov(iv2)*hvv(ir,ivv)
          enddo
       enddo
    enddo
  end subroutine rism1d_closure_calcXvv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!If Xvv_dT has not been calculated since our last solution, allocate memory and
!!!calculate it
!!!IN:
!!!   this    : rism1d closure object
!!!   hvv_dT  : temperature derivative total correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine rism1d_closure_calcXvv_DT(this,hvv_dT)
    implicit none
    type(rism1d_closure),intent(inout) :: this
    _REAL_, intent(in) :: hvv_dT(:,:)
    integer :: iv1, iv2, ivv, ir
    if(associated(this%xvv_dT)) return
    this%xvv_dT => safemem_realloc(this%xvv_dT,this%pot%nr,this%pot%nv,this%pot%nv)
    this%xvv_dT=0
    !.................. getting Xvv(k)=Wvv(k)+RhoV*Hvv(k) ..................
    ivv = 0
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          ivv = ivv + 1
          do ir=1,this%pot%nr
             this%xvv_dT(ir,iv1,iv2) = this%pot%rhov(iv1)*hvv_dT(ir,ivv)
             if (iv1 /= iv2) &
                  this%xvv_dT(ir,iv2,iv1) = this%pot%rhov(iv2)*hvv_dT(ir,ivv)
          enddo
       enddo
    enddo
  end subroutine rism1d_closure_calcXvv_DT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns an NV X NV array of total excess coordination numbers excluding
!!!multiplicity
!!!IN:
!!!   this : rism1d closure object
!!!   hvv  : total correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getExNumber(this,hvv) result(exvv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: hvv(:,:)
    _REAL_ :: exvv(this%pot%nv,this%pot%nv)
    integer :: ivv, iv1, iv2
    !......... excess number (excluding Mult) of sites 2 around 1 ..........
    ivv = 0
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          ivv = ivv + 1
          exvv(iv1,iv2) = this%pot%rhov(iv2)/this%pot%mtv(iv2) * hvv(1,ivv)
          exvv(iv2,iv1) = this%pot%rhov(iv1)/this%pot%mtv(iv1) * hvv(1,ivv)
       enddo
    enddo
  end function rism1d_closure_getExNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a pointer to an NR X NVV array of structure factors.  This memory 
!!!must be freed (preferably with safemem_dealloc) as it is not freed locally or
!!!after the object instance is destroyed.
!!!IN:
!!!   this : rism1d closure object
!!!   hvv  : total correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getStructFactor(this,hvv) result(svv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: hvv(:,:)
    _REAL_, pointer :: svv(:,:)
    _REAL_ :: rhotot
    integer :: ivv, iv1, iv2, ir
    nullify(svv)
    svv=>safemem_realloc(svv, this%pot%nr, this%pot%nvv)
    rhotot = sum(this%pot%rhosp)
    ivv = 0
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          ivv = ivv + 1
          do ir=1,this%pot%nr
             svv(ir,ivv) = hvv(ir,ivv) &
                  * this%pot%rhov(iv1)/this%pot%mtv(iv1) * this%pot%rhov(iv2)/this%pot%mtv(iv2) &
                  / rhotot
          enddo
       enddo
    enddo
  end function rism1d_closure_getStructFactor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a pointer to an NR X NVV array of the running site-site excess number.
!!!This is the excess number of a site within a given radius. The memory for 
!!!this pointer must be freed (preferably with safemem_dealloc) as it is not 
!!!freed locally or after the object instance is destroyed.
!!!IN:
!!!   this : rism1d closure object
!!!   gvv  : pair distribution function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getRunExNumber(this,gvv) result(exnvv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:)
    _REAL_, pointer :: exnvv(:,:,:)
    nullify(exnvv)
    exnvv => exNumber(this,gvv,.true.)
  end function rism1d_closure_getRunExNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a pointer to an NR X NVV array of the running site-site number.
!!!This is the number of a site within a given radius. The memory for 
!!!this pointer must be freed (preferably with safemem_dealloc) as it is not 
!!!freed locally or after the object instance is destroyed.
!!!IN:
!!!   this : rism1d closure object
!!!   gvv  : pair distribution function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getRunNumber(this,gvv) result(nvv)
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:)
    _REAL_, pointer :: nvv(:,:,:)
    nullify(nvv)
    nvv => exNumber(this,gvv,.false.)
  end function rism1d_closure_getRunNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the pressure in [kT / A^3] of the system using the virial 
!!!path.  To convert to Pacals, for example, multiply by 
!!!1.d30 * kb * temperature
!!!where kb is Boltzmann's constant [j/K] and temperature is in [K]
!!!IN:
!!!   this : rism1d closure object
!!!   gvv  : pair distribution function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getPressureVirial(this,gvv) result(pressure)
    use constants, only :pi
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: gvv(:,:)
    _REAL_ :: pressure
    _REAL_ :: r
    integer :: ir, ivv, iv1, iv2, cnt

    !make sure the derivative of the potential is calculated
    call rism1d_potential_derivative(this%pot)

    !calculate the integral for each site
    pressure=0
    ivv=0
    do iv1 = 1, this%pot%nv
       do iv2 = 1, iv1 
          ivv = ivv + 1
          do ir = 2, this%pot%nr
             r = (ir -1)*this%pot%dr
             !the full eqn. is a double sum over sites. so, if the
             !sites are different, we need to double count the
             !contribution.
             if (iv1 == iv2)  then
                cnt = 1
             else
                cnt = 2
             endif
             pressure = pressure + &
                  cnt*this%pot%rhov(iv1)*this%pot%rhov(iv2)*&
                  gvv(ir, ivv)*this%pot%duvv(ir,ivv)*r**3
          end do
       end do
    end do
    pressure = sum(this%pot%rhov) - this%pot%dr*PI*2d0/3d0*pressure
  end function rism1d_closure_getPressureVirial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the pressure in [kT / A^3] of the system using the free energy 
!!!path.  To convert to Pacals, for example, multiply by 
!!!1.d30 * kb * temperature
!!!where kb is Boltzmann's constant [j/K] and temperature is in [K]
!!!IN:
!!!   this : rism1d closure object
!!!   gvv  : pair distribution function
!!!   cvv  : direct correlation function
!!!   rhotrgt : (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. HNC and KH)
!!!OUT:
!!!    pressure from the free energy route if defined for the given closure.  If not, returns HUGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getPressureFE(this,gvv,cvv,rhotrgt) result(pressure)
    use safemem
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    _REAL_, optional :: rhotrgt 
    _REAL_ :: pressure
    !pressr : r-space pressure contribution
    !rtrgt  : target density
    _REAL_ :: pressr, rtrgt

    if(present(rhotrgt))then
       rtrgt = rhotrgt
    else
       rtrgt = sum(this%pot%rhosp)
    end if

    !r-space contribution. This is closure dependent
    if(associated(this%KH))then
       pressr = rism1d_kh_pressureFE(this%kh,gvv,cvv,this%pot%mtv,this%pot%rhov,this%pot%dr)
    elseif(associated(this%PSEN))then
       pressr = rism1d_psen_pressureFE(this%psen,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%rhov,this%pot%dr)
    elseif(associated(this%HNC))then
       pressr = rism1d_hnc_pressureFE(this%hnc,gvv,cvv,this%pot%mtv,this%pot%rhov,this%pot%dr)
    else
       pressure = HUGE(1d0)
       return
    end if

    !k-space contribution. This is not closure dependent
    call fe_press_k(this,cvv)

    !check that the k-space contribution was correctly calculated
    if(this%pressk == tiny(0d0))then
       pressure = 0
       return
    end if

    pressure = rtrgt*(1.d0 + (rtrgt/sum(this%pot%rhosp))*pressr + this%pressk)
  end function rism1d_closure_getPressureFE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the freeEnergy in kT
!!!IN:
!!!   this : rism1d object
!!!   gvv  : pair distribution function
!!!   cvv  : direct correlation function
!!!   rhotrgt : (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. HNC and KH)
!!!OUT:
!!!    free energy if defined for the given closure.  If not, returns HUGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getFreeEnergy(this,gvv,cvv,rhotrgt) result(fe)
    use safemem
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    _REAL_, optional :: rhotrgt 
    _REAL_ :: fe
    !fer : r-space fe contribution
    !rtrgt  : target density
    _REAL_ :: fer, rtrgt

    if(present(rhotrgt))then
       rtrgt = rhotrgt
    else
       rtrgt = sum(this%pot%rhosp)
    end if

    !r-space contribution. This is closure dependent
    if(associated(this%KH))then
       fer = rism1d_kh_freeEnergy(this%kh,gvv,cvv,this%pot%mtv,this%pot%rhov,this%pot%dr)
    elseif(associated(this%PSEN))then
       fer = rism1d_psen_freeEnergy(this%psen,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%rhov,this%pot%dr)
    elseif(associated(this%HNC))then
       fer = rism1d_hnc_freeEnergy(this%hnc,gvv,cvv,this%pot%mtv,this%pot%rhov,this%pot%dr)
    else
       fe = HUGE(1d0)
       return
    end if

    !k-space contribution. This is not closure dependent
    call fe_press_k(this,cvv)

    !check that the k-space contribution was correctly calculated
    if(this%fek == tiny(0d0))then
       fe = 0
       return
    end if
    fe = (this%fek + fer)/sum(this%pot%rhosp)
  end function rism1d_closure_getFreeEnergy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the partial molar volume of each species in A^3
!!!IN:
!!!   this : rism1d closure object
!!!   cvv  : direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getPMV(this,cvv) result(pmv)
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: cvv(:,:)
    _REAL_ :: pmv(this%pot%nsp)
    _REAL_ :: compressibility, r, cvvk0, cvvk0r
    integer :: ivv, iv1, iv2, iv21, iv22, ir, isp, iat
    compressibility = rism1d_closure_getCompressibility(this,cvv)
    iv22 = 0
    do isp=1,this%pot%nsp
       iv21 = iv22 + 1
       do iat=1,this%pot%nat(isp)
          iv22 = iv22 + 1
       enddo
       !calculate Cvv(k=0)
       cvvk0 = 0.d0
       do ir=2,this%pot%nr
          r = (ir-1)*this%pot%dr
          cvvk0r = 0.d0
          do iv2=iv21,iv22
             do iv1=1,this%pot%nv
                ivv = this%pot%jvv(iv1,iv2)
                cvvk0r = cvvk0r + this%pot%rhov(iv1)*this%pot%mtv(iv2)*cvv(ir,ivv)
             enddo
          enddo
          cvvk0 = cvvk0 + r**2*cvvk0r
       enddo
       cvvk0 = cvvk0 * 4.d0*pi*this%pot%dr
       pmv(isp) = compressibility*(1.d0-cvvk0)
    enddo
  end function rism1d_closure_getPMV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the excess chemical potential in kT for each site.
!!!IN:
!!!   this : rism1d closure object
!!!   gvv  : pair distribution function
!!!   cvv  : direct correlation function
!!!   rhotrgt : (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. HNC and KH)
!!!OUT:
!!!    excess chemical potential if defined for the given closure.  If not, returns HUGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getExChem(this,gvv,cvv,rhotrgt) result(exchem)
    use safemem
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    _REAL_, optional :: rhotrgt 
    _REAL_ :: exchem(this%pot%nv)
    !rhotrgtv  : target density
    _REAL_ :: rhotrgtv(this%pot%nv)
    integer :: isp, iat, iv

    if(present(rhotrgt))then
       iv = 0
       do isp=1,this%pot%nsp
          do iat=1,this%pot%nat(isp)
             iv = iv + 1
             rhotrgtv(iv) = this%pot%rhov(iv) *rhotrgt/sum(this%pot%rhosp)
          enddo
       enddo
    else
       rhotrgtv = this%pot%rhov
    end if
    if(associated(this%KH))then
       exchem = rism1d_kh_exChem(this%kh,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%rhov,rhotrgtv,this%pot%dr)
    elseif(associated(this%PSEN))then
       exchem = rism1d_psen_exChem(this%psen,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%rhov,rhotrgtv,this%pot%dr)
    elseif(associated(this%HNC))then
       exchem = rism1d_hnc_exChem(this%hnc,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%rhov,rhotrgtv,this%pot%dr)
    else
       exchem = HUGE(1d0)
    end if
  end function rism1d_closure_getExChem

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the ionic(?) excess chemical potential in kT for each site.
!!!IN:
!!!   this : rism1d closure object
!!!   gvv  : pair distribution function
!!!   cvv  : direct correlation function
!!!   rhotrgt : (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. HNC and KH)
!!!OUT:
!!!    ionic excess chemical potential if defined for the given closure.  If not, returns HUGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getExChemIon(this,gvv,cvv,rhotrgt) result(exchem)
    use safemem
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:), cvv(:,:)
    _REAL_, optional :: rhotrgt 
    _REAL_ :: exchem(this%pot%nv)
    !rhotrgtv  : target density
    _REAL_ :: rhotrgtv(this%pot%nv)
    integer :: isp, iat, iv

    if(present(rhotrgt))then
       iv = 0
       do isp=1,this%pot%nsp
          do iat=1,this%pot%nat(isp)
             iv = iv + 1
             rhotrgtv(iv) = this%pot%rhov(iv) *rhotrgt/sum(this%pot%rhosp)
          enddo
       enddo
    else
       rhotrgtv = this%pot%rhov
    end if
    if(associated(this%KH))then
       exchem = rism1d_kh_exChemIon(this%kh,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%rhov,rhotrgtv,this%pot%dr)
    elseif(associated(this%PSEN))then
       exchem = rism1d_psen_exChemIon(this%psen,gvv,this%pot%uvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%rhov,rhotrgtv,this%pot%dr)
    elseif(associated(this%HNC))then
       exchem = rism1d_hnc_exChemIon(this%hnc,gvv,cvv,this%pot%mtv,this%pot%jvv,&
            this%pot%rhov,rhotrgtv,this%pot%dr)
    else
       exchem = HUGE(1d0)
    end if
  end function rism1d_closure_getExChemIon

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculates the solvation energy in kT for each site.
!!!IN:
!!!   this   : rism1d closure object
!!!   gvv    : pair distribution function
!!!   uvv    : site-site potential
!!!   cvv    : direct correlation function
!!!   gvv_dT : temperature derivative pair distribution function
!!!   cvv_dT : temperature derivative direct correlation function
!!!   rhotrgt: (optional) The final, physical density for thermodynamics.  This
!!!             can be used as an effective correction for some closures that
!!!             over estimate the pressure (e.g. HNC and KH)
!!!OUT:
!!!    solvation energy if defined for the given closure.  If not, returns HUGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function rism1d_closure_getSolvene(this,gvv,uvv,cvv,gvv_dT,cvv_dT,rhotrgt) result(solvene)
    use safemem
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:), uvv(:,:), cvv(:,:), gvv_dT(:,:), cvv_dT(:,:)
    _REAL_, optional :: rhotrgt 
    _REAL_ :: solvene(this%pot%nv)
    !rhotrgtv  : target density
    _REAL_ :: rhotrgtv(this%pot%nv)
    integer :: isp, iat, iv

    if(present(rhotrgt))then
       iv = 0
       do isp=1,this%pot%nsp
          do iat=1,this%pot%nat(isp)
             iv = iv + 1
             rhotrgtv(iv) = this%pot%rhov(iv) *rhotrgt/sum(this%pot%rhosp)
          enddo
       enddo
    else
       rhotrgtv = this%pot%rhov
    end if
    if(associated(this%KH))then
       solvene = rism1d_kh_Solvene(this%kh,gvv,cvv,gvv_dT,cvv_dT,this%pot%mtv,this%pot%jvv,&
            this%pot%rhov,rhotrgtv,this%pot%dr)
    elseif(associated(this%PSEN))then
       solvene = rism1d_psen_Solvene(this%psen,gvv,uvv,cvv,gvv_dT,cvv_dT,this%pot%mtv,this%pot%jvv,&
            this%pot%rhov,rhotrgtv,this%pot%dr)
    elseif(associated(this%HNC))then
       solvene = rism1d_hnc_Solvene(this%hnc,gvv,cvv,gvv_dT,cvv_dT,this%pot%mtv,this%pot%jvv,&
            this%pot%rhov,rhotrgtv,this%pot%dr)
    else
       solvene = HUGE(1d0)
   end if
  end function rism1d_closure_getSolvene

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!PRIVATE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Returns a pointer to an NR X NV X NV array of the running (excess) site-site number.
!!!This is the (excess) number of a site within a given radius. The memory for 
!!!this pointer must be freed (preferably with safemem_dealloc) as it is not 
!!!freed locally or after the object instance is destroyed.
!!!IN:
!!!   this : rism1d object
!!!   gvv  : pair distribution function
!!!   excess : .true. for excess, .false. for total
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function exNumber(this,gvv,excess) result(exnvv)
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(in) :: this
    _REAL_, intent(in) :: gvv(:,:)
    logical, intent(in) ::excess
    _REAL_, pointer :: exnvv(:,:,:)
    _REAL_ :: tempvv0, tempvv1
    integer :: ivv, iv1, iv2, ir
    nullify(exnvv)
    exnvv => safemem_realloc(exnvv, this%pot%nr, this%pot%nv, this%pot%nv)

    !the cummulative integral is done inplace using the trapazoidal rule
    !for either case of hvv or gvv.  This means using temporary variables to
    !hold values of hvv and gvv before they are overwritten
    ivv = 0
    do iv2=1,this%pot%nv
       do iv1=1,iv2
          ivv = ivv + 1
          !select hvv or gvv
          if(excess)then
             exnvv(:,iv1,iv2) = gvv(:,ivv)-1d0
          else
             exnvv(:,iv1,iv2) = gvv(:,ivv)
          end if
          !integrate
          tempvv0 = exnvv(1,iv1,iv2)
          exnvv(1,iv1,iv2) = 0.d0
          do ir=2,this%pot%nr
             tempvv1 = exnvv(ir,iv1,iv2)
             exnvv(ir,iv1,iv2) = exnvv(ir-1,iv1,iv2) &
                  + 2.d0*pi*this%pot%dr &
                  *(((ir-2)*this%pot%dr)**2*tempvv0 &
                  + ((ir-1)*this%pot%dr)**2*exnvv(ir,iv1,iv2) )
             tempvv0=tempvv1
          enddo
       end do
    end do
    do iv2=1,this%pot%nv
       do iv1=1,iv2
!          exnvv(:,iv2,iv1) = this%pot%rhov(iv1)/this%pot%mtv(iv1) * exnvv(:,iv1,iv2)
          exnvv(:,iv2,iv1) = this%pot%rhov(iv2) * exnvv(:,iv1,iv2)
          if(iv1 /= iv2)&
!               exnvv(:,iv1,iv2) = this%pot%rhov(iv2)/this%pot%mtv(iv2) * exnvv(:,iv1,iv2)
               exnvv(:,iv1,iv2) = this%pot%rhov(iv1) * exnvv(:,iv1,iv2)
       end do
    end do
  end function exNumber

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!Calculate the k-space contributions to the free energy and pressure along the
!!!free energy path
!!!IN:
!!!   this : rism1d closure object
!!!   cvv  : direct correlation function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fe_press_k(this,cvv)
    use constants, only : pi
    implicit none
    type(rism1d_closure), intent(inout) :: this
    _REAL_, intent(in) :: cvv(:,:)
    _REAL_, pointer :: cvvk(:,:)=>NULL()
    _REAL_ :: bk(this%pot%nv,this%pot%nv), wck(this%pot%nv,this%pot%nv), ak(this%pot%nv,this%pot%nv),&
         det,detln
    _REAL_ :: k, r, presskv, fekv
    integer :: indx(this%pot%nv)
    integer :: ir, ivv, iv1, iv2, iv3, err
    if(this%fek /=HUGE(1d0))then
       return
    end if
    cvvk => safemem_realloc(cvvk, this%pot%nr, this%pot%nvv)
    !.............................. FFT>K Cvv ..............................
    do ivv=1,this%pot%nvv
       do ir=1,this%pot%nr
          cvvk(ir,ivv) = (ir-1)*this%pot%dr*cvv(ir,ivv)
       enddo
    enddo

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          cvvk(ir,ivv) = cvvk(ir,ivv) + this%pot%ulrvv(ir,ivv)
       enddo
    enddo

    do ivv=1,this%pot%nvv
       call  sinfti (cvvk(2,ivv),this%pot%nr-1, this%pot%dr, +1)
    enddo

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          cvvk(ir,ivv) = cvvk(ir,ivv) - this%pot%ulkvv(ir,ivv)
       enddo
    enddo

    do ivv=1,this%pot%nvv
       do ir=2,this%pot%nr
          cvvk(ir,ivv) = cvvk(ir,ivv) / ((ir-1)*this%pot%dk)
       enddo
    enddo

    this%fek = 0.d0
    this%pressk = 0.d0
    do ir=2,this%pot%nr
       k = (ir-1)*this%pot%dk

       !........................... loading Cvv(k) ............................
       ivv = 0
       do iv2=1,this%pot%nv
          do iv1=1,iv2
             ivv = ivv + 1
             bk(iv1,iv2) = cvvk(ir,ivv)
             bk(iv2,iv1) = bk(iv1,iv2)
          enddo
       enddo

       !..................... getting (Wvv+Zvv)^(tr)*Cvv ......................
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             wck(iv1,iv2) = 0.d0
             do iv3=1,this%pot%nv
                wck(iv1,iv2) = wck(iv1,iv2) + this%pot%wzvv(ir,iv3,iv1)*bk(iv3,iv2)
             enddo
          enddo
       enddo

       !................. getting (Wvv+Zvv)^(tr)*Cvv*Rho_v ..................
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             wck(iv1,iv2) = wck(iv1,iv2)*this%pot%rhov(iv2)
          enddo
       enddo

       !........... loading Tr [(Wvv+Zvv)^(tr)*k*Cvv*Rho_v] into xA ...........
       fekv = 0.d0
       do iv1=1,this%pot%nv
          fekv = fekv + wck(iv1,iv1)
       enddo

       !............... getting Ak=(1-(Wvv+Zvv)^(tr)*Cvv*Rho_v) ...............
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             ak(iv1,iv2) = - wck(iv1,iv2)
          enddo
       enddo

       do iv1=1,this%pot%nv
          ak(iv1,iv1) = 1.d0 + ak(iv1,iv1)
       enddo

       !............................ inverting Ak .............................
       !LAPACK LU-factorization
       !     call DGETRF(M,N,A,LDA,IPIV,INFO)
       call DGETRF(this%pot%nv,this%pot%nv,ak,this%pot%nv,indx,err)
       if(err > 0)then
          call rism_report_error("LU-factorization failed.  U = 0")
       elseif(err<0)then
          err = err*(-1)
          call rism_report_error("LU-factorization failed.")
       endif
       !DGETRF does not calculate the determinant of the permumation
       !maxtrix (as the Numerical Recipes subroutine does) so we
       !calculate it here IPIV is 'compressed' permutation matrix
       !where the index of the exchange row is record for each row
       !index.  So we count the number of rows that are not
       !exchanged with themselves and divide by two to get the total
       !number of exchanges.  The detminant is then
       !(-1)^number_of_exchanges
       det =1
       do iv1=1,this%pot%nv
          if(indx(iv1) /= iv1) det = det*(-1)
       end do

       !...................... checking for Det[Ak] > 0 ......................

       !we only care about the sign of the determinant so we use the
       !SIGN intrinsic and avoid possible overflows
       do iv1=1,this%pot%nv
          det = det * sign(1d0,ak(iv1,iv1))
       enddo
       if (det <= 0.d0)  then
          call rism_report_warn("fe_press_k: non-positive Det[1-Wvv*Cvv*Rho]")
          this%pressk = tiny(0.d0)
          this%fek = tiny(0.d0)
          return
       endif

       !................. calculating Ln(Det[1-Wvv*Cvv*Rho]) ..................
       detln = 0.d0
       do iv1=1,this%pot%nv
          detln = detln + log( abs( ak(iv1,iv1)) )
       enddo

       !....... adding Ln(Det[1-Wvv*Cvv*Rho]) to xA and loading into xP .......
       fekv = fekv + detln
       presskv = detln

       !............ calculating (1-Wvv*Cvv*Rho)^(-1)*Wvv*Cvv*Rho .............
       !        call DGETRS(TRANS, N, NRHS,A,LDA,IPIV,B,LDB,INFO)
       call DGETRS('N', this%pot%nv, this%pot%nv,ak,this%pot%nv,indx,wck,this%pot%nv,err)
       if(err < 0)then
          err = err*(-1)
          call rism_report_error("fe_press_k: Linear equation solver failed.")
       endif

       !............................ adding kZvv/k ............................
       do iv2=1,this%pot%nv
          do iv1=1,this%pot%nv
             wck(iv1,iv2) = wck(iv1,iv2) + this%pot%zkvv(ir,iv1,iv2)/k
          enddo
       enddo

       !............ adding Tr [(Wvv+Zvv)^(tr)*k*Cvv*Rho_v] to xP .............
       do iv1=1,this%pot%nv
          presskv = presskv + wck(iv1,iv1)
       enddo

       !.................. maintaining K-parts of xA and xP ...................
       this%fek = this%fek + fekv*k**2
       this%pressk = this%pressk + presskv*k**2
    enddo
    this%fek = this%fek/sum(this%pot%rhosp) * this%pot%dk/(2.d0*pi)**2
    this%pressk = -this%pressk/sum(this%pot%rhosp) * this%pot%dk/(2.d0*pi)**2
    if(safemem_dealloc(cvvk) /= 0)then
       call rism_report_error("dealloc failed in fe_press_k")
    end if
  end subroutine fe_press_k

end module rism1d_closure_c
