! **********************************************************************
! **                Coordinate transformation: main unit              **
! **                                                                  **
! **   Subunits:                                                      **
! **     dlf_neb.f90                                                  **
! **     dlf_dimer.f90                                                **
! **     dlf_hdlc_interface.f90                                       **
! **                                                                  **
! **                                                                  **
! **                                                                  **
! **********************************************************************
!!****h* DL-FIND/coords
!!
!! NAME
!! coords
!!
!! FUNCTION
!! Coordinate transformation: main unit
!!
!! Weight transformation is calculated at coordinate initialisation
!!
!! DATA
!! $Date: 2010-05-07 19:30:43 $                   
!! $Rev: 307 $                                                              
!! $Author: gberan $                                                         
!! $URL: http://ccpforge.cse.rl.ac.uk/svn/dl-find/branches/release_chemsh3.3/dlf_coords.f90 $   
!! $Id: dlf_coords.f90,v 1.1 2010-05-07 19:30:43 gberan Exp $                      
!!
!! COPYRIGHT
!!
!!  Copyright 2007 Johannes Kaestner (j.kaestner@dl.ac.uk),
!!  Tom Keal (keal@mpi-muelheim.mpg.de)
!!
!!  This file is part of DL-FIND.
!!
!!  DL-FIND is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as 
!!  published by the Free Software Foundation, either version 3 of the 
!!  License, or (at your option) any later version.
!!
!!  DL-FIND is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public 
!!  License along with DL-FIND.  If not, see 
!!  <http://www.gnu.org/licenses/>.
!!
!! SOURCE
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_coords_init
!!
!! FUNCTION
!! Initialise coordinate transformation
!!
!! In case of a direct coordinate transformation (glob%icoord<10):
!! * Find the number of degrees of freedom
!! * (Initialise HDLCs)
!! * Allocate global arrays (glob%icoords, glob%igradient, glob%step)
!!
!! In multiple image methods (NEB, dimer) call their respective init routines
!!
!! SYNOPSIS
subroutine dlf_coords_init
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,stdout,printl
  use dlf_allocate, only: allocate, deallocate
  implicit none
  integer               :: iat,ivar
  real(rk)              :: svar
  real(rk), allocatable :: tmp_icoords(:)
  real(RK), external    :: ddot
  logical               :: tok
! **********************************************************************

  ! ====================================================================
  ! CHECK FOR A CONSISTENT RESIDUE SPECIFICATION
  ! ====================================================================
  if(mod(glob%icoord,10)==3.or.mod(glob%icoord,10)==4) then
    ! The only allowed residue specification is: all atoms are 
    ! optimised, and all are part of one residue
    if(minval(glob%spec(:)) < 0) call dlf_fail("Frozen atoms are not &
        &permitted when using pure internals, use HDLC instead")
    if(minval(glob%spec(:)) /= maxval(glob%spec(:))) &
        call dlf_fail("All atoms must belong to the same residue for &
        &pure internals")
    if(maxval(glob%spec(:))==0) glob%spec(:)=1
  end if

  select case (glob%icoord)

! ======================================================================
! Cartesians
! ======================================================================
  case (0, 10)
    call dlf_cartesian_get_nivar(glob%nivar)
    if(printl >=4) then
      do iat=1,glob%nat
        select case (glob%spec(iat))
        case (0)
          write(stdout,1000) iat,"free"
        case (-1)
          write(stdout,1000) iat,"frozen"
        case (-2)
          write(stdout,1000) iat,"x frozen"
        case (-3)
          write(stdout,1000) iat,"y frozen"
        case (-4)
          write(stdout,1000) iat,"z frozen"
        case (-23)
          write(stdout,1000) iat,"x and y frozen"
        case (-24)
          write(stdout,1000) iat,"x and z frozen"
        case (-34)
          write(stdout,1000) iat,"y and z frozen"
        case default
          write(stdout,1001) iat,"free. Spec. residue:",glob%spec(iat)
        end select
      end do
    end if

    ! Extra coordinates for Lagrange-Newton
    if (glob%icoord == 10) glob%nivar = glob%nivar + 2

    ! allocate arrays
    call allocate( glob%icoords,glob%nivar)
    call allocate( glob%igradient,glob%nivar)
    call allocate( glob%step,glob%nivar)
    call allocate( glob%iweight,glob%nivar)
    glob%icoords(:)=0.D0
    glob%igradient(:)=0.D0
    glob%step(:)=0.D0
    glob%iweight(:)=0.D0

    ! Allocate conint arrays
    ! and initialise extra Lagrange-Newton coordinates
    if (glob%icoord == 10) call dlf_ln_allocate

    ! set weight
    ivar=1
    do iat=1,glob%nat
      if(glob%spec(iat)>=0) then
        glob%iweight(ivar:ivar+2)=glob%weight(iat)
        ivar=ivar+3
      else if(glob%spec(iat)==-1) then
      else if(glob%spec(iat)>=-4) then
        glob%iweight(ivar:ivar+1)=glob%weight(iat)
        ivar=ivar+2
      else
        glob%iweight(ivar)=glob%weight(iat)
        ivar=ivar+1
      end if
    end do
    ! Set some sensible weights for the extra coordinates
    if (glob%icoord == 10) then
       glob%iweight(glob%nivar - 1:glob%nivar) = 1.0d0
       ivar = ivar + 2
    endif
    if(ivar-1/=glob%nivar) then
      call dlf_fail("Error in cartesian iweight calculation")
    end if

! ======================================================================
! HDLC/DLC
! ======================================================================
  case (1:4, 11:14)

    call dlf_hdlc_init(glob%nat,glob%spec,mod(glob%icoord,10),glob%ncons, &
        glob%icons,glob%nconn,glob%iconn)

    ! calculate weights here
    call dlf_hdlc_create(glob%nat,glob%spec,glob%znuc,1,glob%xcoords, &
        glob%weight,glob%mass)

    call dlf_hdlc_get_nivar(glob%nivar)

    ! Extra coordinates for Lagrange-Newton
    if ((glob%icoord / 10) == 1) glob%nivar = glob%nivar + 2

    !allocate arrays
    call allocate( glob%icoords,glob%nivar)
    call allocate( glob%igradient,glob%nivar)
    call allocate( glob%step,glob%nivar)
    glob%icoords(:)=0.D0
    glob%igradient(:)=0.D0
    glob%step(:)=0.D0

    ! Allocate conint arrays
    ! and initialise extra Lagrange-Newton coordinates
    if ((glob%icoord / 10) == 1) call dlf_ln_allocate

    ! get weights
    call allocate( glob%iweight,glob%nivar)

    if ((glob%icoord / 10) == 1) then
       call dlf_hdlc_getweight(glob%nat, glob%nivar - 2, glob%weight, &
            glob%iweight(1:glob%nivar - 2))
       ! Set some sensible weights for the extra coordinates
       glob%iweight(glob%nivar - 1:glob%nivar) = 1.0d0
    else
       call dlf_hdlc_getweight(glob%nat,glob%nivar,glob%weight,glob%iweight)
    endif

! ======================================================================
! NEB 
! ======================================================================
  case (100:199)
    call dlf_neb_init(glob%nimage,glob%icoord)

! ======================================================================
! Dimer Method
! ======================================================================
  case (200:299)
    call dlf_dimer_init(glob%icoord)

! ======================================================================
! Wrong coordinate setting
! ======================================================================
  case default
    write(stderr,*) "Coordinate type",glob%icoord,"not implemented"
    call dlf_fail("Coordinate type error")

  end select

  ! distort if requested and a single-image method is used
  if(glob%icoord<100 .and. glob%tcoords2 .and. abs(glob%distort)>0.D0) then
    ! make xcoords2 relative
    if(.not.glob%tsrelative) glob%xcoords2(:,:,1)=glob%xcoords2(:,:,1)-glob%xcoords
    ! guess at distance distort in xyz space
    svar=dsqrt(ddot(glob%nvar,glob%xcoords2(:,:,1),1,glob%xcoords2(:,:,1),1))
    glob%xcoords2(:,:,1)=glob%xcoords2(:,:,1) / svar * glob%distort ! SIGN ENTERS HERE
    ! make xcoords2 absolute
    glob%xcoords2(:,:,1)=glob%xcoords2(:,:,1)+glob%xcoords

    ! calculate internals for coords
    glob%xgradient(:,:)=0.D0
    call dlf_direct_xtoi(glob%nvar,glob%nivar,glob%xcoords,glob%xgradient, &
        glob%icoords,glob%igradient)
    ! calculate internals for coords2 (stored in tmp_icoords)
    call allocate(tmp_icoords, glob%nivar)
    call dlf_direct_xtoi(glob%nvar,glob%nivar,glob%xcoords2(:,:,1),glob%xgradient, &
        tmp_icoords,glob%igradient)
    ! make tmp_icoords relative
    tmp_icoords=tmp_icoords-glob%icoords
    svar=dsqrt(ddot(glob%nivar,tmp_icoords,1,tmp_icoords,1))
    ! now distort icoords
    glob%icoords=glob%icoords+tmp_icoords/svar*abs(glob%distort)
    call deallocate(tmp_icoords)
    ! transform icoords back to xcoords 
    call dlf_direct_itox(glob%nvar,glob%nivar, &
        glob%icoords,glob%xcoords,tok)
    if(.not.tok) then
      call dlf_fail("Back transformation after distort failed. Use a&
          & smaller value for distort.")
    end if
  end if

! formats
1000 format ("Atom ",i6,2x,a)
1001 format ("Atom ",i6,2x,a,2x,i6)

end subroutine dlf_coords_init
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_coords_destroy
!!
!! FUNCTION
!! Deallocate global arrays belonging to internal coordinates
!!
!! SYNOPSIS
subroutine dlf_coords_destroy
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  use dlf_allocate, only: deallocate
  implicit none
! **********************************************************************
  select case (glob%icoord)

! ======================================================================
! Cartesians
! ======================================================================
  case (0, 10)

    ! deallocate arrays
    call deallocate( glob%icoords)
    call deallocate( glob%igradient)
    call deallocate( glob%step)
    call deallocate( glob%iweight)

! ======================================================================
! HDLC/DLC
! ======================================================================
  case (1:4, 11:14)

    ! deallocate arrays
    call deallocate( glob%icoords)
    call deallocate( glob%igradient)
    call deallocate( glob%step)
    call deallocate( glob%iweight)

    call dlf_hdlc_destroy

! ======================================================================
! NEB 
! ======================================================================
  case (100:199)
    call dlf_neb_destroy

! ======================================================================
! Dimer Method
! ======================================================================
  case (200:299)
    call dlf_dimer_destroy

! ======================================================================
! Wrong coordinate setting
! ======================================================================
  case default
    write(stderr,*) "Coordinate type",glob%icoord,"not implemented"
    call dlf_fail("Coordinate type error")

  end select

end subroutine dlf_coords_destroy
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_coords_xtoi
!!
!! FUNCTION
!! Transform Cartesian coordinates and gradient to internals
!!
!! In case of a direct coordinate transformation, call dlf_direct_xtoi
!!
!! In multiple image methods (NEB, dimer) call their respective xtoi routines
!!
!! INPUTS
!!
!! glob%nvar, glob%nivar, glob%xcoords, glob%xgradient
!!
!! OUTPUTS
!! 
!! glob%icoords, glob%igradient, trerun_energy, testconv
!!
!! SYNOPSIS
subroutine dlf_coords_xtoi(trerun_energy,testconv,iimage)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob, stdout, stderr, printl
  implicit none
  logical, intent(out) :: trerun_energy
  logical, intent(out) :: testconv ! is convergence tested here?
  integer, intent(out) :: iimage ! which image to be calculated next?
! **********************************************************************

  testconv=.false.
  iimage=1
  
  if(printl>=6) print*,"Transforming X to I"
  trerun_energy=.false.

  select case (glob%icoord/10)

! ======================================================================
! Direct coordinate transform
! ======================================================================
  case (0)

    call dlf_direct_xtoi(glob%nvar,glob%nivar,glob%xcoords, &
          glob%xgradient,glob%icoords,glob%igradient)

! ======================================================================
! Lagrange-Newton coordinates
! ======================================================================
  case (1)

    call dlf_ln_xtoi

! ======================================================================
! NEB 
! ======================================================================
  case (10:19) ! 100-199

    call dlf_neb_xtoi(trerun_energy,iimage)

! ======================================================================
! Dimer 
! ======================================================================
  case (20:29) ! 200-299

    call dlf_dimer_xtoi(trerun_energy,testconv)

! ======================================================================
! Wrong coordinate setting
! ======================================================================
  case default
    write(stderr,*) "Coordinate type",glob%icoord,"not implemented"
    call dlf_fail("Coordinate type error")

  end select

end subroutine dlf_coords_xtoi
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_coords_itox
!!
!! FUNCTION
!! Transform internal coordinates to xyz coordinates, handle breakdowns.
!! Gradient is not back-transformed.
!!
!! In case of a direct coordinate transformation, call dlf_direct_itox.
!!
!! In multiple image methods (NEB, dimer) call their respective itox routines.
!!
!! INPUTS
!!
!! glob%nvar, glob%nivar, glob%icoords
!!
!! OUTPUTS
!! 
!! glob%xcoords
!!
!! SYNOPSIS
subroutine dlf_coords_itox(iimage)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,stdout,printl
  implicit none
  integer, intent(out) :: iimage ! which image to be calculated next?  
  logical :: tok
! **********************************************************************
  iimage=1

  if(printl>=6) print*,"Transforming I to X"

  select case (glob%icoord/10)

! ======================================================================
! Direct coordinate transform
! ======================================================================
  case (0)

    call dlf_direct_itox(glob%nvar,glob%nivar,glob%icoords,glob%xcoords,tok)
    if(.not.tok.and. mod(glob%icoord,10)>0 .and. mod(glob%icoord,10)<=4 ) then
      if(printl>=4) write(stdout, &
          "('HDLC coordinate breakdown. Recalculating HDLCs and &
          &restarting optimiser.')")
      call dlf_hdlc_reset
      ! the arrays glob%spec,glob%znuc have to be changed when using 
      !  more instances of hdlc
      call dlf_hdlc_create(glob%nvar/3,glob%spec,glob%znuc,1, &
          glob%xcoords,glob%weight,glob%mass)
      ! recalculate iweight
      call dlf_hdlc_getweight(glob%nat,glob%nivar,glob%weight,glob%iweight)
      call dlf_formstep_restart
      glob%havehessian=.false.
    end if

! ======================================================================
! Lagrange-Newton coordinates
! ======================================================================
  case (1)

    call dlf_direct_itox(glob%nvar,glob%nivar - 2, &
         glob%icoords(1:glob%nivar - 2),glob%xcoords,tok)
    if(.not.tok.and. mod(glob%icoord,10)>0 .and. mod(glob%icoord,10)<=4 ) then
      if(printl>=4) write(stdout, &
          "('HDLC coordinate breakdown. Recalculating HDLCs and &
          &restarting optimiser.')")
      call dlf_hdlc_reset
      ! the arrays glob%spec,glob%znuc have to be changed when using 
      !  more instances of hdlc
      call dlf_hdlc_create(glob%nvar/3,glob%spec,glob%znuc,1, &
          glob%xcoords,glob%weight,glob%mass)
      ! recalculate iweight
      call dlf_hdlc_getweight(glob%nat,glob%nivar - 2, glob%weight, &
           glob%iweight(1:glob%nivar - 2))
      call dlf_formstep_restart
      glob%havehessian=.false.
    end if

! ======================================================================
! NEB 
! ======================================================================
  case (10:19) ! 100-199

    call dlf_neb_itox(iimage)

! ======================================================================
! Dimer 
! ======================================================================
  case (20:29) ! 200-299

    call dlf_dimer_itox

! ======================================================================
! Wrong coordinate setting
! ======================================================================
  case default
    write(stderr,*) "Coordinate type",glob%icoord,"not implemented"
    call dlf_fail("Coordinate type error")

  end select

end subroutine dlf_coords_itox
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_direct_xtoi
!!
!! FUNCTION
!! Transform xyz coordinates to internal coordinates
!! Only one set of coordinates is transformed, not a path 
!! (as possible by dlf_coords_xtoi)
!!
!! x-weights are transformed to i-weights
!!
!! INPUTS
!!
!! only local variables
!!
!! OUTPUTS
!! 
!! only local variables except for glob%massweight
!!
!! SYNOPSIS
subroutine dlf_direct_xtoi(nvar,nivar,xcoords,xgradient,icoords,&
    igradient)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  implicit none
  integer ,intent(in)  :: nvar,nivar
  real(rk),intent(in)  :: xcoords(nvar)
  real(rk),intent(in)  :: xgradient(nvar)
  real(rk),intent(out) :: icoords(nivar)
  real(rk),intent(out) :: igradient(nivar)
! **********************************************************************
  select case (mod(glob%icoord,10))
! ======================================================================
! Cartesian coordinates
! ======================================================================
  case (0)
    if(glob%tatoms) then
      call dlf_cartesian_xtoi(nvar/3,nivar,glob%massweight,xcoords, &
          xgradient,icoords,igradient)
    else
      icoords(:)=xcoords(:)
      igradient(:)=xgradient(:)
    end if
! ======================================================================
! HDLC/DLC
! ======================================================================
  case (1:4)
    call dlf_hdlc_xtoi(nvar/3,nivar,xcoords,xgradient,icoords,igradient)

! ======================================================================
! Wrong coordinate setting
! ======================================================================
  case default
    write(stderr,*) "Coordinate type",glob%icoord,"not implemented"
    call dlf_fail("Coordinate type error (direct)")

  end select
end subroutine dlf_direct_xtoi
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_direct_itox
!!
!! FUNCTION
!! Transform internal coordinates to xyz coordinates, do not handle breakdowns
!! Only one set of coordinates is transformed, not a path 
!! (as possible by dlf_coords_itox)
!!
!! INPUTS
!!
!! only local variables
!!
!! OUTPUTS
!! 
!! only local variables
!!
!! SYNOPSIS
subroutine dlf_direct_itox(nvar,nivar,icoords,xcoords,tok)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,stdout
  implicit none
  integer ,intent(in)   :: nvar,nivar
  real(rk),intent(in)   :: icoords(nivar)
  real(rk),intent(inout):: xcoords(nvar)
  logical ,intent(out)  :: tok
! **********************************************************************
  select case (mod(glob%icoord,10))
! ======================================================================
! Cartesian coordinates
! ======================================================================
  case (0)
    if(glob%tatoms) then
      call dlf_cartesian_itox(nvar/3,nivar,glob%massweight,icoords,xcoords)
    else
      xcoords(:)=icoords(:)
    end if
    tok=.true.

! ======================================================================
! HDLC/DLC
! ======================================================================
  case (1:4)
    call dlf_hdlc_itox(nvar/3,nivar,icoords,xcoords,tok)

! ======================================================================
! Wrong coordinate setting
! ======================================================================
  case default
    write(stderr,*) "Coordinate type",glob%icoord,"not implemented"
    call dlf_fail("Coordinate type error (direct)")

  end select
end subroutine dlf_direct_itox
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_direct_get_nivar
!!
!! FUNCTION
!! return the number of internal degrees of freedom to optimise for 
!! one image
!!
!! INPUTS
!!
!! glob%spec
!!
!! OUTPUTS
!! 
!! local
!!
!! SYNOPSIS
subroutine dlf_direct_get_nivar(nivar)
!! SOURCE
  use dlf_global, only: glob,stderr
  implicit none
  integer, intent(out) :: nivar ! number of internal variables
! **********************************************************************
  select case (mod(glob%icoord,10))
  ! Cartesian coordinates
  case (0)
    call dlf_cartesian_get_nivar(nivar)
  ! HDLC/DLC
  case (1:4)
    call dlf_hdlc_get_nivar(nivar)
! ======================================================================
! Wrong coordinate setting
! ======================================================================
  case default
    write(stderr,*) "Coordinate type",glob%icoord,"not implemented"
    call dlf_fail("Coordinate type error (direct_get_nivar)")

  end select
end subroutine dlf_direct_get_nivar
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_cartesian_get_nivar
!!
!! FUNCTION
!! return the number of degrees of freedom to optimise in Cartesian 
!! coordinates
!!
!! INPUTS
!!
!! glob%spec
!!
!! OUTPUTS
!! 
!! local
!!
!! SYNOPSIS
subroutine dlf_cartesian_get_nivar(nivar)
!! SOURCE
  use dlf_global, only: glob,stderr
  implicit none
  integer, intent(out) :: nivar ! number of internal variables
  integer :: iat
  logical :: warned
  ! negative spec values:
  ! -1  x,y,z frozen
  ! -2  x frozen
  ! -3  y frozen
  ! -4  z frozen
  ! -23 x and y frozen
  ! -24 x and z frozen
  ! -34 y and z frozen
! **********************************************************************
  warned=.false.
  if(.not.glob%tatoms) then
    ! no Cartesian constraints if no atoms as input
    nivar=glob%nvar
    return
  end if
  nivar=0
  do iat=1,glob%nat
    if(glob%spec(iat)>0) then
      if(.not.warned) then
        !print*,"Warning: fragments not used when Cartesian coordinates are requested!"
        warned=.true.
      end if
      nivar=nivar+3
    else if(glob%spec(iat)==0) then
      ! not frozen
      nivar=nivar+3
    else if(glob%spec(iat)==-1) then
    else if(glob%spec(iat)>=-4) then
      nivar=nivar+2
    else
      nivar=nivar+1
    end if
    ! invalid setting will be recognised in the conversion
  end do
end subroutine dlf_cartesian_get_nivar
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_cartesian_xtoi
!!
!! FUNCTION
!! Transform xyz coordinates to "internal" Cartesian coordinates.
!! Do mass weighting eventually
!!
!! Meaning of negative glob%spec values:
!!   -1  x,y,z frozen
!!   -2  x frozen
!!   -3  y frozen
!!   -4  z frozen
!!   -23 x and y frozen
!!   -24 x and z frozen
!!   -34 y and z frozen
!!
!! INPUTS
!!
!! glob%spec, (glob%mass)
!!
!! OUTPUTS
!! 
!! only local variables
!!
!! SYNOPSIS
subroutine dlf_cartesian_xtoi(nat,nivar,massweight,xcoords,xgradient,&
    icoords,igradient)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  implicit none
  integer ,intent(in)  :: nat,nivar
  logical ,intent(in)  :: massweight
  real(rk),intent(in)  :: xcoords(3,nat)
  real(rk),intent(in)  :: xgradient(3,nat)
  real(rk),intent(out) :: icoords(nivar)
  real(rk),intent(out) :: igradient(nivar)
  integer              :: iat,iivar
  logical              :: warned
  real(rk)             :: massf
! **********************************************************************
  if(.not.glob%tatoms) then
    call dlf_fail("dlf_cartesian_xtoi can only be used for atom input")
  end if
  warned=.false.
  iivar=1
  do iat=1,glob%nat
    if(massweight) then
      massf=dsqrt(glob%mass(iat))
    else
      massf=1.D0
    end if
    if(glob%spec(iat)>0) then
      if(.not.warned) then
        !print*,"Warning: fragments not used when Cartesian coordinates are requested!"
        warned=.true.
      end if
      icoords(iivar:iivar+2)=massf*xcoords(:,iat)
      igradient(iivar:iivar+2)=xgradient(:,iat)/massf
      iivar=iivar+3
    else if(glob%spec(iat)==0) then
      ! not frozen
      icoords(iivar:iivar+2)=massf*xcoords(:,iat)
      igradient(iivar:iivar+2)=xgradient(:,iat)/massf
      iivar=iivar+3
    else if(glob%spec(iat)==-1) then
    else if(glob%spec(iat)==-2) then
      icoords(iivar:iivar+1)=massf*xcoords(2:3,iat)
      igradient(iivar:iivar+1)=xgradient(2:3,iat)/massf
      iivar=iivar+2
    else if(glob%spec(iat)==-3) then
      icoords(iivar:iivar+1)=massf*xcoords(1:3:2,iat)
      igradient(iivar:iivar+1)=xgradient(1:3:2,iat)/massf
      iivar=iivar+2
    else if(glob%spec(iat)==-4) then
      icoords(iivar:iivar+1)=massf*xcoords(1:2,iat)
      igradient(iivar:iivar+1)=xgradient(1:2,iat)/massf
      iivar=iivar+2
    else if(glob%spec(iat)==-23) then
      icoords(iivar)=massf*xcoords(3,iat)
      igradient(iivar)=xgradient(3,iat)/massf
      iivar=iivar+1
    else if(glob%spec(iat)==-24) then
      icoords(iivar)=massf*xcoords(2,iat)
      igradient(iivar)=xgradient(2,iat)/massf
      iivar=iivar+1
    else if(glob%spec(iat)==-34) then
      icoords(iivar)=massf*xcoords(1,iat)
      igradient(iivar)=xgradient(1,iat)/massf
      iivar=iivar+1
    else
      write(stderr,"('Spec setting of atom',i5,' is wrong:',i5)") &
          iat,glob%spec(iat)
      call dlf_fail("Wrong spec setting")
    end if
  end do
  if(iivar/=nivar+1) then
    print*,iivar-1,nivar
    call dlf_fail("Error in the transformation cartesian_xtoi")
  end if
end subroutine dlf_cartesian_xtoi
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_cartesian_itox
!!
!! FUNCTION
!! Convert internal cartesian coordinats (frozen components missing)
!! to the full set of cartesians. Do mass re-weighting eventually.
!!
!! IMPORTANT: this routine relies on the fact that the frozen components 
!! of glob%xcoords are not modified by any other routine than this one.
!!
!! INPUTS
!!
!! glob%spec, (glob%mass)
!!
!! OUTPUTS
!! 
!! only local variables
!!
!! SYNOPSIS
subroutine dlf_cartesian_itox(nat,nivar,massweight,icoords,xcoords)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  implicit none
  integer ,intent(in)   :: nat,nivar
  logical ,intent(in)   :: massweight
  real(rk),intent(in)   :: icoords(nivar)
  real(rk),intent(inout):: xcoords(3,nat)
  integer               :: iat,iivar
  logical               :: warned
  real(rk)              :: massf
! **********************************************************************
  if(.not.glob%tatoms) then
    call dlf_fail("dlf_cartesian_itox can only be used for atom input")
  end if
  warned=.false.
  iivar=1
  do iat=1,glob%nat
    if(massweight) then
      massf=dsqrt(glob%mass(iat))
    else
      massf=1.D0
    end if
    if(glob%spec(iat)>0) then
      if(.not.warned) then
        !print*,"Warning: fragments not used when cartesian coordinates are requested!"
        warned=.true.
      end if
      xcoords(:,iat)=icoords(iivar:iivar+2)/massf
      iivar=iivar+3
    else if(glob%spec(iat)==0) then
      ! not frozen
      xcoords(:,iat)=icoords(iivar:iivar+2)/massf
      iivar=iivar+3
    else if(glob%spec(iat)==-1) then
    else if(glob%spec(iat)==-2) then
      xcoords(2:3,iat)=icoords(iivar:iivar+1)/massf
      iivar=iivar+2
    else if(glob%spec(iat)==-3) then
      xcoords(1:3:2,iat)=icoords(iivar:iivar+1)/massf
      iivar=iivar+2
    else if(glob%spec(iat)==-4) then
      xcoords(1:2,iat)=icoords(iivar:iivar+1)/massf
      iivar=iivar+2
    else if(glob%spec(iat)==-23) then
      xcoords(3,iat)=icoords(iivar)/massf
      iivar=iivar+1
    else if(glob%spec(iat)==-24) then
      xcoords(2,iat)=icoords(iivar)/massf
      iivar=iivar+1
    else if(glob%spec(iat)==-34) then
      xcoords(1,iat)=icoords(iivar)/massf
      iivar=iivar+1
    else
      write(stderr,"('Spec setting of atom',i5,' is wrong:',i5)") &
          iat,glob%spec(iat)
      call dlf_fail("Wrong spec setting")
    end if
  end do
  if(iivar/=nivar+1) then
    call dlf_fail("Error in the transformation cartesian_itox")
  end if
end subroutine dlf_cartesian_itox
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_coords_hessian_xtoi
!!
!! FUNCTION
!! Transform Cartesian Hessian into internal Hessian
!!
!! INPUTS
!!
!! local vars (+ glob%spec, glob%mass)
!!
!! OUTPUTS
!! 
!! glob%ihessian
!!
!! SYNOPSIS
subroutine dlf_coords_hessian_xtoi(nvar,xhessian)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  implicit none
  integer, intent(in) :: nvar
  real(rk),intent(in) :: xhessian(nvar,nvar)
! **********************************************************************
  select case (glob%icoord)
  ! Cartesians
  case (0)
    call dlf_cartesian_hessian_xtoi(nvar,glob%nivar,glob%massweight,xhessian,&
        glob%ihessian)
  ! HDLC/DLC
  case (1:4)
    call dlf_hdlc_hessian_xtoi(nvar/3,glob%nivar,glob%xcoords,xhessian,glob%ihessian)
  case default
    write(stderr,*) "Hessian transformation for coordinate type", &
        glob%icoord,"not implemented"
    call dlf_fail("Hessian transformation error")
  end select
end subroutine dlf_coords_hessian_xtoi
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_cartesian_hessian_xtoi
!!
!! FUNCTION
!! transform a Cartesian Hessian to an internal Cartesian (taking
!! frozen atoms into account. This routine may be improved for speed
!! and memory usage.
!!
!! INPUTS
!!
!! local vars (+ glob%spec, glob%mass)
!!
!! OUTPUTS
!! 
!! local vars
!!
!! SYNOPSIS
subroutine dlf_cartesian_hessian_xtoi(nvar,nivar,massweight,xhessian,&
    ihessian)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  implicit none
  integer ,intent(in)  :: nvar,nivar
  logical ,intent(in)  :: massweight
  real(rk),intent(in)  :: xhessian(nvar,nvar)
  real(rk),intent(out) :: ihessian(nivar,nivar)
  integer              :: iat,iivar,ivar,nat
  real(rk)             :: umat(nvar,nivar)
  real(rk)             :: umatt(nivar,nvar)
  real(rk)             :: tmpmat(nvar,nivar)
  real(rk)             :: massf
! **********************************************************************
  if(.not.glob%tatoms) then
    if(nvar/=nivar) call dlf_fail("Wrong number of DOF in dlf_cartesian_hessian_xtoi")
    ihessian=xhessian
    return
  end if
  ! set up umat
  nat=nvar/3
  umat=0.D0
  iivar=1
  do iat=1,glob%nat
    if(massweight) then
      massf=1.D0/dsqrt(glob%mass(iat))
    else
      massf=1.D0
    end if
    ivar=(iat-1)*3+1
    if(glob%spec(iat)>=0) then
      ! not frozen
      umat(ivar  ,iivar  )=massf
      umat(ivar+1,iivar+1)=massf
      umat(ivar+2,iivar+2)=massf
      iivar=iivar+3
    else if(glob%spec(iat)==-1) then
    else if(glob%spec(iat)==-2) then
      umat(ivar+1,iivar  )=massf
      umat(ivar+2,iivar+1)=massf
      iivar=iivar+2
    else if(glob%spec(iat)==-3) then
      umat(ivar  ,iivar  )=massf
      umat(ivar+2,iivar+1)=massf
      iivar=iivar+2
    else if(glob%spec(iat)==-4) then
      umat(ivar  ,iivar  )=massf
      umat(ivar+1,iivar+1)=massf
      iivar=iivar+2
    else if(glob%spec(iat)==-23) then
      umat(ivar+2,iivar  )=massf
      iivar=iivar+1
    else if(glob%spec(iat)==-24) then
      umat(ivar+1,iivar  )=massf
      iivar=iivar+1
    else if(glob%spec(iat)==-34) then
      umat(ivar  ,iivar  )=massf
      iivar=iivar+1
    else
      write(stderr,"('Spec setting of atom',i5,' is wrong:',i5)") &
          iat,glob%spec(iat)
      call dlf_fail("Wrong spec setting")
    end if
  end do
  if(iivar/=nivar+1) then
    print*,iivar-1,nivar
    call dlf_fail("Error in the transformation cartesian_xtoi")
  end if
  ! now transform the hessian
  call dlf_matrix_multiply(3*nat,nivar,3*nat,1.D0,xhessian,umat,0.D0,tmpmat)
  umatt=transpose(umat)
  call dlf_matrix_multiply(nivar,nivar,3*nat,1.D0,umatt,tmpmat,0.D0,ihessian)

end subroutine dlf_cartesian_hessian_xtoi
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_coords_tranrot
!!
!! FUNCTION
!! Project out the translation and rotation component of a vector in the 
!! current internal coordinate system
!! at the moment, this is only done in Cartesian coordinates if no atom 
!! is frozen
!! vector may be a gradient (which should not contain rotation or 
!! translation anyway) or the dimer axis in the dimer method
!!
!! Rotation not yet implemented!
!!
!! INPUTS
!!
!! only local variables
!!
!! OUTPUTS
!! 
!! only local variables
!!
!! SYNOPSIS
subroutine dlf_coords_tranrot(nvar,vector)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,printl,stdout
  implicit none
  integer , intent(in)       :: nvar
  real(rk), intent(inout)    :: vector(nvar)
  real(rk)                   :: svar
  integer                    :: nat
! **********************************************************************
  if(mod(glob%icoord,10)/=0) return
  if(glob%massweight) return
  if(.not.glob%tatoms) return
  if(minval(glob%spec(:)) < 0) then
    print*,"Warning: removal of rotation and translation not possible&
        & for frozen atoms"
    return
  end if
  nat=nvar/3
  svar=sum(vector(1:nvar:3))/dble(nat)
  if(printl>=6) write(stdout,'("Removing x-translation:",es12.4)') svar
  vector(1:nvar:3)=vector(1:nvar:3)-svar

  svar=sum(vector(2:nvar:3))/dble(nat)
  if(printl>=6) write(stdout,'("Removing y-translation:",es12.4)') svar
  vector(2:nvar:3)=vector(2:nvar:3)-svar

  svar=sum(vector(3:nvar:3))/dble(nat)
  if(printl>=6) write(stdout,'("Removing z-translation:",es12.4)') svar
  vector(3:nvar:3)=vector(3:nvar:3)-svar

  ! rotation not yet implemented - IMPROVE!
end subroutine dlf_coords_tranrot
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_cartesian_align
!!
!! FUNCTION
!! Translate and rotate coords2 to minimum square distance from coords1.
!! Both sets have to contain Cartesian coordinates.
!!
!! If one atom is frozen, rotation is done around this atom (no translation
!! is done).
!! If more than one atom is frozen, or other Cartesian constraints are
!! present, nothing is done.
!!
!! See W. Kabsch, Acta Cryst. A 32, p 922 (1976).
!! This follows an RMS best fit procedure.
!!
!! Attention: the input coordinates have to be Cartesians, not 
!! mass-weighted Cartesians!
!!
!! INPUTS
!!
!! only local variables
!!
!! OUTPUTS
!! 
!! only local variables
!!
!! SYNOPSIS
subroutine dlf_cartesian_align(nat,coords1,coords2)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,printl,stdout
  implicit none
  integer , intent(in)       :: nat
  real(rk), intent(inout)       :: coords1(3,nat)
  real(rk), intent(inout)    :: coords2(3,nat)
  !
  real(rk)                   :: weight(nat) ! now local, can be set on input
  real(rk)                   :: svar
  integer                    :: ivar,rotat,iat,i,j
  real(rk)                   :: center(3),rmat(3,3),rsmat(3,3)
  real(rk)                   :: eigvec(3,3),eigval(3)
  real(rk)                   :: trans(3),rotmat(3,3)
! **********************************************************************
  if(.not.glob%tatoms) return
  if(minval(glob%spec(:)) < -1) return
  rotat=0
  if(minval(glob%spec(:)) < 0) then
    ivar=0
    do iat=1,nat
      if(glob%spec(iat)==-1) then
        ivar=ivar+1
        if(ivar>1) exit
        rotat=iat
      end if
    end do
    if(ivar>1) return
  end if
  weight(:)=1.D0
  ! remove translation
  ! get centre of coords1 and translational difference
  center(:)=0.D0
  trans(:)=0.D0
  if(rotat==0) then
    do iat=1,nat
      center(:)=center(:)+weight(iat)*coords1(:,iat)
      trans(:)=trans(:)+weight(iat)*(coords1(:,iat)-coords2(:,iat))
    end do
    center(:)=center(:)/sum(weight)
    trans(:)=trans(:)/sum(weight)
  else
    center(:)=coords1(:,rotat)
    trans(:)=coords1(:,rotat)-coords2(:,rotat)
  end if
  ! translate them to common centre
  do iat=1,nat
    coords2(:,iat)=coords2(:,iat)+trans(:)
  end do
  if(printl>=6) write(stdout,"('Translating by ',3f10.5)") trans

  ! now get rotation
  ! following W. Kabsch, Acta Cryst. A 32, p 922 (1976)
  rmat=0.D0
  do iat=1,nat
    do i=1,3
      do j=1,3
        rmat(i,j)=rmat(i,j)+weight(iat)*(coords1(i,iat)-center(i))* &
            (coords2(j,iat)-center(j))
      end do
    end do
  end do
  rmat=rmat/sum(weight)
  !write(*,"('R   ',3f10.3)") rmat
  rsmat=transpose(rmat)
  eigvec=matmul(rsmat,rmat)
  rsmat=eigvec

  !write(stdout,"('RtR ',3f10.3)") rsmat
  call dlf_matrix_diagonalise(3,rsmat,eigval,eigvec)
  !do i=1,3
  !  write(*,"('Eigval, vec a_k ',f10.3,5x,3f10.3)") eigval(i),eigvec(:,i)
  !end do

  ! rsmat are the vectors b:
  j=-1
  do i=1,3
    if(eigval(i)<1.D-8) then
      if(i>1) call dlf_fail("Error in dlf_cartesian_align")
      j=1
    else
      rsmat(:,i)=1.d0/dsqrt(eigval(i)) * matmul(rmat,eigvec(:,i))
    end if
  end do
  if(j==1) then
    ! one eigenvalue was zero, the system is planar
    rsmat(1,1)=rsmat(2,2)*rsmat(3,3)-rsmat(3,2)*rsmat(2,3)
    rsmat(2,1)=rsmat(3,2)*rsmat(1,3)-rsmat(1,2)*rsmat(3,3)
    rsmat(3,1)=rsmat(1,2)*rsmat(2,3)-rsmat(2,2)*rsmat(1,3)
  end if

  do i=1,3
    do j=1,3
      rotmat(i,j)=sum(rsmat(i,:)*eigvec(j,:))
    end do
  end do
  !write(*,"('rotmat ',3f10.3)") rotmat

  do iat=1,nat
    coords2(:,iat)= coords2(:,iat)-center
    coords2(:,iat)=matmul(rotmat,coords2(:,iat))
    coords2(:,iat)= coords2(:,iat)+center
  end do
end subroutine dlf_cartesian_align
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_checkpoint_coords_write
!!
!! FUNCTION
!! Write checkpoint information
!!
!! SYNOPSIS
subroutine dlf_checkpoint_coords_write
!! SOURCE
  use dlf_global, only: glob,stderr
  implicit none
! **********************************************************************

  ! hdlc
  if(mod(glob%icoord,10)>=1 .and. mod(glob%icoord,10)<=4) then
    call dlf_checkpoint_hdlc_write
  end if

  ! NEB
  if(glob%icoord>=100 .and. glob%icoord<200) then
    call dlf_checkpoint_neb_write
  end if

  ! DIMER
  if(glob%icoord>=200 .and. glob%icoord<300) then
    call dlf_checkpoint_dimer_write
    ! make sure L-BFGS checkpoint is written if dimer is used but no
    ! global l-bfgs
    if(glob%iopt/=3) call DLF_checkpoint_LBFGS_write
  end if

end subroutine dlf_checkpoint_coords_write
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* coords/dlf_checkpoint_coords_read
!!
!! FUNCTION
!! Read checkpoint information
!!
!! SYNOPSIS
subroutine dlf_checkpoint_coords_read(tok)
!! SOURCE
  use dlf_global, only: glob,stderr
  implicit none
  logical,intent(out) :: tok
! **********************************************************************
  tok=.true.

  ! hdlc
  if(mod(glob%icoord,10)>=1 .and. mod(glob%icoord,10)<=4) then
    call dlf_checkpoint_hdlc_read(tok)
    if(.not.tok) return
  end if

  ! NEB
  if(glob%icoord>=100 .and. glob%icoord<200) then
    call dlf_checkpoint_neb_read(tok)
    if(.not.tok) return
  end if

  ! DIMER
  if(glob%icoord>=200 .and. glob%icoord<300) then
    call dlf_checkpoint_dimer_read(tok)
    if(.not.tok) return

    ! make sure L-BFGS checkpoint is written if dimer is used but no
    ! global l-bfgs
    if(glob%iopt/=3) then
      call DLF_checkpoint_LBFGS_read(tok)
    end if

  end if

end subroutine dlf_checkpoint_coords_read
!!****
