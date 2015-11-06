! **********************************************************************
! **                Optimisation Algorithms: main unit                **
! **********************************************************************
!!****h* DL-FIND/formstep
!!
!! NAME
!! formstep
!!
!! FUNCTION
!! Optimisation algorithms: determine a search direction. 
!!
!! The variable
!!              IOPT
!! determines which algorithm is to be used. The routines called by this file are
!! supposed to do unconstrained optimisation in the corresponding internal
!! coordinates.
!! The file also contains routines for calculating and updating a Hessian.
!!
!! Inputs
!!    glob%icoords
!!    glob%igradient
!!    glob%toldenergy
!!     Initialisation routines may require additional parameters.
!! 
!! Outputs
!!    glob%step
!!
!! DATA
!! $Date: 2010-05-07 19:30:43 $
!! $Rev: 390 $
!! $Author: gberan $
!! $URL: http://ccpforge.cse.rl.ac.uk/svn/dl-find/branches/release_chemsh3.3/dlf_formstep.f90 $
!! $Id: dlf_formstep.f90,v 1.1 2010-05-07 19:30:43 gberan Exp $
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
!!****
module dlf_formstep_module
  use dlf_parameter_module, only: rk
  ! conjugate gradient data
  real(rk), allocatable, save :: oldg1(:)     ! glob%nivar
  real(rk), allocatable, save :: g1(:)        ! glob%nivar
  real(rk), allocatable, save :: oldcoords(:) ! glob%nivar
  integer, save               :: cgstep       ! number of CG steps taken
  integer, parameter          :: maxcgstep=10
  ! damped dynamics data
  real(rk)                    :: fricm        ! decreasing friction in damped dynamics
  ! data for dlf_formstep_set_tsmode
  real(rk), allocatable, save :: tscoords(:) ! x-coords
  real(rk), allocatable, save :: tsmode_r(:) ! x-coords, relative
  real(rk), save              :: energy      ! ts energy
  logical , save              :: tenergy     ! is energy set?
  logical , save              :: tsc_ok      ! are TS coords ok to use?
  logical , save              :: tsm_ok      ! is TS-mode ok to use?
end module dlf_formstep_module

! module for Optimisers using the Hessian
module dlf_hessian
  use dlf_parameter_module, only: rk
  logical ,save        :: fd_hess_running ! true if FD Hessian is currently running
                                   ! initially set F in dlf_formstep_init
  integer ,save        :: nivar
  integer ,save        :: iivar,direction
  real(rk),save        :: soft ! Hessian eigenvalues absolutely smaller 
                               ! than "soft" are ignored in P-RFO
  integer ,save        :: follow ! Type of mode following:
                          ! 0: no mode following: TS mode has the lowest eigenvalue
                          ! 1: specify direction by input - not yet implemented
                          ! 2: determine direction at first P-RFO step
                          ! 3: update direction at each P-RFO step
  logical ,save        :: tsvectorset ! is tsverctor defined?
  integer ,save        :: tsmode ! number of mode to maximise
  logical ,save        :: twopoint ! type of finite difference Hessian
  real(rk),save        :: storeenergy
  logical ,save        :: carthessian ! should the Hessian be updated in 
                                      ! Cartesian coordinates (T) or internals (F)
  real(rk),allocatable,save :: eigvec(:,:)  ! (nivar,nivar) Hessian eigenmodes
  real(rk),allocatable,save :: eigval(:)    ! (nivar) Hessian eigenvalues
  real(rk),allocatable,save :: storegrad(:) ! (nivar) old gradient in fdhessian
  integer             ,save :: iupd ! actual number of Hessian updates
  real(rk),allocatable,save :: tsvector(:) ! (nivar) Vector to follow in P-RFO
  ! The old arrays are set in formstep. Used there and in hessian_update
  real(rk),allocatable,save :: oldc(:)      ! (nivar) old i-coords
  real(rk),allocatable,save :: oldgrad(:)   ! (nivar) old i-coords
  real(rk)            ,save :: minstep      ! minimum step length for Hessian update to be performed
                                       ! set in formstep_init
  logical,save              :: fracrec ! recalculate a fraction of the hessian?
end module dlf_hessian

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* formstep/dlf_formstep
!!
!! FUNCTION
!!
!! The optimisation algorithm: calculate a step vector from the gradient
!! and possibly some gradient history information or hessian.
!!
!! Some optimisers are implemented here, some in their own files (L-BFGS,
!! P-RFO).
!!
!! INPUTS
!!
!! glob%igradient
!!
!!
!! OUTPUTS
!! 
!! glob%step
!!
!! SYNOPSIS
subroutine dlf_formstep
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,stdout,printl
  use dlf_formstep_module
  use dlf_hessian
  use dlf_allocate
  implicit none
  !
  real(rk)       :: oldstep(glob%nivar)
  real(rk)       :: svar,gamma
  logical        :: trestart
  integer        :: ivar,jvar
  real(rk),allocatable :: hess(:,:) ! for Newton-Raphson
  real(rk),external :: ddot
  real(rk)       :: friction
  real(rk)       :: eigvalThreshold = 1.0d-10
  real(rk)       :: projGradient
! **********************************************************************
  select case (glob%iopt)

! ======================================================================
! Steepest descent
! ======================================================================
  case (0)

    glob%step(:)= -1.D0 * glob%igradient(:)
    if(glob%iline/=0)  glob%step(:)= glob%step(:)* glob%scalestep

! ======================================================================
! Conjugate gradient following Polak and Ribiere
! ======================================================================
  case (1)
    
    if(glob%toldenergy) then

      !Powell-Beale Restarts
      trestart= (&
          dot_product( oldg1(:) , glob%igradient(:)) .gt. &
          0.2D0 * dot_product( glob%igradient(:), glob%igradient(:)) )

      ! get old step from coordinate differences
      oldstep(:)=glob%icoords(:)-oldcoords(:)
      ! normalise?

      if(trestart) then
        gamma = 0.D0
        if(printl >= 2) print*,"Restarting CG algorithm"
      else

        ! should old gradient and old step be stored in the global module
        ! or internally here?
        gamma=dot_product( &
            glob%igradient(:)-oldg1(:) , glob%igradient(:) ) / &
            dot_product( oldg1(:) , oldg1(:) )
        if(printl>=6) print*,"Continuing CG algorithm"
      end if

      glob%step(:)= -1.D0 * glob%igradient(:) + gamma * oldstep(:)

    else
      ! first step: do an SD step
      glob%step(:)= -1.D0 * glob%igradient(:)
    end if
    oldg1(:)=glob%igradient(:)
    oldcoords(:)=glob%icoords(:)

! ======================================================================
! Conjugate gradient following Polak and Ribiere - better implementation
! ======================================================================
  case (2)
    ! oldcoords is the old search direction here!
    if(glob%toldenergy) then

      if(cgstep < maxcgstep) then

        gamma=ddot(glob%nivar,glob%igradient(:)-oldg1(:),1,&
            glob%igradient(:),1) / &
            ddot(glob%nivar,oldg1(:),1,oldg1(:),1)
        cgstep=cgstep+1
        if(gamma < 0.D0) gamma=0.D0
      else
        gamma = 0.D0
        oldcoords(:)=0.D0
        cgstep=0
        if(printl >= 4) write(stdout,'("Restarting CG algorithm")')

      end if

      glob%step(:)= -1.D0 * glob%igradient(:) + gamma * oldcoords(:)

    else
      ! first step: do an SD step
      glob%step(:)= -1.D0 * glob%igradient(:)
      cgstep=0
    end if
    
    oldg1(:)=glob%igradient(:)
    oldcoords(:)=glob%step(:)

    glob%step(:)=glob%step(:)*glob%scalestep

! ======================================================================
! L-BFGS
! ======================================================================
  CASE (3)

    CALL DLF_LBFGS_STEP(GLOB%ICOORDS,GLOB%IGRADIENT,GLOB%STEP)

! ======================================================================
! P-RFO
! ======================================================================
  CASE (10)

    ! store old coords and grad
    if(sum((oldc(:)-glob%icoords(:))**2) > minstep ) then
      oldc(:)=glob%icoords(:)
      oldgrad(:)=glob%igradient(:)
    end if
    ! send information to set_tsmode
    call dlf_formstep_set_tsmode(1,-1,glob%energy) ! send energy
    call dlf_formstep_set_tsmode(glob%nvar,0,glob%xcoords) ! TS-geometry
    call dlf_prfo_step(glob%nivar,GLOB%ICOORDS,GLOB%IGRADIENT, &
        glob%ihessian,GLOB%STEP)

! ======================================================================
! Hessian and thermal analysis only
! ======================================================================
  CASE (11)
    ! do nothing

! ======================================================================
! Newton-Raphson
! ======================================================================
  CASE (20)
    if(.not.glob%havehessian) call dlf_fail("No Hessian present in Newton-Raphson")

    ! store Hessian
    call allocate(hess,glob%nivar,glob%nivar)
    hess(:,:)=glob%ihessian(:,:)

!    ivar = array_invert(hess,svar,.true.,glob%nivar)
    call dlf_matrix_invert(glob%nivar,.true.,hess,svar)
    if(printl>=6) write(stdout,"('Determinant of H in NR step: ',es10.3)") svar

!    ! multiply inverse hessian with gradient:
!    do ivar=1,glob%nivar
!      glob%step(ivar)=-sum(hess(ivar,:)*glob%igradient(:))
!    end do
    call dlf_matrix_multiply(glob%nivar,1,glob%nivar,-1.D0,hess, &
        glob%igradient, 0.D0, glob%step)

    call deallocate(hess)
    
    ! store old coords and grad
    if(sum((oldc(:)-glob%icoords(:))**2) > minstep ) then
      oldc(:)=glob%icoords(:)
      oldgrad(:)=glob%igradient(:)
    end if

! ======================================================================
! Damped dynamics
! ======================================================================
  CASE (30)
    if(glob%toldenergy) then
      if(glob%energy<glob%oldenergy) then
        friction=fricm*glob%fricfac
        fricm=friction
      else
        if(printl>=4) write(stdout,"('Energy increasing, using high &
            &friction')")
        friction=glob%fricp
      end if
    else
      friction=glob%fric0
      fricm=friction
      oldcoords(:)=glob%icoords(:)
    end if
    glob%step(:)=1.D0/(1.D0+friction)*( &
        (1.D0-friction)*glob%icoords(:) - &
        (1.D0-friction)*oldcoords(:) - &
        glob%igradient(:)*glob%timestep**2 )

    !svar=sum(glob%step**2)*0.5D0/timestep**2+glob%energy
    !print*,"Total energy:",svar,"friction",friction

    oldcoords(:)=glob%icoords(:)

! ======================================================================
! Lagrange-Newton
! ======================================================================
  CASE (40)
    if(.not.glob%havehessian) call dlf_fail("No Hessian present in Newton-Raphson")

    ! Solving the equation for the step,
    !
    !    H s = -g
    ! =>   s = - H(-1) g
    !
    ! But matrix inversion is prone to numerical instabilities when the
    ! interstate coupling gradient is very small.
    !
    ! Therefore instead we diagonalise the Hessian matrix (following MNDO).
    !
    !      D = X(-1) H X
    !
    ! The step is then,
    !
    !      s = - X D(-1) X(-1) g
    !
    ! X, the matrix of eigenvectors, is simple to invert as X(-1) = X(t)
    ! D, the matrix of eigenvalues, is also straightforward as it is diagonal
    !
    call dlf_matrix_diagonalise(glob%nivar, glob%ihessian, eigval, eigvec)
    glob%step = 0.0d0
    do ivar = 1, glob%nivar
       if (abs(eigval(ivar)) > eigvalThreshold) then 
          projGradient = ddot(glob%nivar, eigvec(1, ivar), 1, glob%igradient, 1)
          do jvar = 1, glob%nivar
             glob%step(jvar) = glob%step(jvar) - &
                  (projGradient / eigval(ivar)) * eigvec(jvar, ivar)
          enddo
       endif
    enddo

    if (printl >= 6) then
       write(stdout, '(a)') "Hessian eigenvalues:"
       write(stdout, '(12f10.5)') eigval(:)
       write(stdout, '(a)') "Hessian eigenvectors:"
       do ivar = 1, glob%nivar
          write(stdout, '(12f10.5)') eigvec(ivar, :)
       end do
    end if

    ! store old coords and grad
    if(sum((oldc(:)-glob%icoords(:))**2) > minstep ) then
      oldc(:)=glob%icoords(:)
      call dlf_ln_savegrads
      ! Note oldgrad(:) is not used for LN, as the 'old' gradient
      ! for updating the Hessian has to be rebuilt with knowledge
      ! of the current Lagrange multipliers     
    end if

! ======================================================================
! Wrong optimisation type setting
! ======================================================================
  case default
    write(stderr,*) "Optimisation algorithm",glob%iopt,"not implemented"
    call dlf_fail("Optimisation algorithm error")

  end select
end subroutine dlf_formstep
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* formstep/dlf_formstep_restart
!!
!! FUNCTION
!!
!! Restart the optimisation algorithm.
!!
!! SYNOPSIS
subroutine dlf_formstep_restart
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  use dlf_formstep_module, only: cgstep,maxcgstep,fricm,oldcoords
  implicit none
! **********************************************************************
  select case (glob%iopt)
  case(2)
    cgstep = maxcgstep + 1 ! make sure CG is restarted
! L-BFGS
  case (3)
    call dlf_lbfgs_restart
  case (30)
    fricm=glob%fric0
    oldcoords(:)=glob%icoords(:)
  end select
  glob%toldenergy=.false.

  ! call an external routine (currently only used by ChemShell)
  call dlf_update()

end subroutine dlf_formstep_restart
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* formstep/dlf_formstep_init
!!
!! FUNCTION
!!
!! Allocate arrays for the optimisation algorithm.
!!
!!
!! SYNOPSIS
subroutine dlf_formstep_init(needhessian)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  use dlf_formstep_module
  use dlf_hessian
  use dlf_allocate, only: allocate
  implicit none
  logical,intent(out)     :: needhessian
! **********************************************************************
  fd_hess_running=.false.
  needhessian=.false.

  ! get Hessian values from global module
  soft=glob%soft
  ! To reproduce the old default "twopoint" behaviour, we default to 
  ! a two point FD Hessian if an external Hessian is selected but not
  ! available
  twopoint = (glob%inithessian == 2 .or. glob%inithessian == 0)
  carthessian=(glob%carthessian==1)
  follow=0 ! get from global eventually - IMPROVE
  tsmode=1
  minstep=glob%minstep

  ! energy not set in dlf_formstep_set_tsmode
  tenergy=.false.

  select case (glob%iopt)
! steepest descent
  case (1,2,30)
    call allocate(oldg1,glob%nivar)
    call allocate(g1,glob%nivar)
    call allocate(oldcoords,glob%nivar)
! L-BFGS
  case (3)
    call dlf_lbfgs_init(glob%nivar,glob%lbfgs_mem)
  case (10)
! P-RFO
    needhessian=.true.
    call allocate(tsvector,glob%nivar)
    tsvector(:)=0.D0
    tsvectorset=.false.
  case (11)
! Hessian and thermal analysis only
    needhessian=.true.
! Newton-Raphson
  case (20, 40)
    needhessian=.true.
  end select

  ! allocate the global Hessian if required
  if(needhessian) then
    ! part of module dlf_hessian
    nivar=glob%nivar
    call allocate(glob%ihessian,nivar,nivar)
    glob%ihessian=0.D0
    glob%havehessian=.false.

    ! update and PRFO
    call allocate(oldc,nivar)
    oldc(:)=0.D0 ! initialise, as the question of update may depend on them
    call allocate(oldgrad,nivar)
    iupd=0
    ! FD HESSIAN calculation
    call allocate(storegrad,nivar)

    call allocate(eigval,nivar)
    call allocate(eigvec,nivar,nivar)

  end if

end subroutine dlf_formstep_init
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* formstep/dlf_formstep_destroy
!!
!! FUNCTION
!!
!! Deallocate arrays for the optimisation algorithm.
!!
!!
!! SYNOPSIS
subroutine dlf_formstep_destroy
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  use dlf_formstep_module
  use dlf_hessian
  use dlf_allocate, only: deallocate
  implicit none
! **********************************************************************

  select case (glob%iopt)
! steepest descent
  case (1,2,30)
    call deallocate(oldg1)
    call deallocate(g1)
    call deallocate(oldcoords)

! L-BFGS
  case (3)
    call dlf_lbfgs_destroy
! P-RFO
  case (10)
    call deallocate(tsvector)
  end select

  if(allocated(glob%ihessian)) then
    call deallocate(glob%ihessian)
    glob%havehessian=.false.

    call deallocate(oldc)
    call deallocate(oldgrad)

    ! FD HESSIAN calculation
    call deallocate(storegrad)

    call deallocate(eigval)
    call deallocate(eigvec)

  end if
end subroutine dlf_formstep_destroy
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_checkpoint_formstep_write
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr
  use dlf_formstep_module
  use dlf_hessian
  use dlf_checkpoint, only: tchkform,write_separator
  implicit none
! **********************************************************************

! Open the checkpoint file (it may result in an empty file)
  if(tchkform) then
    open(unit=100,file="dlf_formstep.chk",form="formatted")
  else
    open(unit=100,file="dlf_formstep.chk",form="unformatted")
  end if

! Data for dlf_formstep_set_tsmode
  if(allocated(tscoords).or.allocated(tsmode_r).or.tenergy) then
    call write_separator(100,"TSMODE")

    ! write allocation status
    if(tchkform) then
      write(100,*) allocated(tscoords), allocated(tsmode_r),tenergy
    else
      write(100) allocated(tscoords), allocated(tsmode_r),tenergy
    end if

    if(allocated(tscoords)) then
      if(tchkform) then
        write(100,*) tscoords(:)
      else
        write(100) tscoords(:)
      end if
    end if

    if(allocated(tsmode_r)) then
      if(tchkform) then
        write(100,*) tsmode_r(:)
      else
        write(100) tsmode_r(:)
      end if
    end if

    if(tenergy) then
      if(tchkform) then
        write(100,*) energy
      else
        write(100) energy
      end if
    end if

    call write_separator(100,"END TSM")
  end if

  select case (glob%iopt)

! Algorithms with no checkpoint:
  case (0, 11, 20, 40, 51, 52)

! Conjugate gradient, damped dyn
  case (1:2,30)

    if(tchkform) then
      call write_separator(100,"CG-Arrays")
      write(100,*) cgstep,fricm
      write(100,*) oldg1,g1,oldcoords
      call write_separator(100,"END")
    else
      call write_separator(100,"CG-Arrays")
      write(100) cgstep,fricm
      write(100) oldg1,g1,oldcoords
      call write_separator(100,"END")
    end if

! L-BFGS
  CASE (3)

    CALL DLF_checkpoint_LBFGS_write

! P-RFO
  case (10)
    call write_separator(100,"TS-vectorset")
    if(tchkform) then
      write(100,*) tsvectorset
    else
      write(100) tsvectorset
    end if
    if(tsvectorset) then
      call write_separator(100,"TS-vector")
      if(tchkform) then
        write(100,*) tsmode,tsvector
      else
        write(100) tsmode,tsvector
      end if
    end if
    call write_separator(100,"END")

! ======================================================================
! Wrong optimisation type setting
! ======================================================================
  case default
    write(stderr,'(a,i4,a)') "Optimisation algorithm",glob%iopt,"not implemented"
    call dlf_fail("Optimisation algorithm error")

  end select
  
  ! close dlf_formstep.chk
  close(100)

! ======================================================================
! Write Hessian Data
! ======================================================================
  if(allocated(glob%ihessian)) then
    if(tchkform) then
      open(unit=100,file="dlf_hessian.chk",form="formatted")
      call write_separator(100,"Hessian size")
      write(100,*) nivar
      call write_separator(100,"Hessian data")
      write(100,*) glob%havehessian,fd_hess_running,iivar,direction,storeenergy,iupd,fracrec
      call write_separator(100,"Hessian arrays")
      write(100,*) glob%ihessian,oldc,oldgrad
      if(allocated(storegrad)) write(100,*) storegrad
      call write_separator(100,"END")
      close(100)
    else
      open(unit=100,file="dlf_hessian.chk",form="unformatted")
      call write_separator(100,"Hessian size")
      write(100) nivar
      call write_separator(100,"Hessian data")
      write(100) glob%havehessian,fd_hess_running,iivar,direction,storeenergy,iupd,fracrec
      call write_separator(100,"Hessian arrays")
      write(100) glob%ihessian,oldc,oldgrad
      if(allocated(storegrad)) write(100) storegrad
      call write_separator(100,"END")
      close(100)
    end if
  end if

end subroutine dlf_checkpoint_formstep_write

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_checkpoint_formstep_read(tok)
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stdout,stderr
  use dlf_formstep_module
  use dlf_hessian
  use dlf_allocate, only: allocate, deallocate
  use dlf_checkpoint, only: tchkform, read_separator
  implicit none
  logical,intent(out) :: tok
  logical             :: tchk
  integer             :: var
  logical             :: alloc_tscoords,alloc_tsmode_r
! **********************************************************************
  tok=.false.

  ! check if checkpoint file exists
  INQUIRE(FILE="dlf_formstep.chk",EXIST=tchk)
  if(.not.tchk) then
    write(stdout,10) "File dlf_formstep.chk not found"
    return
  end if
  
  ! open the checkpoint file
  if(tchkform) then
    open(unit=100,file="dlf_formstep.chk",form="formatted")
  else
    open(unit=100,file="dlf_formstep.chk",form="unformatted")
  end if

! Data for dlf_formstep_set_tsmode
  call read_separator(100,"TSMODE",tchk)
  if(tchk) then

    ! read allocation status and data
    if(tchkform) then
      read(100,*,end=201,err=200) alloc_tscoords,alloc_tsmode_r,tenergy
    else
      read(100,end=201,err=200) alloc_tscoords,alloc_tsmode_r,tenergy
    end if

    if(alloc_tscoords) then
      if(.not.allocated(tscoords)) call allocate(tscoords,glob%nvar)
      if(tchkform) then
        read(100,*,end=201,err=200) tscoords(:)
      else
        read(100,end=201,err=200) tscoords(:)
      end if
    else
      if(allocated(tscoords)) call deallocate(tscoords)
    end if

    if(alloc_tsmode_r) then
      if(.not.allocated(tsmode_r)) call allocate(tsmode_r,glob%nvar)
      if(tchkform) then
        read(100,*,end=201,err=200) tsmode_r(:)
      else
        read(100,end=201,err=200) tsmode_r(:)
      end if
    else
      if (allocated(tsmode_r)) call deallocate(tsmode_r)
    end if

    if(tenergy) then
      if(tchkform) then
        read(100,*,end=201,err=200) energy
      else
        read(100,end=201,err=200) energy
      end if
    end if

    call read_separator(100,"END TSM",tchk)
    if(.not.tchk) return

  else
    ! the checkpoint file has to be reopened
    if(tchkform) then
      open(unit=100,file="dlf_formstep.chk",form="formatted")
    else
      open(unit=100,file="dlf_formstep.chk",form="unformatted")
    end if
  end if

  select case (glob%iopt)

! ======================================================================
! Algorithms with no checkpoint
! ======================================================================
  case (0,11,20, 40, 51, 52)

    tok=.true.

! ======================================================================
! Algorithms with checkpoint handled here: Conjugate gradient, damped
! dynamics
! ======================================================================
  case (1:2,30)

    call read_separator(100,"CG-Arrays",tchk)
    if(.not.tchk) return 
    
    if(tchkform) then
      read(100,*,end=201,err=200) cgstep,fricm
      read(100,*,end=201,err=200) oldg1,g1,oldcoords
    else
      read(100,end=201,err=200) cgstep,fricm
      read(100,end=201,err=200) oldg1,g1,oldcoords
    end if

    call read_separator(100,"END",tchk)
    if(.not.tchk) return

! ======================================================================
! L-BFGS
! ======================================================================
  CASE (3)

    CALL DLF_checkpoint_LBFGS_read(tok)
    return

! P-RFO
  case (10)

    call read_separator(100,"TS-vectorset",tchk)
    if(.not.tchk) return 

    if(tchkform) then
      read(100,*,end=201,err=200) tsvectorset
    else
      read(100,end=201,err=200) tsvectorset
    end if

    if(tsvectorset) then
      call read_separator(100,"TS-vector",tchk)
      if(.not.tchk) return 

      if(tchkform) then
        read(100,*,end=201,err=200) tsmode,tsvector
      else
        read(100,end=201,err=200) tsmode,tsvector
      end if
    end if
    call read_separator(100,"END",tchk)
    if(.not.tchk) return

! ======================================================================
! Wrong optimisation type setting
! ======================================================================
  case default
    write(stderr,'(a,i4,a)') "Optimisation algorithm",glob%iopt,"not implemented"
    call dlf_fail("Optimisation algorithm error")

  end select

  ! close dlf_formstep.chk
  close(100)

! ======================================================================
! Read Hessian Data
! ======================================================================
  if(allocated(glob%ihessian)) then
    tok=.false.
    ! check if checkpoint file exists
    INQUIRE(FILE="dlf_hessian.chk",EXIST=tchk)
    if(.not.tchk) then
      write(stdout,10) "File dlf_hessian.chk not found"
      return
    end if
    if(tchkform) then
      open(unit=100,file="dlf_hessian.chk",form="formatted")
    else
      open(unit=100,file="dlf_hessian.chk",form="unformatted")
    end if

    call read_separator(100,"Hessian size",tchk)
    if(.not.tchk) return

    if(tchkform) then
      read(100,*,end=201,err=200) var
    else
      read(100,end=201,err=200) var
    end if
    if(var/=nivar) then
      print*,var,nivar
      write(stdout,10) "Inconsistent Hessian size"
      close(100)
      return
    end if

    call read_separator(100,"Hessian data",tchk)
    if(.not.tchk) return 

    if(tchkform) then
      read(100,*,end=201,err=200) glob%havehessian,fd_hess_running,iivar, &
          direction,storeenergy,iupd,fracrec
    else
      read(100,end=201,err=200) glob%havehessian,fd_hess_running,iivar, &
          direction,storeenergy,iupd,fracrec
    end if

    call read_separator(100,"Hessian arrays",tchk)
    if(.not.tchk) return 

    if(tchkform) then
      read(100,*,end=201,err=200) glob%ihessian,oldc,oldgrad
    else
      read(100,end=201,err=200) glob%ihessian,oldc,oldgrad
    end if
    if(allocated(storegrad)) then
      if(tchkform) then
        read(100,*,end=201,err=200) storegrad
      else
        read(100,end=201,err=200) storegrad
      end if
    end if
    call read_separator(100,"END",tchk)
    if(.not.tchk) return 

    close(100)
  end if

  tok=.true.
  return

  ! return on error
200 continue
  write(stdout,10) "Error reading CG/Hessian checkpoint file"
  return
201 continue
  write(stdout,10) "Error (EOF) reading CG/Hessian checkpoint file"
  return

10 format("Checkpoint reading WARNING: ",a)


end subroutine dlf_checkpoint_formstep_read

! **********************************************************************
! **                                                                  **
! **        DL-FIND Hessian Routines (including P-RFO)                **
! **                                                                  **
! **                                                                  **
! **                                                                  **
! **                                                                  **
! **                                                                  **
! **********************************************************************

!!****h* formstep/hessian
!!
!! NAME
!! hessian
!!
!! FUNCTION
!! Calculate (and use) the hessian. Includes P-RFO
!!
!!****


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* hessian/dlf_makehessian
!!
!! FUNCTION
!!
!! called from the main cycle in dl_find
!! Build up a Hessian by either:
!! * Reading it in from external sources (dlf_get_hessian)
!! * Updating an existing Hessian
!! * Building it from scratch by finite-difference
!! * Setting it to the identity matrix
!! * Improving an existing Hessian by finite-difference in its eigenmodes
!!
!! SYNOPSIS
subroutine dlf_makehessian(trerun_energy,tconv)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stdout,printl
  use dlf_stat, only: stat
  use dlf_hessian
  implicit none
  logical, intent(inout) :: trerun_energy
  logical, intent(inout) :: tconv
  integer                :: status
  integer                :: iimage ! not used here
  ! this should be improved ...
  real(rk) :: hess_loc(glob%nvar,glob%nvar)
  integer :: i
! **********************************************************************
  if (glob%icoord >= 10 .and. glob%icoord <= 19) then
     call dlf_conint_make_ln_hess(trerun_energy,tconv)
     return
  endif

  call dlf_hessian_update(glob%nivar, glob%icoords, oldc, glob%igradient, &
       oldgrad, glob%ihessian, glob%havehessian, fracrec)

  if(.not.glob%havehessian) then

    if(.not. fd_hess_running) then
      ! test for convergence before calculating the Hessian
      call convergence_test(stat%ccycle,.false.,tconv)
      if(tconv) return

      if (glob%inithessian == 4) then
         ! Initial Hessian is the identity matrix
         glob%ihessian = 0.0d0
         do i = 1, glob%nivar
            glob%ihessian(i, i) = 1.0d0
         end do
         glob%havehessian = .true.
      end if

      if (glob%inithessian == 0) then
         ! try to get an analytic Hessian - care about allocation business
         ! call to an external routine ...
         call dlf_get_hessian(glob%nvar,glob%xcoords,hess_loc,status)
         !switch:
         !status=1
         
         if(status==0) then
            !        write(*,"('HESS',2F10.2)") HESS_LOC
            ! convert it into internals
            call clock_start("COORDS")
            !     write(*,"('xHESS',12F10.4)") HESS_LOC
            call dlf_coords_hessian_xtoi(glob%nvar,hess_loc)
            !     write(*,"('iHESS',8F10.4)") glob%ihessian
            call clock_stop("COORDS")
            glob%havehessian=.true.
         else
            if(printl>=2) write(stdout,'(a)') &
                 "External Hessian not available, using two point FD."
         end if
      end if
    end if

    if(.not.glob%havehessian) then

      if (glob%inithessian == 3) then
         ! Simple diagonal Hessian a la MNDO
         call dlf_diaghessian(glob%nivar,glob%energy,glob%icoords, &
              glob%igradient,glob%ihessian,glob%havehessian)
      else
         ! Finite Difference Hessian calculation in internal coordinates
         call dlf_fdhessian(glob%nivar,fracrec,glob%energy,glob%icoords, &
              glob%igradient,glob%ihessian,glob%havehessian)
      end if

      ! check if FD-Hessian calculation currently running
      trerun_energy=(fd_hess_running) 
      if(trerun_energy) then
        call clock_start("COORDS")
        call dlf_coords_itox(iimage)
        call clock_stop("COORDS")
      end if

    end if

    ! tmp print JK
    if(glob%havehessian) then
      write(*,"('HESS',6F10.4)") glob%ihessian
      call dlf_matrix_diagonalise(glob%nivar,glob%ihessian,eigval,eigvec)
      write(stdout,"('Hessian eigenvalues: ',12f10.4)") eigval

    end if
  end if

end subroutine dlf_makehessian
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* hessian/dlf_diaghessian
!!
!! FUNCTION
!!
!! Calculate a simple diagonal Hessian by one-point finite difference,
!! much like the default behaviour of the standard MNDO optimiser.
!!
!! Fracrecalc is not applicable here as it makes no sense to partially update
!! an already available Hessian with a diagonal Hessian
!!
!! SYNOPSIS
subroutine dlf_diaghessian(nvar_,energy,coords,gradient,hess,havehessian)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stdout,printl
  use dlf_hessian
  implicit none
  !
  integer,intent(in)     :: nvar_ ! Number of variables 
  real(rk),intent(inout) :: energy          ! at the 2nd step out
  real(rk),intent(inout) :: coords(nvar_)   ! always changed
  real(rk),intent(inout) :: gradient(nvar_) ! at the 2nd step out
  real(rk),intent(inout) :: hess(nvar_,nvar_)
  logical,intent(out)    :: havehessian
  integer                :: i
  real(rk)               :: dx
  real(rk)               :: svar
! **********************************************************************

  if(.not.fd_hess_running) then
    ! First step - initialise
     if(printl >= 4) then
        write(stdout,'(A)') &
             "Finite-difference calculation for a diagonal Hessian"
     end if

     ! backup energy and gradient so it can be restored after Hessian calc
     storeenergy = energy
     storegrad(1:nvar_) = gradient(:)

     hess = 0.0d0
     ! move to the 2nd point, backwards along the gradient
     do i = 1, nvar_
        coords(i) = coords(i) - sign(glob%delta, gradient(i))
     end do

     fd_hess_running = .true.
     havehessian = .false.
   
     if(printl >= 4) then
        write(stdout,'("Delta : ",es10.2)') glob%delta
     end if

     return
  end if ! End of initialisation step


  ! Second step - calculate diagonal Hessian
  if(printl>=4) then
    write(stdout,"('Finite-difference diagonal Hessian: 2nd point')")
    write(stdout,"('Energy difference to 1st point:        ',es10.2,' H')") energy-storeenergy
    svar=sqrt(sum((gradient(:)-storegrad(1:nvar_))**2))
    write(stdout,"('Abs. Gradient difference to 1st point: ',es10.2)") svar
  end if

  ! calculate Hessian and restore coordinates
  do i = 1, nvar_
     dx = sign(glob%delta, storegrad(i))
     hess(i, i) = (storegrad(i) - gradient(i)) / dx
     ! Hessian should be positive
     if (hess(i, i) < 0.0d0) hess(i, i) = abs(storegrad(i)) / glob%delta
     ! Minimum threshold set to identity matrix element
     hess(i, i) = max(hess(i, i), 1.0d0)
     coords(i) = coords(i) + dx   
  end do

  ! restore energy and gradient
  energy = storeenergy
  gradient(:) = storegrad(1:nvar_)
  
  fd_hess_running=.false.
  havehessian=.true.

end subroutine dlf_diaghessian
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* hessian/dlf_fdhessian
!!
!! FUNCTION
!!
!! Calculate the Hessian by finite differencing.
!!
!! One-point and Two-point formulas are available.
!!
!! If Fracrecalc is true, only part of the Hessian is recalculated:
!! The Hessian is expressed in its eigenmodes. Finite-difference elongations
!! along the lowest eigenmodes are calculated. The gradient is transformed
!! back into the eigenmode-space. Finally, the Hessian is transformed back
!! into coordinate space. Frac_recalc is only used if a Hessian is
!! already available, and if do_partial_fd (in _update) is true.
!!
!! SYNOPSIS
subroutine dlf_fdhessian(nvar_,fracrecalc,energy,coords,gradient,hess,havehessian)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stdout,printl
  use dlf_allocate, only: allocate, deallocate
  use dlf_hessian
  implicit none
  !
  integer,intent(in)     :: nvar_ ! Number of variables 
  logical,intent(inout)  :: fracrecalc ! Recalculate only fraction
  real(rk),intent(inout) :: energy          ! at the last step out
  real(rk),intent(inout) :: coords(nvar_)   ! always changed
  real(rk),intent(inout) :: gradient(nvar_) ! at the last step out
  real(rk),intent(inout) :: hess(nvar_,nvar_)
  logical,intent(out)    :: havehessian
  integer                :: fail,ivar,jvar
  real(rk)               :: svar
  real(rk),allocatable   :: tmpmat(:,:)
! **********************************************************************
  if(.not.fd_hess_running) then
    ! First step - initialise
    if(fracrecalc) then

      ! all negative, "zero", and the first positive (i.e. > soft) modes
      ! will be recalculated.

      ! Diagonalise the present Hessian
      call dlf_matrix_diagonalise(nvar_,hess,eigval(1:nvar_),eigvec(1:nvar_,1:nvar_))
      !write(stdout,"('eigval-pa',12f10.4)") eigval

      do ivar=1,nvar_
        if(eigval(ivar)>max(soft,0.D0)) exit
      end do

      if(printl >= 4) then
        write(stdout,'(A,i3,a)') "Finite-difference calculation of the", &
            ivar," lowest Hessian eigenmodes"
      end if
      nivar=ivar
    else
      if(printl >= 4) then
        write(stdout,'(A)') &
            "Finite-difference calculation of the whole Hessian"
      end if
      nivar=nvar_
    end if
    iivar=1
    direction=1 ! do first 1, then -1
    storeenergy=energy

    ! keep gradient and restore it when Hessian has finished
    storegrad(1:nvar_)=gradient(:)

    fd_hess_running=.true.

    if(fracrecalc) then
      hess=0.D0
      do ivar=1,nvar_
        hess(ivar,ivar)=eigval(ivar)
      end do

      ! do the step
      coords(:)=coords(:) + glob%delta * eigvec(1:nvar_,1)

    else
      hess=0.D0

      ! do the step
      coords(1)=coords(1)+glob%delta

    end if

    

    havehessian=.false.

    if(printl >= 4) then
      write(stdout,'("Delta : ",es10.2)') glob%delta
    end if

    return
    
  end if !(.not.fd_hess_running)

  if(printl>=4) then
    write(stdout,"('Finite-difference Hessian: variable ',i4,'/',i4,&
        &' direction=',i2)") iivar,nivar,direction
    write(stdout,"('Energy difference to midpoint:        ',es10.2,' H')") energy-storeenergy
    svar=sqrt(sum((gradient(:)-storegrad(1:nvar_))**2))
    write(stdout,"('Abs. Gradient difference to midpoint: ',es10.2)") svar
  end if

  ! The values of iivar and direction are those for which the gradient is currently available

  if(iivar<nivar.or. (twopoint.and.direction==1) ) then
    ! general step in the course of Hessian calculation

    if(direction==1.and.twopoint) then
      hess(iivar,:)=gradient(:)
      direction=-1
      if(fracrecalc) then
        coords(:)=coords(:) - 2.D0 * glob%delta * eigvec(1:nvar_,iivar)
      else
        coords(iivar)=coords(iivar)-2.D0*glob%delta
      end if
    else
      if(twopoint) then
        !set back coordinates
        if(fracrecalc) then
          gradient(:)=(hess(iivar,:)-gradient(:))/(2.D0*glob%delta)
          call dlf_matrix_multiply(1,nvar_,nvar_,1.D0,gradient,eigvec(1:nvar_,1:nvar_), &
               0.D0,hess(iivar,:))
          coords(:)=coords(:) + glob%delta * eigvec(1:nvar_,iivar)
        else
          hess(iivar,:)=(hess(iivar,:)-gradient(:))/(2.D0*glob%delta)
          coords(iivar)=coords(iivar)+glob%delta
        end if
        direction=1
      else
        if(fracrecalc) then
          gradient(:)=(gradient(:)-storegrad(1:nvar_)) / glob%delta
          call dlf_matrix_multiply(1,nvar_,nvar_,1.D0,gradient,eigvec(1:nvar_,1:nvar_), &
               0.D0,hess(iivar,:))
          coords(:)=coords(:) - glob%delta * eigvec(1:nvar_,iivar)
        else
          hess(iivar,:)=(gradient(:)-storegrad(1:nvar_)) / glob%delta
          coords(iivar)=coords(iivar)-glob%delta
        end if
      end if
      iivar=iivar+1
      if(fracrecalc) then
        coords(:)=coords(:) + glob%delta * eigvec(1:nvar_,iivar)
      else
        coords(iivar)=coords(iivar)+glob%delta
      end if
    end if
      
    havehessian=.false.
      
  else
    ! final step: calculate Hessian and deallocate
    !set back coordinates
    if(twopoint) then
      if(fracrecalc) then
        gradient(:)=(hess(iivar,:)-gradient(:))/(2.D0*glob%delta)
        call dlf_matrix_multiply(1,nvar_,nvar_,1.D0,gradient,eigvec(1:nvar_,1:nvar_), &
             0.D0,hess(iivar,:))
        coords(:)=coords(:) + glob%delta * eigvec(1:nvar_,iivar)
      else
        hess(iivar,:)=(hess(iivar,:)-gradient(:))/(2.D0*glob%delta)
        coords(iivar)=coords(iivar)+glob%delta
      end if
    else
      if(fracrecalc) then
        gradient(:)=(gradient(:)-storegrad(1:nvar_)) / glob%delta
        call dlf_matrix_multiply(1,nvar_,nvar_,1.D0,gradient,eigvec(1:nvar_,1:nvar_), &
             0.D0,hess(iivar,:))
        coords(:)=coords(:) - glob%delta * eigvec(1:nvar_,iivar)
      else
        hess(iivar,:)=(gradient(:)-storegrad(1:nvar_)) / glob%delta
        coords(iivar)=coords(iivar)-glob%delta
      end if
    end if

    ! Symmetrise the Hessian
    if(fracrecalc) then

      ! build up the Hessian

      ! symmetrise
      do ivar=1,nivar
        do jvar=ivar+1,nvar_
          if(jvar<=nivar) then
            hess(ivar,jvar)=0.5D0*(hess(ivar,jvar)+hess(jvar,ivar))
            hess(jvar,ivar)=hess(ivar,jvar)
          else
            hess(jvar,ivar)=hess(ivar,jvar)
          end if
        end do
      end do

      ! Multiply with the eigenvectors to restore the real Hessian
      call allocate(tmpmat,nvar_,nvar_)
      eigvec(1:nvar_,1:nvar_) = transpose(eigvec(1:nvar_,1:nvar_))
      call dlf_matrix_multiply(nvar_,nvar_,nvar_,1.D0,hess,eigvec(1:nvar_,1:nvar_),0.D0,tmpmat)
      eigvec(1:nvar_,1:nvar_) = transpose(eigvec(1:nvar_,1:nvar_))
      call dlf_matrix_multiply(nvar_,nvar_,nvar_,1.D0,eigvec(1:nvar_,1:nvar_),tmpmat,0.D0,hess)
      call deallocate(tmpmat)
      
    else
      do ivar=1,nivar
        do jvar=ivar+1,nivar
          hess(ivar,jvar)=0.5D0*(hess(ivar,jvar)+hess(jvar,ivar))
          hess(jvar,ivar)=hess(ivar,jvar)
        end do
      end do
    end if

    ! restore energy and gradient
    energy=storeenergy
    gradient(:)=storegrad(1:nvar_)

    fd_hess_running=.false.

    havehessian=.true.

  end if

end subroutine dlf_fdhessian
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* hessian/dlf_hessian_update
!!
!! FUNCTION
!!
!! Update the global Hessian
!!
!! If we do not have a Hessian, there is nothing to update. Return.
!!
!! If we have one, and are allowed to update it, update it. If too many
!!  updates have been made, destroy and recalculate it.
!!
!! COMMENTS
!!
!! The gradient, coordinates and hessian are passed to and from the
!! routine through the arguments and are therefore independent of the
!! data stored in glob. This is useful for the Lagrange-Newton method
!! where the optimiser gradient is not necessarily the same as the 
!! gradient used for the update of the Hessian.
!!
!! INPUTS
!!
!! nvar - dimensions of gradient, coordinates, Hessian
!! coords(nvar) - current coordinates
!! oldcoords(nvar) - coordinates of previous step
!! gradient(nvar) - current gradient
!! oldgradient(nvar) - gradient of previous step
!! havehessian - if false return as no update is possible
!! fracrecalc
!! glob%update - specifies which update algorithm to use
!! glob%maxupd - check that max no. of updates is not exceeded
!! glob%nivar - check that nvar is consistent
!! glob%icoord - for check of consistency of nvar
!! (dlf_hessian) fd_hess_running
!! (dlf_hessian) iupd - no. of updates since last reset
!!
!! OUTPUTS
!! 
!! hess(nvar, nvar) - the updated Hessian in internal coordinates
!! havehessian - false if no update requested or max updates reached
!! fracrecalc
!! (dlf_hessian) iupd - no. of updates since last reset
!!
!! SYNOPSIS
subroutine dlf_hessian_update(nvar, coords, oldcoords, gradient, &
     oldgradient, hess, havehessian, fracrecalc)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stdout, stderr, printl
  use dlf_hessian
  implicit none
  integer,intent(in)    :: nvar ! used for temporary storage arrays
  real(rk),intent(in) :: coords(nvar)   
  real(rk),intent(in) :: oldcoords(nvar)
  real(rk),intent(in) :: gradient(nvar)
  real(rk),intent(in) :: oldgradient(nvar)
  real(rk),intent(inout) :: hess(nvar,nvar)
  logical,intent(inout) :: havehessian
  logical,intent(inout) :: fracrecalc ! recalculate fraction of Hessian
  ! temporary arrays
  real(rk) :: fvec(nvar),step(nvar), tvec(nvar)
  real(rk) :: fx_xx,xx,svar,bof, dds, ddtd
  real(RK) ,external :: ddot
  integer  :: ivar,jvar
  logical,parameter :: do_partial_fd=.true. ! Current main switch!
! **********************************************************************
  if(.not.fd_hess_running) fracrecalc=.false.
  if(.not.havehessian.or.glob%update==0) then
    havehessian=.false.
    return
  end if

  ! Check for maximum number of updates reached
  ! In case of partial finite-difference, do an update first, then return and 
  ! Recalculate the lower modes
  if(iupd>=glob%maxupd .and. .not. do_partial_fd) then
    havehessian=.false.
    iupd=0
    hess = -1.D0
    fracrecalc=.false.
    return
  end if

  if (glob%icoord >= 10 .and. glob%icoord <= 19) then
     if (nvar /= glob%nivar - 2) &
          call dlf_fail("Inconsistent Lagrange-Newton nvar in dlf_hessian_update")
  elseif(nvar/=glob%nivar) then
      call dlf_fail("Inconsistent nvar in dlf_hessian_update")
  endif

  ! Useful variables for updating
  fvec(:) = gradient(:) - oldgradient(:)
  step(:) = coords(:) - oldcoords(:)

  xx=ddot(nvar,step,1,step,1)
  if(xx <= minstep ) then
    if(printl>=2) write(stdout,"('Step too small. Skipping hessian update')")
    return
  end if

  iupd=iupd+1

  if(printl>=2) then
     select case (glob%update)
     case(1)
        write(stdout,"('Updating Hessian with the Powell update, No ',i5)") iupd
     case(2)
        write(stdout,"('Updating Hessian with the Bofill update, No ',i5)") iupd
     case(3)
        write(stdout,"('Updating Hessian with the BFGS update, No ',i5)") iupd
     end select
  end if
    
  select case (glob%update)
  case(1,2)
     ! Powell/Bofill updates

     ! fvec = fvec - hessian x step
     call dlf_matrix_multiply(nvar,1,nvar,-1.D0, hess,step,1.D0,fvec)

     fx_xx=ddot(nvar,fvec,1,step,1) / xx

     if(glob%update==2) then
        svar=ddot(nvar,fvec,1,fvec,1)
        bof=fx_xx**2 * xx / svar
        if(printl>=6) write(stdout,'("Bof=",es10.3)') bof
     end if

     do ivar=1,nvar
        do jvar=ivar,nvar

           ! Powell
           svar=fvec(ivar)*step(jvar) + step(ivar)*fvec(jvar) - fx_xx* step(ivar)*step(jvar)

           if(glob%update==2) then
              ! Bofill
              svar=svar*(1.D0-bof) + bof/fx_xx * fvec(ivar)*fvec(jvar)
           end if

           hess(ivar,jvar) = hess(ivar,jvar) + svar/xx
           hess(jvar,ivar) = hess(ivar,jvar)
        end do
     end do
     
  case(3)
     ! BFGS update

     ! tvec is hessian x step
     tvec = 0.0d0
     call dlf_matrix_multiply(nvar, 1, nvar, 1.0d0, hess, step, 0.0d0, tvec)

     dds = ddot(nvar, fvec, 1, step, 1)
     ddtd = ddot(nvar, step, 1, tvec, 1)

     do ivar = 1, nvar
        do jvar = ivar, nvar
           svar = (fvec(ivar) * fvec(jvar)) / dds - &
                  (tvec(ivar) * tvec(jvar)) / ddtd
           hess(ivar, jvar) = hess(ivar, jvar) + svar
           hess(jvar, ivar) = hess(ivar, jvar)
        end do
     end do

  case default
     ! Update mechanism not recognised
     write(stderr,*) "Hessian update", glob%update, "not implemented"
     call dlf_fail("Hessian update error")

  end select

  ! Check for maximum number of updates reached
  ! In case of partial finite-difference, do an update first, then return and 
  ! Recalculate the lower modes
  if(iupd>=glob%maxupd .and. do_partial_fd) then
    havehessian=.false.
    iupd=0
    fracrecalc=.true.

  end if

end subroutine dlf_hessian_update
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* hessian/dlf_prfo_step
!!
!! FUNCTION
!!
!! Calculate a P-RFO step
!! This routine does only require the current Hessian and gradient,
!! not an old one.
!!
!! Notes to the Hessian update:
!! Paul updates some modes of the Hessian, whenever the number of positive 
!! eigenvalues (non-positive is below "soft") changes. He recalculates all
!! negative and soft modes plus one.
!!
!! SYNOPSIS
subroutine dlf_prfo_step(nvar,coords,gradient,hessian,step)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: stdout,printl,stderr,glob,pi
  use dlf_hessian
  implicit none
  integer ,intent(in)   :: nvar 
  real(rk),intent(in)   :: coords(nvar)
  real(rk),intent(in)   :: gradient(nvar)
  real(rk),intent(in)   :: hessian(nvar,nvar)
  real(rk),intent(out)  :: step(nvar)
  !
  ! these may be made variable in the future:
  integer   :: maxit=100 ! maximum number of iterations to find lambda
  real(rk)  :: tol=1.D-5 ! tolerance in iterations
  real(rk)  :: bracket=1.D-4  ! threshold for using bracket scheme - why not 0.D0?
  real(rk)  :: delta=0.05D0, big=1000.D0 ! for bracketing 
  !
  real(rk)              :: lamts,lamu,laml,lam
  real(rk)              :: ug(nvar) ! eigvec * gradient
  real(rk)              :: tmpvec(nvar)
  integer               :: ivar,iter,nmode
  real(rk)              :: svar,lowev,maxov
  real(rk)              :: soft_tmp,ev2(nvar)
  real(RK) ,external    :: ddot
  logical               :: conv=.false.,err=.false.
  logical :: dbg=.true.
! **********************************************************************

  call dlf_matrix_diagonalise(nvar,hessian,eigval,eigvec)

  if(printl >= 2) then
    write(stdout,"('Eigenvalues of the Hessian:')") 
    write(stdout,"(9f10.4)") eigval
  end if

  ! Determine mode to follow
  if(follow==0) then
    tsmode=1
  else if(follow==1) then
    call dlf_fail("Hessian mode following 1 not implemented")
  else if(follow==2.or.follow==3) then
    if(tsvectorset) then
      maxov= dabs(ddot(nvar,eigvec(:,tsmode),1,tsvector,1))
      nmode=tsmode
      if(printl>=4) write(stdout,"('Overlap of current TS mode with &
          &previous one ',f6.3)") maxov
      do ivar=1,nvar
        if(ivar==tsmode) cycle
        svar= dabs(ddot(nvar,eigvec(:,ivar),1,tsvector,1))
        if(svar > maxov) then
          if(printl>=6) write(stdout,"('Overlap of mode',i4,' with &
              & TS-vector is ',f6.3,', larger than TS-mode',f6.3)") &
              ivar,svar,maxov
          maxov=svar
          nmode=ivar
        end if
      end do
      if(nmode /= tsmode) then
        !mode switching!
        if(printl>=2) write(stdout,"('Switching TS mode from mode',i4,&
            &' to mode',i4)") tsmode,nmode
        tsmode=nmode
        if(printl>=4) write(stdout,"('Overlap of current TS mode with &
            &previous one ',f6.3)") maxov
      end if
      if(follow==3) tsvector=eigvec(:,tsmode)
    else
      ! first step: use vector 1
      tsvector(:)=eigvec(:,1)
      tsvectorset=.true.
    end if
  else
    write(stderr,"('Wrong setting of follow:',i5)") follow
    call dlf_fail("Hessian mode following wrong")
  end if

  ! print frequency in case of mass-weighted coordinates
  if(glob%massweight .and. printl>=2) then
    ! sqrt(H/u)/a_B/2/pi/c / 100
    svar=sqrt( 4.35974417D-18/ 1.66053886D-27 ) / ( 2.D0 * pi * &
        0.5291772108D-10 * 299792458.D0) / 100.D0
    svar=sqrt(abs(eigval(tsmode))) * svar
    if(eigval(tsmode)<0.D0) svar=-svar
    write(stdout,"('Frequency of transition mode',f10.3,' cm^-1 &
        &(negative value denotes imaginary frequency)')") &
        svar
  end if

  call dlf_formstep_set_tsmode(nvar,11,eigvec(:,tsmode))

  ! calculate eigvec*gradient
  do ivar=1,nvar
    ug(ivar) = ddot(nvar,eigvec(:,ivar),1,gradient(:),1)
  end do

  ! calculate Lambda that minimises along the TS-mode:
  lamts=0.5D0 * ( eigval(tsmode) + dsqrt( eigval(tsmode)**2 + 4.D0 * ug(tsmode)**2) )

  if(printl >= 2 .or. dbg) then
    write(stdout,'("Lambda for maximising TS mode:     ",es12.4," Eigenvalue:",es12.4)') lamts,eigval(tsmode)
  end if

  ! Calculate the number of modes considered "soft"
  nmode=0
  soft_tmp=soft
  ev2=0.D0
  do ivar=1,nvar
    if(ivar==tsmode) cycle
    if(abs(eigval(ivar)) < soft ) then
      nmode=nmode+1
      ev2(ivar)=dabs(eigval(ivar))
    end if
  end do

  ! Check that at most 6 modes are considered "soft"
  if(nmode>6) then
    do ivar=nmode-1,6,-1
      soft_tmp=maxval(ev2)
      ev2(maxloc(ev2))=0.D0
    end do
    if(printl>=4) write(stdout,'("Criterion for soft modes tightened to &
        &",es12.4)') soft_tmp
    ! recalculate nmode
    nmode=0
    do ivar=1,nvar
      if(ivar==tsmode) cycle
      if(abs(eigval(ivar)) < soft_tmp ) nmode=nmode+1
    end do
  end if

  if(nmode>0.and.printl>=2) &
      write(stdout,'("Ignoring ",i3," soft modes")') nmode

  ! find lowest eigenvalue that is not TS-mode and not soft
  !   i.e. the lowest eigenmode that is minimised
  do ivar=1,nvar
    if(ivar==tsmode) cycle
    if(abs(eigval(ivar)) < soft_tmp ) cycle
    lowev=eigval(ivar)
    exit
  end do

  lamu=0.D0
  laml=0.D0
  lam=0.D0

  if(lowev < bracket) then
    lam=lowev-delta
    lamu=lowev
    laml=-big
  end if

  do iter=1,maxit

    svar=0.D0

    do ivar=1,nvar
      if(ivar==tsmode) cycle
      if(abs(eigval(ivar)) < soft_tmp ) cycle
      if(dabs(lam - eigval(ivar)) > 1.D-14) &
          svar=svar+ ug(ivar)**2 / (lam - eigval(ivar) )
    end do

    if(abs(svar-lam) < tol) then
      ! we are converged

      if(lam>lowev) then
        print*,"Lambda > lowest non-TS eigenvalue, bad Hessian?"
        err=.true.
      end if

      if(lam>0.D0 .and. lowev>0.D0) then
        print*,"Lambda and lowest non-TS eigenvalue >0. Bad Hessian?"
        print*,"Lambda:",lam
        !err=.true.
      end if

      if(dbg.and..not.err) then
        print*,"Lambda converged in",iter,"iterations"
      end if
      
      conv=.true.

      exit

    end if

    ! we are not converged. Next iteration:

    !write(*,'("A",4f15.8)') svar,lam,lamu,laml
    if(lowev < bracket ) then
      if(svar < lam) lamu=lam
      if(svar > lam) laml=lam
      if(laml > -big) then
        lam=0.5D0 * (lamu + laml)
      else
        lam = lam-delta
      end if
    else
      lam=svar
    end if
    !write(*,'("B",4f15.8)') svar,lam,lamu,laml

  end do

  !if(.not.conv) call dlf_fail("P-RFO loop not converged")
  if(err) call dlf_fail("P-RFO error")

  if(printl >= 2 .or. dbg) &
      write(stdout,'("Lambda for minimising other modes: ",es12.4)') lam

  ! calculate step:
  step=0.D0
  do ivar=1,nvar
    if(ivar==tsmode) then
      if( abs(lamts-eigval(ivar)) < 1.D-5 ) then
        ug(ivar)=1.D0
      else
        ug(ivar)=ug(ivar) / (eigval(ivar) - lamts)
      end if
    else
      if(abs(eigval(ivar)) < soft_tmp ) then
        if(printl>=6) write(stdout,'("Mode ",i4," ignored, as &
            &|eigenvalue| ",es10.3," < soft =",es10.3)') &
            ivar,eigval(ivar),soft_tmp
        cycle
      end if
      if( abs(lam-eigval(ivar)) < 1.D-5 ) then
        print*,"WARNING: lam-eigval(ivar) small for non-TS mode",ivar,"!"
        ug(ivar)= -1.D0 / eigval(ivar) ! take a newton-raphson step for this one
      else
        ug(ivar)=ug(ivar) / (eigval(ivar) - lam)
      end if
    end if
    if(printl>=6) write(stdout,'("Mode ",i4," Length ",es10.3)') ivar,ug(ivar)
    step(:) = step(:) - ug(ivar) * eigvec(:,ivar)
  end do

  if(printl >= 2 .or. dbg) &
      write(stdout,'("P-RFO step length:                 ",es12.4)') sqrt(sum(step(:)**2))

end subroutine dlf_prfo_step
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* formstep/dlf_formstep_set_tsmode
!!
!! FUNCTION
!!
!! Write the transition mode to a file
!!
!! mode:
!!   -2   deallocate all arrays
!!   -1   coords contains only energy (nvar=1)
!!   00   TS structure in x-coordinates
!!   01   TS structure in i-coordinates (similar x-coords have to be provided priorly)
!!   10   TS mode relative to TS structure in x-coordinates
!!   11   TS mode relative to TS structure in i-coordinates
!!   20   TS mode absolute in x-coordinates (not yet implemented)
!!   21   TS mode absolute in i-coordinates (not yet implemented)
!!
!! This routine does not produce errors in case of wrong input,
!! it just does nothing in this case.
!!
!! SYNOPSIS
subroutine dlf_formstep_set_tsmode(nvar,mode,coords)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,stdout,printl,printf
  use dlf_formstep_module, only: tscoords, tsmode_r, energy, tenergy, &
      tsc_ok, tsm_ok
  use dlf_allocate, only: allocate, deallocate
  implicit none
  integer ,intent(in) :: nvar
  integer ,intent(in) :: mode
  real(rk),intent(in) :: coords(nvar)
  real(rk),allocatable :: store1(:),store2(:),store3(:)
  logical              :: tok
  real(rk)             :: svar
  real(rk),external    :: ddot
! **********************************************************************
  if(mode==-1) then
    if(nvar/=1) return
    tenergy=.true.
    energy=coords(1)
  else if(mode==0) then
    if(allocated(tscoords)) call deallocate(tscoords)
    if(nvar/=glob%nvar) return
    call allocate(tscoords,nvar)
    tscoords=coords
    tsc_ok=.true.
  else if(mode==1) then
    if(.not.(allocated(tscoords).and.tsc_ok)) return
    call dlf_direct_itox(glob%nvar,nvar,coords,tscoords,tok)
    if(.not.tok) then
      call deallocate(tscoords)
      write(stdout,'(a)') "HDLC breakdown in set_tsmode ignored in mode 1"
      tsc_ok=.false.
      return
    end if
  else if(mode==10) then
    if(allocated(tsmode_r)) call deallocate(tsmode_r)
    if(nvar/=glob%nvar) return
    call allocate(tsmode_r,nvar)
    tsmode_r=coords
    tsm_ok=.true.
  else if(mode==11) then
    ! transform tscoords to icoords
    if(.not.(allocated(tscoords).and.tsc_ok)) return
    if(allocated(tsmode_r)) call deallocate(tsmode_r)
    ! nvar= number of internal coordinates
    call allocate(tsmode_r,glob%nvar)
    tsmode_r=tscoords
    call allocate(store1,glob%nvar) ! xgradient
    store1=0.D0
    call allocate(store2,nvar) ! igradient
    store2=0.D0
    call allocate(store3,nvar) ! icoords
    store3=0.D0
    call dlf_direct_xtoi(glob%nvar,nvar,tscoords,store1,store3,store2)
    ! add coords
    ! make sure the relative coords are short:
    svar=ddot(nvar,coords,1,coords,1)
    store3=store3+coords/sqrt(svar)*0.05D0
    ! transform to x-coords
    call dlf_direct_itox(glob%nvar,nvar,store3,tsmode_r,tok)
    call deallocate(store1)
    call deallocate(store2)
    call deallocate(store3)
    if(.not.tok) then
      write(stdout,'(a)') "HDLC breakdown in set_tsmode ignored in mode 11"
      call deallocate(tsmode_r)
      tsm_ok=.false.
      return
    end if
    ! subtract tscoords
    tsmode_r=tsmode_r-tscoords
    tsm_ok=.true.
  else if(mode==-2) then
    if(allocated(tsmode_r)) call deallocate(tsmode_r)
    if(allocated(tscoords)) call deallocate(tscoords)
    tsc_ok=.false.
    tsm_ok=.false.
  end if
  if(tsc_ok.and.tsm_ok.and. tenergy) then
    if(printf>=2) then
      if(printl>=4) write(stdout,"('Writing TS-mode')")
      call dlf_put_coords(glob%nvar,1,energy,tscoords,glob%iam)
      if(glob%tsrelative) then
        call dlf_put_coords(glob%nvar,2,energy,tsmode_r,glob%iam) 
      else
        call dlf_put_coords(glob%nvar,2,energy,tscoords+tsmode_r,glob%iam) 
      end if
    end if
    !call deallocate(tscoords)
    !call deallocate(tsmode_r)
    tsc_ok=.false.
    tsm_ok=.false.
    tenergy=.false.
  end if

end subroutine dlf_formstep_set_tsmode
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* formstep/dlf_formstep_get_ra
!!
!! FUNCTION
!!
!! Get a real-number array from the formstep module
!!
!! SYNOPSIS
subroutine dlf_formstep_get_ra(label,array_size,array,tok)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_formstep_module, only: tscoords, tsmode_r
  implicit none
  character(*), intent(in) :: label
  integer     , intent(in) :: array_size
  real(rk)    , intent(out):: array(array_size)
  logical     , intent(out):: tok
! **********************************************************************
  tok=.false.
  if (label=="TSCOORDS") then
    if(.not.allocated(tscoords)) return
    if(size(tscoords) /= array_size) return
    array(:)=tscoords(:)
    tok=.true.
  else if (label=="TSMODE_R") then
    if(.not.allocated(tsmode_r)) return
    if(size(tsmode_r) /= array_size) return
    array(:)=tsmode_r(:)
    tok=.true.
  else
    call dlf_fail("Wrong label in dlf_formstep_get_ra")
  end if
end subroutine dlf_formstep_get_ra
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* hessian/dlf_thermal
!!
!! FUNCTION
!!
!! * diagonalise the (mass-weighted) hessian
!! * calculate vibration frequencies
!! * calculate thermal and entropic contributions
!!
!! SYNOPSIS
subroutine dlf_thermal()
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stdout,printl,pi
  use dlf_stat, only: stat
  use dlf_hessian
  implicit none
  real(rk)        :: frequency_factor,svar
  real(rk)        :: temperature,planck,boltz,wavenumber
  real(rk)        :: speed_of_light, angstrom,hartree,avog
  real(rk)        :: evib, frequ, svib, vibtemp, zpe 
  real(rk)        :: zpetot,evibtot,svibtot
  integer         :: ival
! **********************************************************************
  if (.not.glob%massweight .or. glob%icoord /= 0) then
     call dlf_fail("Thermal analysis only possible in mass-weighted &
         &internal coordinates")
  endif
  if (.not.glob%havehessian) then
    call dlf_fail("No hessian present for thermal analysis")
  end if
  ! only act on first node, but print analysis even for printl=0
  if(printl<0) return


!!$  ! tmp print JK
!!$  if(glob%havehessian) then
!!$    write(*,"('HESS',6F10.4)") glob%ihessian
!!$    call dlf_matrix_diagonalise(glob%nivar,glob%ihessian,eigval,eigvec)
!!$    write(stdout,"('Hessian eigenvalues: ',12f10.4)") eigval
!!$    
!!$  end if


  ! constants from NIST (2008)
  speed_of_light=299792458.D0
  planck=6.62606896D-34
  boltz=1.3806504D-23
  avog=6.02214179D23

  ! derived values
  angstrom=0.5291772108D-10
  hartree=4.35974417D-18

  ! sqrt(H/u)/a_B/2/pi/c / 100
  frequency_factor=dsqrt( 4.35974417D-18/ 1.66053886D-27 ) / ( 2.D0 * pi * &
      angstrom * speed_of_light) / 100.D0

  temperature=glob%temperature

  write(stdout,*)
  write(stdout,"('Thermochemical analysis')")
  write(stdout,"('Temperature: ',f10.2,' Kelvin')") temperature
  write(stdout,"(' Mode     Eigenvalue Frequency Vib.T.(K)      &
      &  ZPE (H)   Vib. Ene.(H)      - T*S (H)')")
  zpetot=0.D0
  evibtot=0.D0
  svibtot=0.D0

  ! print out frequencies:
  do ival=1,glob%nivar
    wavenumber = sqrt(abs(eigval(ival))) * frequency_factor

    if(eigval(ival)<0.D0) then
      wavenumber=-wavenumber
      write(stdout,"(i5,f15.10,f10.3)") ival,eigval(ival),wavenumber
    else

      ! convert eig from cm-1 to Hz 
      frequ = wavenumber * 1.D2 * speed_of_light
      ! vibrational temperature in K
      vibtemp = frequ*planck/boltz
      ! zero point vib. (J)
      zpe = planck*frequ*0.5D0
      ! vibrational energy
      svar = dexp(vibtemp/temperature)-1.D0
      evib = planck*frequ / svar
      ! vibrational entropy
      svar = dexp(-vibtemp/temperature)
      svib = vibtemp/temperature / (1.D0/svar -1.D0)
      svib = svib - dlog(1.D0 - svar)
      svib = svib * boltz
      
      write(stdout,"(i5,f15.10,2f10.3,3f15.10)") ival,eigval(ival),wavenumber, &
          vibtemp,zpe/hartree,evib/hartree, -svib/hartree*temperature
      
      zpetot=zpetot+zpe
      evibtot=evibtot+evib
      svibtot=svibtot+svib
    end if
  end do
  write(*,"('total',35x,3f15.10)") &
      zpetot/hartree,evibtot/hartree,-svibtot/hartree*temperature
  write(*,"('total vibrational energy correction to E_electronic',&
      &f15.10,' H')") (zpetot+evibtot-svibtot*temperature)/hartree
  write(*,"('total ZPE  ',f15.5,' J/mol')") zpetot*avog
  write(*,"('total E vib',f15.5,' J/mol')") evibtot*avog
  write(*,"('total S vib',f15.5,' J/mol/K')") svibtot*avog
  
end subroutine dlf_thermal
!!****


