! **********************************************************************
! **                  Utility functions to DL-FIND                    **
! **                                                                  **
! **  Several sub-files:                                              **
! **   dlf_allocate                                                   **
! **   dlf_checkpoint                                                 **
! **   dlf_linalg                                                     **
! **   dlf_time                                                       **
! **                                                                  **
! **********************************************************************
!!****h* DL-FIND/utilities
!!
!! NAME
!! utilities
!!
!! FUNCTION
!! Utility functions to DL-FIND
!!
!! DATA
!! $Date: 2010-05-07 19:30:43 $
!! $Rev: 325 $
!! $Author: gberan $
!! $URL: http://ccpforge.cse.rl.ac.uk/svn/dl-find/branches/release_chemsh3.3/dlf_util.f90 $
!! $Id: dlf_util.f90,v 1.1 2010-05-07 19:30:43 gberan Exp $
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

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
character(2) function get_atom_symbol(atomic_number)
  implicit none
  integer, intent(in) :: atomic_number
  character(2), parameter :: elements(111) = &
       (/ 'H ','He', &
          'Li','Be','B ','C ','N ','O ','F ','Ne', &
          'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
          'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu', &
          'Zn','Ga','Ge','As','Se','Br','Kr', &
          'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag', &
          'Cd','In','Sn','Sb','Te','I ','Xe', &
          'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy', &
          'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt', &
          'Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
          'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf', &
          'Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
          'Rg' /)
! **********************************************************************
  if (atomic_number >= 1 .and. atomic_number <= size(elements)) then
     get_atom_symbol = elements(atomic_number)
  else
     get_atom_symbol = 'XX'
  endif
end function get_atom_symbol 

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine write_xyz(unit,nat,znuc,coords)
  use dlf_parameter_module, only: rk
  implicit none
  integer,intent(in) :: unit
  integer,intent(in) :: nat
  integer,intent(in) :: znuc(nat)
  real(rk),intent(in):: coords(3,nat)
  integer            :: iat
  character(2)       :: str2
  character(2), external :: get_atom_symbol
! **********************************************************************
  write(unit,*) nat
  write(unit,*)
  do iat=1,nat
    str2 = get_atom_symbol(znuc(iat))
    write(unit,'(a2,3f12.7)') str2,coords(:,iat)*0.5291772d0
  end do
  call flush(unit)
end subroutine write_xyz

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine write_xyz_active(unit,nat,znuc,spec,coords)
  use dlf_parameter_module, only: rk
  implicit none
  integer,intent(in) :: unit
  integer,intent(in) :: nat
  integer,intent(in) :: znuc(nat)
  integer,intent(in) :: spec(nat)
  real(rk),intent(in):: coords(3,nat)
  integer            :: iat,nact
  character(2)       :: str2
  character(2), external :: get_atom_symbol
! **********************************************************************
  nact=0
  do iat=1,nat
    if(spec(iat)/=-1) nact=nact+1
  end do
  write(unit,*) nact
  write(unit,*)
  do iat=1,nat
    if(spec(iat)==-1) cycle
    str2 = get_atom_symbol(znuc(iat))
    write(unit,'(a2,3f12.7)') str2,coords(:,iat)*0.5291772d0
  end do
  call flush(unit)
end subroutine write_xyz_active

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module bspline

  ! creates a set of cubic spline functions with the same number of 
  ! grid points (length).
  ! Their x-values may be different
  !
  ! Calling-Sequence:
  !   call spline_init
  !   call spline_create once for each function (each ifunc)
  !   call spline_get for each interpolated value
  !   call spline_destroy

  use dlf_parameter_module, only: rk
  use dlf_allocate, only: allocate,deallocate
  implicit none

  public :: spline_init, spline_create, spline_get, spline_destroy

  interface spline_init
    module procedure spline_init
  end interface

  interface spline_create
    module procedure spline_create
  end interface

  interface spline_get
    module procedure spline_get
  end interface

  interface spline_destroy
    module procedure spline_destroy
  end interface

  private
  logical :: tinit
  integer :: nfunc ! number of functions to be interpolated
  integer :: length ! number of values to be interpolated between
  logical, allocatable, save :: created(:) ! (nfunc)
  real(rk),allocatable, save :: gridx(:,:) ! (length,nfunc) !x values of grid
  real(rk),allocatable, save :: gridy(:,:) ! (length,nfunc) !y values of grid
  real(rk),allocatable, save :: grid_d2ydx2(:,:) ! (length,nfunc)
contains

  subroutine spline_init(length_in,nfunc_in)
    integer, intent(in) :: length_in 
    integer, intent(in) :: nfunc_in 
    ! ******************************************************************

    ! re-initialise if necessary
    if(allocated(created)) call spline_destroy

    nfunc=nfunc_in
    length=length_in
    if(nfunc*length<=0) call dlf_fail("wrong parameters to spline_init")

    tinit=.true.
    call allocate(created,nfunc)
    call allocate(gridx,length,nfunc)
    call allocate(gridy,length,nfunc)
    call allocate(grid_d2ydx2,length,nfunc)

    created(:)=.false.

  end subroutine spline_init

  subroutine spline_destroy
    ! ******************************************************************
    tinit=.false.
    call deallocate(created)
    call deallocate(gridx)
    call deallocate(gridy)
    call deallocate(grid_d2ydx2)

  end subroutine spline_destroy

  SUBROUTINE spline_create(ifunc,x_in,y_in)
    ! creates the spline function number ifunc, i.e. calculates its 
    ! second derivatives at the grid points
    ! x_in must be monotonically increasing
    ! using natural boundary conditions
    ! ******************************************************************
    integer , intent(in) :: ifunc
    real(rk), intent(in) :: x_in(length),y_in(length)
    !
    real(rk) :: store(length)! u
    real(rk) :: svar, svar2  ! sig, p
    integer  :: ival
    ! ******************************************************************

    ! check integrity
    if(.not.tinit) call dlf_fail("spline_create must not be called &
        &before spline_init!")
    if(ifunc<1) call dlf_fail("ifunc < 1 in spline_create")
    if(ifunc>nfunc) call dlf_fail("ifunc > nfunc in spline_create")
    
    gridx(:,ifunc)=x_in(:)
    gridy(:,ifunc)=y_in(:)

    ! natural boundaries
    grid_d2ydx2(1,ifunc)=0.D0
    store(1)=0.D0

    do ival=2,length-1
      svar=(gridx(ival,ifunc)-gridx(ival-1,ifunc))/(gridx(ival+1,ifunc)-gridx(ival-1,ifunc))
      svar2=svar*grid_d2ydx2(ival-1,ifunc)+2.D0
      grid_d2ydx2(ival,ifunc)=(svar-1.D0)/svar2
      store(ival)=(6.D0*((gridy(ival+1,ifunc)-gridy(ival,ifunc))/(gridx(ival+1,ifunc)- &
          gridx(ival,ifunc))-(gridy(ival,ifunc)-gridy(ival-1,ifunc)) /(gridx(ival,ifunc)- &
          gridx(ival-1,ifunc)))/(gridx(ival+1,ifunc)-gridx(ival-1,ifunc))- &
          svar*store(ival-1))/svar2
    enddo

    ! natural boundaries
    grid_d2ydx2(length,ifunc)=0.D0

    do ival=length-1,1,-1 
      grid_d2ydx2(ival,ifunc)=grid_d2ydx2(ival,ifunc)*grid_d2ydx2(ival+1,ifunc)+store(ival) 
    enddo

    created(ifunc)=.true.

  END SUBROUTINE spline_create

  SUBROUTINE spline_get(ifunc,xval,yval,dyval)
    ! calculates a cubic-spline interpolated value (yval) and its 
    ! derivative (dyval) at a position xval
    ! this is done for the interpolation ifunc
    ! ******************************************************************
    integer ,intent(in) :: ifunc
    real(rk),intent(in) :: xval
    real(rk),intent(out):: yval
    real(rk),intent(out):: dyval
    ! 
    integer  :: low,high,ivar
    real(rk) :: aval,bval,delta
    ! ******************************************************************

    ! check integrity
    if(.not.tinit) call dlf_fail("spline_get must not be called before &
        &spline_init!")
    if(ifunc<1) call dlf_fail("ifunc < 1 in spline_get")
    if(ifunc>nfunc) call dlf_fail("ifunc > nfunc in spline_get")
    if(.not.created(ifunc)) call dlf_fail("spline_get must not be called&
        & before spline_create!")
    
    low=1 
    ! bisection on the grid
    high=length
    do while (high-low>1) 
      ivar=(high+low)/2
      if(gridx(ivar,ifunc).gt.xval)then
        high=ivar
      else
        low=ivar
      endif
    end do
    !low and high now bracket the input value of xval.

    delta=gridx(high,ifunc)-gridx(low,ifunc)

    if (delta<=0.D0) call dlf_fail("grid points for spline not distinct")

    aval=(gridx(high,ifunc)-xval)/delta !Cubic spline polynomial is now evaluated.
    bval=(xval-gridx(low,ifunc))/delta
    yval=aval*gridy(low,ifunc)+bval*gridy(high,ifunc) + ((aval**3-aval)* &
        grid_d2ydx2(low,ifunc)+(bval**3-bval)*grid_d2ydx2(high,ifunc)) &
        *(delta**2)/6.D0
    dyval=1.D0/delta * ( gridy(high,ifunc) - gridy(low,ifunc) + delta**2/6.D0 * ( &
        (3.D0*bval**2-1.D0) * grid_d2ydx2(high,ifunc) - &
        (3.D0*aval**2-1.D0) * grid_d2ydx2(low,ifunc) ) ) 
  END SUBROUTINE spline_get

end module bspline
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module dlf_store

  ! Allocates bunches of memory, to store (mainly real number) data
  ! The data are organised in a linked list. They can be set, read out
  ! and deleted
  !
  use dlf_parameter_module, only: rk
  use dlf_global, only: stderr,stdout
  use dlf_allocate, only: allocate,deallocate
  implicit none

  public :: store_initialise, store_allocate, store_set, store_get, &
      store_delete, store_delete_all

  interface store_initialise
    module procedure store_initialise
  end interface

  interface store_allocate
    module procedure store_allocate
  end interface

  interface store_set
    module procedure store_set
    module procedure store_set_a2
    module procedure store_set_a3
  end interface

  interface store_get
    module procedure store_get
    module procedure store_get_a2
    module procedure store_get_a3
  end interface

  interface store_delete
    module procedure store_delete
  end interface

  interface store_delete_all
    module procedure store_delete_all
  end interface

  private

  type store_type_R
    character(40)              :: tag
    integer                    :: size
    real(rk),pointer           :: array(:)
    type(store_type_R),pointer :: next
  end type store_type_R

  type(store_type_R),pointer,save :: first_R
  logical                   ,save :: tinit=.false.

contains

  subroutine store_initialise
    ! ******************************************************************
    if(tinit) call dlf_fail("store is already initialised")
    allocate(first_R)
    nullify(first_R%next)
    nullify(first_R%array)
    first_R%size=0
    first_R%tag=""
    tinit=.true.
  end subroutine store_initialise

  subroutine store_allocate(tag,size)
    character(*)      ,intent(in) :: tag
    integer           ,intent(in) :: size
    type(store_type_R),pointer    :: this
    ! ******************************************************************
    this => first_R
    !print*,"This%size 1",This%size
    ! find last used entry
    do while (associated(this%next))
      !print*," allocate search: this%tag ",this%tag
      if(this%tag == tag) call dlf_fail("Store tag aleady allocated")
      this => this%next
    end do
    ! this is now last set entry
    !print*,"This%size 2",This%size

    ! if this is the first entry, and it is empty, use it, otherwhise use the next
    if(this%size > 0) then
      allocate(this%next)
      this => this%next
      nullify(this%next)
      nullify(this%array)
    end if

    ! set tag and size
    this%tag=tag
    this%size=size
    !print*,"allocated(this%array)",associated(this%array)
    if(associated(this%array)) print*,"shape(this%array)",shape(this%array)
    allocate (this%array(size)) ! pointer

  end subroutine store_allocate

  subroutine store_set(tag,size,array)
    character(*)      ,intent(in) :: tag
    integer           ,intent(in) :: size
    real(rk)          ,intent(in) :: array(size)
    type(store_type_R),pointer    :: this
    ! ******************************************************************

    ! find correct entry
    this => first_R
    do while (associated(this))
      !print*," set search: this%tag ",this%tag
      if(this%tag == tag) exit
      this => this%next
    end do

    if(.not.associated(this)) then
      write(stdout,*) "Storage tag ",tag," not found!"
      call dlf_fail("Storage tag to set not found")
    end if

    if(this%size /= size) call dlf_fail("Storage set size inconsistent")

    this%array(:)=array(:)

    !print*,"Tag ",tag," now set"

  end subroutine store_set

  subroutine store_get(tag,size,array)
    character(*)      ,intent(in) :: tag
    integer           ,intent(in) :: size
    real(rk)          ,intent(out):: array(size)
    type(store_type_R),pointer    :: this
    ! ******************************************************************

    ! find correct entry
    this => first_R
    do while (associated(this))
      !print*," get search: this%tag ",this%tag
      if(this%tag == tag) exit
      this => this%next
    end do

    if(.not.associated(this)) then
      write(stdout,*) "Storage tag ",tag," not found!"
      call dlf_fail("Storage tag to get not found")
    end if

    if(this%size /= size) call dlf_fail("Storage get size inconsistent")

    array(:)=this%array(:)
    !print*,"Tag ",tag," now got"

  end subroutine store_get

  subroutine store_delete(tag)
    ! delete one tag
    character(*)      ,intent(in) :: tag
    type(store_type_R),pointer    :: this,del

    ! find correct entry
    this => first_R

    ! handle the case of deletion of the first entry separately
    if(this%tag == tag) then
      if(.not.associated(this%next)) then
        ! first entry is the only entry
        !print*," deleting first and only entry ",tag
        if(associated(this%array)) deallocate(this%array) ! pointer
        this%size=0
        this%tag=""
      else
        ! other entries exist
        !print*," deleting first entry ",tag
        first_R => this%next
        if(associated(this%array)) deallocate(this%array) ! pointer
        deallocate(this)
      end if
      return
    end if

    do while (associated(this%next))
      !print*," delete search: this%next%tag ",this%next%tag
      if(this%next%tag == tag) exit
      this => this%next
    end do
    ! this%next points to the tag to be deleted

    if(.not.associated(this%next)) then
      write(stdout,*) "Storage tag ",tag," not found!"
      call dlf_fail("Storage tag to delete not found")
    end if

    ! delete this%next
    this%next%size=0
    if(associated(this%next%array)) then
      deallocate(this%next%array) ! pointer
      !print*,"Deallocating tag ",this%next%tag
    else
      !print*,"Warning: not able to deallocate tag",tag
    end if

    if(associated(this%next%next)) then
      del => this%next%next
      deallocate(this%next)
      this%next => del
    else
      ! the one to be deleted is the last one
      deallocate(this%next)
      nullify(this%next)
    end if

  end subroutine store_delete

  subroutine store_delete_all
    ! delete all tags
    type(store_type_R),pointer    :: this,next

    if(.not.tinit) return

    this => first_R

    do while (associated(this%next))
      next => this%next
      !print*," deleteing this%tag ",this%tag

      if(associated(this%array)) deallocate(this%array) ! pointer
      deallocate(this)

      this => next
    end do

    ! now deallocate the last one
    !print*," deleteing this%tag ",this%tag

    if(associated(this%array)) deallocate(this%array) ! pointer
    deallocate(this)

    tinit=.false.

  end subroutine store_delete_all

  subroutine store_set_a2(tag,size,array)
    ! dummy routine to cover rank 2 arrays
    character(*)      ,intent(in) :: tag
    integer           ,intent(in) :: size
    real(rk)          ,intent(in) :: array(:,:)
    call store_set(tag,size,reshape(array,(/size/)))
  end subroutine store_set_a2

  subroutine store_set_a3(tag,size,array)
    ! dummy routine to cover rank 3 arrays
    character(*)      ,intent(in) :: tag
    integer           ,intent(in) :: size
    real(rk)          ,intent(in) :: array(:,:,:)
    call store_set(tag,size,reshape(array,(/size/)))
  end subroutine store_set_a3

  subroutine store_get_a2(tag,size1,size2,array)
    ! dummy routine to cover rank 2 arrays
    character(*)      ,intent(in) :: tag
    integer           ,intent(in) :: size1
    integer           ,intent(in) :: size2
    real(rk)          ,intent(out):: array(:,:)
    real(rk)  :: tmp_array(size1*size2)
    call store_get(tag,size1*size2,tmp_array)
    array=reshape(tmp_array,(/size1,size2/))
  end subroutine store_get_a2

  subroutine store_get_a3(tag,size1,size2,size3,array)
    ! dummy routine to cover rank 3 arrays
    character(*)      ,intent(in) :: tag
    integer           ,intent(in) :: size1
    integer           ,intent(in) :: size2
    integer           ,intent(in) :: size3
    real(rk)          ,intent(out):: array(:,:,:)
    real(rk)  :: tmp_array(size1*size2*size3)
    call store_get(tag,size1*size2*size3,tmp_array)
    array=reshape(tmp_array,(/size1,size2,size3/))
  end subroutine store_get_a3

end module dlf_store
