! **********************************************************************
! **        Driver file for dl_find                                   **
! ** isystem:                                                         **
! **  1 - Mueller-Brown potential (2 Dim)                             **
! **  2 - 3 Lennard-Jones Atoms (9 Dim)                               **
! **  3 - 100 Lennard-Jones Atoms (300 Dim)                           **
! **                                                                  **
! **                                                                  **
! **********************************************************************

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

module driver_module
  integer,parameter :: isystem=3 ! decide which system to run on
end module driver_module

program main
  use driver_module
  implicit none
  integer :: ivar

  call dlf_mpi_initialize() ! only necessary for a parallel build; 
                            ! can be present for a serial build
  call dlf_output(6,0)

  select case (isystem)
  case (1)
    !call system ("rm -f dimer.xy")
    !do ivar=1,100
      call dl_find(2,1,0,1)
    !  call system ("echo '' >> dimer.xy")
    !end do
  case (2)
   ! call dl_find(12,1,4)
    call dl_find(12,16,4,1) ! 1 frame + masses
  case (3)
    !call dl_find(99,1,33)
    !call dl_find(30,1,15)

    ! LJ-particle surface with one atom hopping on it
    call dl_find(21,35,7,1) ! 1 frame + weigths + masses

  case default
    call dlf_mpi_finalize() ! only necessary for a parallel build;
                            ! can be present for a serial build

    stop "Wrong isystem"
  end select

  call dlf_mpi_finalize() ! only necessary for a parallel build;
                          ! can be present for a serial build


end program main

! **********************************************************************
! subroutines that have to be provided to dl_find from outside
! **********************************************************************

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_params(nvar,nvar2,nspec,coords,coords2,spec,ierr, &
    tolerance,printl,maxcycle,maxene,tatoms,icoord, &
    iopt,iline,maxstep,scalestep,lbfgs_mem,nimage,nebk,dump,restart,&
    nz,ncons,nconn,update,maxupd,delta,soft,inithessian,carthessian,tsrel, &
    maxrot,tolrot,nframe,nmass,nweight,timestep,fric0,fricfac,fricp, &
    imultistate, state_i,state_j,pf_c1,pf_c2,gp_c3,gp_c4,ln_t1,ln_t2, &
    printf,tolerance_e,distort,massweight,minstep,maxdump,task,temperature, &
    po_pop_size,po_radius,po_contraction,po_tolerance_r,po_tolerance_g, &
    po_distribution,po_maxcycle,po_init_pop_size,po_reset,po_mutation_rate, &
    po_death_rate,po_scalefac,po_nsave,ntasks,tdlf_farm,n_po_scaling)

  use dlf_parameter_module, only: rk
  use driver_module
  implicit none
  integer   ,intent(in)      :: nvar 
  integer   ,intent(in)      :: nvar2
  integer   ,intent(in)      :: nspec
  real(rk)  ,intent(inout)   :: coords(nvar) ! start coordinates
  real(rk)  ,intent(inout)   :: coords2(nvar2) ! a real array that can be used
                                               ! depending on the calculation
                                               ! e.g. a second set of coordinates
  integer   ,intent(inout)   :: spec(nspec)  ! specifications like fragment or frozen
  integer   ,intent(out)     :: ierr
  real(rk)  ,intent(inout)   :: tolerance
  real(rk)  ,intent(inout)   :: tolerance_e
  integer   ,intent(inout)   :: printl
  integer   ,intent(inout)   :: maxcycle
  integer   ,intent(inout)   :: maxene
  integer   ,intent(inout)   :: tatoms
  integer   ,intent(inout)   :: icoord
  integer   ,intent(inout)   :: iopt
  integer   ,intent(inout)   :: iline
  real(rk)  ,intent(inout)   :: maxstep
  real(rk)  ,intent(inout)   :: scalestep
  integer   ,intent(inout)   :: lbfgs_mem
  integer   ,intent(inout)   :: nimage
  real(rk)  ,intent(inout)   :: nebk
  integer   ,intent(inout)   :: dump
  integer   ,intent(inout)   :: restart
  integer   ,intent(inout)   :: nz
  integer   ,intent(inout)   :: ncons
  integer   ,intent(inout)   :: nconn
  integer   ,intent(inout)   :: update
  integer   ,intent(inout)   :: maxupd
  real(rk)  ,intent(inout)   :: delta
  real(rk)  ,intent(inout)   :: soft
  integer   ,intent(inout)   :: inithessian
  integer   ,intent(inout)   :: carthessian
  integer   ,intent(inout)   :: tsrel
  integer   ,intent(inout)   :: maxrot
  real(rk)  ,intent(inout)   :: tolrot
  integer   ,intent(inout)   :: nframe
  integer   ,intent(inout)   :: nmass
  integer   ,intent(inout)   :: nweight
  real(rk)  ,intent(inout)   :: timestep
  real(rk)  ,intent(inout)   :: fric0
  real(rk)  ,intent(inout)   :: fricfac
  real(rk)  ,intent(inout)   :: fricp
  integer   ,intent(inout)   :: imultistate
  integer   ,intent(inout)   :: state_i
  integer   ,intent(inout)   :: state_j
  real(rk)  ,intent(inout)   :: pf_c1  
  real(rk)  ,intent(inout)   :: pf_c2  
  real(rk)  ,intent(inout)   :: gp_c3  
  real(rk)  ,intent(inout)   :: gp_c4
  real(rk)  ,intent(inout)   :: ln_t1  
  real(rk)  ,intent(inout)   :: ln_t2  
  integer   ,intent(inout)   :: printf
  real(rk)  ,intent(inout)   :: distort
  integer   ,intent(inout)   :: massweight
  real(rk)  ,intent(inout)   :: minstep
  integer   ,intent(inout)   :: maxdump
  integer   ,intent(inout)   :: task
  real(rk)  ,intent(inout)   :: temperature
  integer   ,intent(inout)   :: po_pop_size
  real(rk)  ,intent(inout)   :: po_radius
  real(rk)  ,intent(inout)   :: po_contraction
  real(rk)  ,intent(inout)   :: po_tolerance_r
  real(rk)  ,intent(inout)   :: po_tolerance_g
  integer   ,intent(inout)   :: po_distribution
  integer   ,intent(inout)   :: po_maxcycle
  integer   ,intent(inout)   :: po_init_pop_size
  integer   ,intent(inout)   :: po_reset
  real(rk)  ,intent(inout)   :: po_mutation_rate
  real(rk)  ,intent(inout)   :: po_death_rate
  real(rk)  ,intent(inout)   :: po_scalefac
  integer   ,intent(inout)   :: po_nsave
  integer   ,intent(inout)   :: ntasks
  integer   ,intent(inout)   :: tdlf_farm
  integer   ,intent(inout)   :: n_po_scaling

  ! local variables
  real(rk)                   :: svar
  integer                    :: iat,jat
  interface
    subroutine read_rand(arr)
      use dlf_parameter_module, only: rk
      use dlf_global, only : glob
      real(rk)  :: arr(:)
    end subroutine read_rand
  end interface
! **********************************************************************
  ierr=0
  tsrel=1
  
  print*,"External sizes:",nvar,nvar2,nspec
  ! Minima of Mueller-Brown potential:
  ! left:   -0.5582236346340204   1.4417258418038705   energy:  -0.6919788547639341
  ! middle: -0.050010822944531706 0.4666941048659066   energy:  -0.38098027419650493
  ! right:   0.623499404927291    0.028037758526434815 energy:  -0.5102203967776056
  ! TS:
  ! main:   -0.8220015634663578   0.624312796202859    energy:  -0.19181529956913862
  ! other:   0.21248658127591744  0.2929883285085865   energy:  -0.34079688732228874

!!$  spec(1)=-1
!!$  spec(2)=-34
!!$  spec(3)=-4
!!$  spec(:)=1
!!$  spec(3)=-34
!!$  spec(1:3)=2
!!$  spec(6)=2
!!$  spec(9:10)=0
  ncons=1
  ! cosntraint
!!$  spec(11)=1
!!$  spec(12)=1
!!$  spec(13)=2
!!$  spec(14:15)=0

  ncons=0
  nconn=0

  coords(:)=-1.D0

  
!!$  ! Near Saddle point:
! coords=(/ -0.6D0, 0.6D0 /) 

  coords2(:)=-1.D0
  if(isystem==1) then
!!$    coords(:)=(/-0.5582236D0,1.4417258D0/)
!!$    !coords(:)=(/-0.05001082294D0, 0.4666941D0/)
!!$    coords2(:)=(/0.6234994D0,  0.0280377585D0/)
!!$
    ! somewhere near minima
    !coords(:)=(/-0.3D0,0.6D0/)
    !coords2(:)=(/-0.31D0,0.590D0/)
    call random_number(coords)
    coords(1)=coords(1)-0.7D0
    coords(2)=coords(2)*1.4D0
    !call random_number(coords2)
    !coords2=coords+(/-0.733D0, 0.680D0/)! coords2-0.5D0
    !coords(:)=(/-0.05001082294D0, 0.4666941D0/)
    coords(:)=(/0.6D0,  0.1D0/)
  end if
  if(isystem==2) then
    ! minimum distance is 1.1224620
!!$    coords(:)=0.D0
!!$    coords(4)=-1.D0
!!$    coords(7)=0.6D0
!!$    coords(8)=1.D0
!!$
!!$    coords2(:)=0.D0
!!$    coords2(4)=-1.D0
!!$    coords2(7)=0.6D0
!!$    coords2(8)=-1.D0

    ! Arrangements of 4 LJ-Particles:
    ! Square with a=1.1126198 is a 2nd order Saddle point
    ! Rombus (60.27 deg) with a=1.1202310 (1.124800) is a 1st order Saddle point

    ! This is a rombus quite near to the TS
    svar=1.1224620D0 ! very good
    svar=1.0D0 
    coords(:)=0.D0
    coords(4:6)  =(/svar,0.D0,0.D0/)
    coords(7:9)  =(/svar*0.5D0, svar*0.5D0*dsqrt(3.D0), 0.D0/)
    coords(10:12)=(/svar*1.5D0, svar*0.5D0*dsqrt(3.D0), 0.D0/)

    coords(12)=1.D0

    ! Two points quite near to a TS:
    svar=1.1224620D0*0.9D0
!!$    coords(:)=0.D0
!!$    coords(4:6)  =(/svar,0.D0,0.D0/)
!!$    coords(7:9)  =(/svar*0.5D0, svar*0.5D0*dsqrt(3.D0), 0.D0/)
!!$    coords(10:12)=(/svar*1.5D0, svar*0.5D0*dsqrt(3.D0), 0.1D0/)
    coords2(:)=0.D0
    coords2(4:6)  =(/svar,0.D0,0.D0/)
    coords2(7:9)  =(/svar*0.5D0, svar*0.5D0*dsqrt(3.D0), 0.D0/)
    coords2(10:12)=(/svar*1.5D0, svar*0.5D0*dsqrt(3.D0), -0.3D0/)
    nframe=1

!!$    coords2(:)=0.D0
!!$    coords2(1:2)=(/-svar,svar*sqrt(3.D0)/)
!!$    coords2(7:8)=(/svar,svar*sqrt(3.D0)/)

!!$    ! four atoms in nearly a square
!!$    coords(:)=0.D0
!!$    coords(3)=0.d0
!!$    coords(4)=1.0D0
!!$    coords(8)=1.0D0
!!$    coords(10)=1.0D0
!!$    coords(11)=1.0D0
!!$    coords(12)=0.610D0

    ! masses
    nmass=4
    coords2(13:14)=1.D0
    coords2(15:16)=10.D0

  end if
  if(isystem==3) then
    ! one hopping atom on a 2*3 surface
    svar=1.1224620D0
    do iat=0,1
      do jat=0,2
        coords(3+iat*9+jat*3+1)=dble(jat)*svar
        coords(3+iat*9+jat*3+2)=dble(iat)*svar
        coords(3+iat*9+jat*3+3)=0.D0
      end do
    end do
    ! TS energy with 4 frozen atoms: -0.01156076358016981
    ! hopping atom (atom 1) 
    coords(1)=svar*0.7D0
    coords(2)=svar*0.5D0
    coords(3)=svar*0.7D0
    !coords(9)=-0.3D0
    !coords(18)=-0.3D0
    nmass=7
    nframe=1
    coords2(:)=1.D0
    coords2(1:nvar)=coords(:)
    !coords2(nvar+1:2*nvar)=coords(:)
    ! hop atom other minimum
    coords2(1)=coords2(1)+svar
    !coords2(9)=coords2(9)+0.1D0
    !coords2(3)=coords2(3)+0.2D0
    ! for dimer:
    !coords2(1)=coords2(1)+0.6D0*svar
    !coords2(3)=coords2(3)+0.2D0
    !coords2(nvar+1)=coords2(nvar+1)+svar
    spec(:)=-1 ! NEB: -1
    spec(1)=1
    spec(3)=1
    spec(6)=1

    ! weights
    nweight=7
    !coords2(22:28) are weights
!    coords2(22:28)=1.0d0
!    coords2(24)=1.d0

    ! masses
    !coords2(29:35) are masses
    coords2(29:35)=10.D0
!    coords2(29:29)=1.D0
  end if

  tolerance=-4.5D-1 ! negative: default settings
  printl=4
  printf=4
  maxcycle=500
  maxene=100000

  !tolrot=5.D-5
  !maxrot=100

  task=0 !1011

  distort=0.0D0 !-0.2D0
  tatoms=0
  icoord=221
  massweight=0
  iopt=11
  iline=0
  maxstep=0.3D0
  scalestep=1.0D0
  lbfgs_mem=40
  nimage=7
  nebk=1.D-2

  ! Hessian
  delta=1.D-4
!  soft=-6.D-4
  update=0
  maxupd=50
  inithessian = 0
  minstep=1.D-5

  ! damped dynamics
  fric0=0.1D0
  fricfac=1.0D0
  fricp=0.1D0

  dump=0
  restart=0
  temperature = 300.0D0 ! K

  ! Parallel optimization

  po_pop_size=25
  po_radius=0.5D0
  po_init_pop_size=50
  po_contraction=0.95D0
  po_tolerance_r=1.0D-8
  po_tolerance_g=1.0D-6
  po_distribution=3
  po_maxcycle=100000
  po_reset=500
  po_mutation_rate=0.15D0
  po_death_rate=0.5D0
  po_scalefac=10.0D0
  po_nsave=10
  n_po_scaling=0 ! meaning that the base radii values for the sampling and tolerance 
                 ! are used for all components of the working coordinate vector.
                 ! Remember to change the second arg in the call to dl_find if a 
                 ! non-zero value for n_po_scaling is desired, and also to add the 
                 ! necessary values to the coords2 array...

  ! Taskfarming 
  ! (the two lines below, with any values assigned, 
  ! may be safely left in place for a serial build) 
  ntasks = 1
  tdlf_farm = 1

  ! PARAMETERS FOR specific systems
  select case (isystem)
  case (1)
    tatoms=0
    !
  case (2)
    tatoms=1
  !  call read_rand(coords)
  !  coords=coords*2.D0
!    coords(1)=coords(1)+10.D0
  case (3)
    tatoms=1
!!$    lbfgs_mem=40
!!$    call read_rand(coords)
!!$    coords=coords*2.D0
  end select

end subroutine dlf_get_params

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_gradient(nvar,coords,energy,gradient,iimage,status)
  !  Mueller-Brown Potential
  !  see K Mueller and L. D. Brown, Theor. Chem. Acta 53, 75 (1979)
  !  taken from JCP 111, 9475 (1999)
  use dlf_parameter_module, only: rk
  use driver_module
  implicit none
  integer   ,intent(in)    :: nvar
  real(rk)  ,intent(in)    :: coords(nvar)
  real(rk)  ,intent(out)   :: energy
  real(rk)  ,intent(out)   :: gradient(nvar)
  integer   ,intent(in)    :: iimage
  integer   ,intent(out)   :: status
  !
  ! variables for Mueller-Brown potential
  real(rk) :: acappar(4),apar(4),bpar(4),cpar(4),x0par(4),y0par(4)
  real(rk) :: ebarr
  real(rk) :: x,y,svar,svar2
  integer  :: icount
  ! variables for Lennard-Jones potentials
  real(rk) :: acoords(3,nvar/3)
  real(rk) :: agrad(3,nvar/3)
  real(rk) :: epsilon=1.D-1
  real(rk) :: sigma=1.D0
  real(rk) :: r
  integer  :: nat,iat,jat
! **********************************************************************
  status=1
  select case (isystem)
  case (1)
    print*,"coords in energy eval",coords
    ! assign parameters (could be moved to something external ...)
    ebarr=0.5D0
    acappar(1)=-200.D0*ebarr/106.D0
    acappar(2)=-100.D0*ebarr/106.D0
    acappar(3)=-170.D0*ebarr/106.D0
    acappar(4)=  15.D0*ebarr/106.D0
    apar(1)=-1.D0
    apar(2)=-1.D0
    apar(3)=-6.5D0
    apar(4)=0.7D0
    bpar(1)=0.D0
    bpar(2)=0.D0
    bpar(3)=11.D0
    bpar(4)=0.6D0
    cpar(1)=-10.D0
    cpar(2)=-10.D0
    cpar(3)=-6.5D0
    cpar(4)=0.7D0
    x0par(1)=1.D0
    x0par(2)=0.D0
    x0par(3)=-0.5D0
    x0par(4)=-1.D0
    y0par(1)=0.D0
    y0par(2)=0.5D0
    y0par(3)=1.5D0
    y0par(4)=1.D0

    x =  coords(1)
    y =  coords(2)

    energy=0.D0
    gradient=0.D0
    do icount=1,4
      svar= apar(icount)*(x-x0par(icount))**2 + &
          bpar(icount)*(x-x0par(icount))*(y-y0par(icount)) + &
          cpar(icount)*(y-y0par(icount))**2 
      svar2= acappar(icount) * dexp(svar)
      energy=energy+ svar2
      gradient(1)=gradient(1) + svar2 * &
          (2.D0* apar(icount)*(x-x0par(icount))+bpar(icount)*(y-y0par(icount)))
      gradient(2)=gradient(2) + svar2 * &
          (2.D0* cpar(icount)*(y-y0par(icount))+bpar(icount)*(x-x0par(icount)))
    end do
    !  write(*,'("x,y,func",2f10.5,es15.7)') x,y,energy
  case (2,3)

    ! one could use a Lennard-Jones particle with slightly different
    ! parameters as "Excited state": epsilon=0.9D-3, sigma=1.1D0 to search for
    ! conical intersections

    acoords=reshape(coords,(/3,nvar/3/))
    energy=0.D0
    agrad(:,:)=0.D0
    do iat=1,nvar/3
      do jat=iat+1,nvar/3
        r=sum((acoords(:,iat)-acoords(:,jat))**2)
        ! Lennard-Jones Potential
        energy=energy+4.D0*epsilon * ((sigma**2/r)**6-(sigma**2/r)**3)
        svar = -4.D0*epsilon * (12.D0*sigma**12/r**7-6.D0*sigma**6/r**4)
        agrad(1,iat)=agrad(1,iat)+ svar * (acoords(1,iat)-acoords(1,jat))
        agrad(2,iat)=agrad(2,iat)+ svar * (acoords(2,iat)-acoords(2,jat))
        agrad(3,iat)=agrad(3,iat)+ svar * (acoords(3,iat)-acoords(3,jat))
        agrad(1,jat)=agrad(1,jat)- svar * (acoords(1,iat)-acoords(1,jat))
        agrad(2,jat)=agrad(2,jat)- svar * (acoords(2,iat)-acoords(2,jat))
        agrad(3,jat)=agrad(3,jat)- svar * (acoords(3,iat)-acoords(3,jat))
      end do
    end do
    gradient=reshape(agrad,(/nvar/))
  end select
  status=0
end subroutine dlf_get_gradient

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_hessian(nvar,coords,hessian,status)
  !  get the hessian at a given geometry
  use dlf_parameter_module
  use driver_module
  implicit none
  integer   ,intent(in)    :: nvar
  real(rk)  ,intent(in)    :: coords(nvar)
  real(rk)  ,intent(out)   :: hessian(nvar,nvar)
  integer   ,intent(out)   :: status
  real(rk) :: epsilon=1.D-1 ! must be the same as in get_gradient
  real(rk) :: sigma=1.D0
  real(rk) :: acoords(3,nvar/3),r,svar,svar2
  integer  :: posi,posj,iat,jat,m,n
  ! variables for Mueller-Brown potential
  real(rk) :: acappar(4),apar(4),bpar(4),cpar(4),x0par(4),y0par(4)
  real(rk) :: ebarr
  real(rk) :: x,y
  integer  :: icount
! **********************************************************************
  hessian(:,:)=0.D0
  status=1
  select case (isystem)
  case (1)
    ebarr=0.5D0
    acappar(1)=-200.D0*ebarr/106.D0
    acappar(2)=-100.D0*ebarr/106.D0
    acappar(3)=-170.D0*ebarr/106.D0
    acappar(4)=  15.D0*ebarr/106.D0
    apar(1)=-1.D0
    apar(2)=-1.D0
    apar(3)=-6.5D0
    apar(4)=0.7D0
    bpar(1)=0.D0
    bpar(2)=0.D0
    bpar(3)=11.D0
    bpar(4)=0.6D0
    cpar(1)=-10.D0
    cpar(2)=-10.D0
    cpar(3)=-6.5D0
    cpar(4)=0.7D0
    x0par(1)=1.D0
    x0par(2)=0.D0
    x0par(3)=-0.5D0
    x0par(4)=-1.D0
    y0par(1)=0.D0
    y0par(2)=0.5D0
    y0par(3)=1.5D0
    y0par(4)=1.D0

    x =  coords(1)
    y =  coords(2)

    hessian=0.D0
    do icount=1,4
      svar= apar(icount)*(x-x0par(icount))**2 + &
          bpar(icount)*(x-x0par(icount))*(y-y0par(icount)) + &
          cpar(icount)*(y-y0par(icount))**2 
      svar2= acappar(icount) * dexp(svar)
      !energy=energy+ svar2
!!$      gradient(1)=gradient(1) + svar2 * &
!!$          (2.D0* apar(icount)*(x-x0par(icount))+bpar(icount)*(y-y0par(icount)))
!!$      gradient(2)=gradient(2) + svar2 * &
!!$          (2.D0* cpar(icount)*(y-y0par(icount))+bpar(icount)*(x-x0par(icount)))
      hessian(1,1)=hessian(1,1)+svar2 * &
          (2.D0* apar(icount)*(x-x0par(icount))+bpar(icount)*(y-y0par(icount)))**2 + &
          svar2 * 2.D0 * apar(icount)
      hessian(2,2)=hessian(2,2)+svar2 * &
          (2.D0* cpar(icount)*(y-y0par(icount))+bpar(icount)*(x-x0par(icount)))**2 + &
          svar2 * 2.D0 * cpar(icount)
      hessian(1,2)=hessian(1,2) + svar2 * &
          (2.D0* apar(icount)*(x-x0par(icount))+bpar(icount)*(y-y0par(icount)))* &
          (2.D0* cpar(icount)*(y-y0par(icount))+bpar(icount)*(x-x0par(icount))) + &
          svar2 * bpar(icount)
    end do
    hessian(2,1)=hessian(1,2)
    status=0

  case (2,3)
    acoords=reshape(coords,(/3,nvar/3/))
    do iat=1,nvar/3
      do jat=iat+1,nvar/3
        r=sum((acoords(:,iat)-acoords(:,jat))**2)
        ! Lennard-Jones Potential
        svar = 96.D0*epsilon * (7.D0*sigma**12/r**8-2.D0*sigma**6/r**5) ! coeff of x1x2
        svar2= epsilon * (-2.D0*sigma**12/r**7+sigma**6/r**4) ! for x1x1
        posi=(iat-1)*3+1
        posj=(jat-1)*3+1
        ! off-diag
        hessian(posi,posi+1)  =hessian(posi,posi+1)  + svar * (acoords(1,iat)-acoords(1,jat)) * (acoords(2,iat)-acoords(2,jat))
        hessian(posi,posi+2)  =hessian(posi,posi+2)  + svar * (acoords(1,iat)-acoords(1,jat)) * (acoords(3,iat)-acoords(3,jat))
        hessian(posi+1,posi)  =hessian(posi+1,posi)  + svar * (acoords(2,iat)-acoords(2,jat)) * (acoords(1,iat)-acoords(1,jat))
        hessian(posi+1,posi+2)=hessian(posi+1,posi+2)+ svar * (acoords(2,iat)-acoords(2,jat)) * (acoords(3,iat)-acoords(3,jat))
        hessian(posi+2,posi)  =hessian(posi+2,posi)  + svar * (acoords(3,iat)-acoords(3,jat)) * (acoords(1,iat)-acoords(1,jat))
        hessian(posi+2,posi+1)=hessian(posi+2,posi+1)+ svar * (acoords(3,iat)-acoords(3,jat)) * (acoords(2,iat)-acoords(2,jat))

        do m=0,2
          do n=0,2
            if(m==n) cycle
  hessian(posi+m,posj+n)=hessian(posi+m,posj+n)- svar * (acoords(M+1,iat)-acoords(M+1,jat)) * (acoords(N+1,iat)-acoords(N+1,jat))
  hessian(posj+m,posi+n)=hessian(posj+m,posi+n)- svar * (acoords(M+1,iat)-acoords(M+1,jat)) * (acoords(N+1,iat)-acoords(N+1,jat))
          end do
        end do
        ! Diag for different atoms ...
        do m=0,2
          hessian(posi+m,posj+m)=hessian(posi+m,posj+m) - 24.D0*(svar2+1.D0/24.D0*svar* (acoords(m+1,iat)-acoords(M+1,jat))**2)
          hessian(posj+m,posi+m)=hessian(posj+m,posi+m) - 24.D0*(svar2+1.D0/24.D0*svar* (acoords(m+1,iat)-acoords(M+1,jat))**2)
        end do

        hessian(posj,posj+1)  =hessian(posj,posj+1)  + svar * (acoords(1,iat)-acoords(1,jat)) * (acoords(2,iat)-acoords(2,jat))
        hessian(posj,posj+2)  =hessian(posj,posj+2)  + svar * (acoords(1,iat)-acoords(1,jat)) * (acoords(3,iat)-acoords(3,jat))
        hessian(posj+1,posj)  =hessian(posj+1,posj)  + svar * (acoords(2,iat)-acoords(2,jat)) * (acoords(1,iat)-acoords(1,jat))
        hessian(posj+1,posj+2)=hessian(posj+1,posj+2)+ svar * (acoords(2,iat)-acoords(2,jat)) * (acoords(3,iat)-acoords(3,jat))
        hessian(posj+2,posj)  =hessian(posj+2,posj)  + svar * (acoords(3,iat)-acoords(3,jat)) * (acoords(1,iat)-acoords(1,jat))
        hessian(posj+2,posj+1)=hessian(posj+2,posj+1)+ svar * (acoords(3,iat)-acoords(3,jat)) * (acoords(2,iat)-acoords(2,jat))
        ! diag
        hessian(posi,posi)    =hessian(posi,posi)    + 24.D0*(svar2+1.D0/24.D0*svar* (acoords(1,iat)-acoords(1,jat))**2)
        hessian(posi+1,posi+1)=hessian(posi+1,posi+1)+ 24.D0*(svar2+1.D0/24.D0*svar* (acoords(2,iat)-acoords(2,jat))**2)
        hessian(posi+2,posi+2)=hessian(posi+2,posi+2)+ 24.D0*(svar2+1.D0/24.D0*svar* (acoords(3,iat)-acoords(3,jat))**2)

        hessian(posj,posj)    =hessian(posj,posj)    + 24.D0*(svar2+1.D0/24.D0*svar* (acoords(1,iat)-acoords(1,jat))**2)
        hessian(posj+1,posj+1)=hessian(posj+1,posj+1)+ 24.D0*(svar2+1.D0/24.D0*svar* (acoords(2,iat)-acoords(2,jat))**2)
        hessian(posj+2,posj+2)=hessian(posj+2,posj+2)+ 24.D0*(svar2+1.D0/24.D0*svar* (acoords(3,iat)-acoords(3,jat))**2)
      end do
    end do
    status=0
  end select
end subroutine dlf_get_hessian

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_put_coords(nvar,mode,energy,coords,iam)
  use dlf_parameter_module
  implicit none
  integer   ,intent(in)    :: nvar
  integer   ,intent(in)    :: mode
  integer   ,intent(in)    :: iam
  real(rk)  ,intent(in)    :: energy
  real(rk)  ,intent(in)    :: coords(nvar)
  integer                  :: iat
! **********************************************************************

! Only do this writing of files if I am the rank-zero processor
  if (iam /= 0) return

  if(mod(nvar,3)==0) then
    !assume coords are atoms
    if(mode==2) then
      open(unit=20,file="tsmode.xyz")
    else
      open(unit=20,file="coords.xyz")
    end if
    write(20,*) nvar/3
    write(20,*) 
    do iat=1,nvar/3
      write(20,'("H ",3f12.7)') coords((iat-1)*3+1:(iat-1)*3+3)
    end do
    close(20)
  else
    !print*,"Coords in put_coords: ",coords
  end if
end subroutine dlf_put_coords

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_error()
  implicit none
! **********************************************************************
  call dlf_mpi_abort() ! only necessary for a parallel build;
                       ! can be present for a serial build
  call exit(1)
end subroutine dlf_error

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_update()
  implicit none
! **********************************************************************
  ! only a dummy routine here.
end subroutine dlf_update


subroutine dlf_get_multistate_gradients(nvar,coords,energy,gradient,iimage,status)
  ! only a dummy routine up to now
  ! for conical intersection search
  use dlf_parameter_module
  implicit none
  integer   ,intent(in)    :: nvar
  integer   ,intent(in)    :: coords(nvar)
  real(rk)  ,intent(in)    :: energy(2)
  real(rk)  ,intent(in)    :: gradient(nvar,2)
  integer   ,intent(in)    :: iimage
  integer   ,intent(in)    :: status
end subroutine dlf_get_multistate_gradients


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_put_procinfo(dlf_nprocs, dlf_iam, dlf_global_comm)

  implicit none

  integer, intent(in) :: dlf_nprocs ! total number of processors
  integer, intent(in) :: dlf_iam ! my rank, from 0, in mpi_comm_world
  integer, intent(in) :: dlf_global_comm ! world-wide communicator
! **********************************************************************

!!! variable in the calling program = corresponding dummy argument

end subroutine dlf_put_procinfo


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_procinfo(dlf_nprocs, dlf_iam, dlf_global_comm)

  implicit none

  integer :: dlf_nprocs ! total number of processors
  integer :: dlf_iam ! my rank, from 0, in mpi_comm_world
  integer :: dlf_global_comm ! world-wide communicator
! **********************************************************************

!!! dummy argument = corresponding variable in the calling program

end subroutine dlf_get_procinfo


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_put_taskfarm(dlf_ntasks, dlf_nprocs_per_task, dlf_iam_in_task, &
                        dlf_mytask, dlf_task_comm, dlf_ax_tasks_comm)

  implicit none

  integer, intent(in) :: dlf_ntasks          ! number of taskfarms
  integer, intent(in) :: dlf_nprocs_per_task ! no of procs per farm
  integer, intent(in) :: dlf_iam_in_task     ! my rank, from 0, in my farm
  integer, intent(in) :: dlf_mytask          ! rank of my farm, from 0
  integer, intent(in) :: dlf_task_comm       ! communicator within each farm
  integer, intent(in) :: dlf_ax_tasks_comm   ! communicator involving the 
                                             ! i-th proc from each farm
! **********************************************************************

!!! variable in the calling program = corresponding dummy argument 

end subroutine dlf_put_taskfarm


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_taskfarm(dlf_ntasks, dlf_nprocs_per_task, dlf_iam_in_task, &
                        dlf_mytask, dlf_task_comm, dlf_ax_tasks_comm)

  implicit none

  integer :: dlf_ntasks          ! number of taskfarms
  integer :: dlf_nprocs_per_task ! no of procs per farm
  integer :: dlf_iam_in_task     ! my rank, from 0, in my farm
  integer :: dlf_mytask          ! rank of my farm, from 0
  integer :: dlf_task_comm       ! communicator within each farm
  integer :: dlf_ax_tasks_comm   ! communicator involving the
                                 ! i-th proc from each farm
! **********************************************************************

!!! dummy argument = corresponding variable in the calling program

end subroutine dlf_get_taskfarm


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_output(dum_stdout, dum_stderr)
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,stdout,keep_alloutput
  implicit none
  integer :: dum_stdout
  integer :: dum_stderr
  integer :: ierr
  logical :: topened
  character(len=10) :: suffix

! sort out output units; particularly important on multiple processors
 
! set unit numbers for main output and error messages
  if (dum_stdout >= 0) stdout = dum_stdout 
  if (dum_stderr >= 0) stderr = dum_stderr

  if (glob%iam /= 0) then
     inquire(unit=stdout, opened=topened, iostat=ierr)
     if (topened .and. ierr == 0) close(stdout)
     if (keep_alloutput) then ! hardwired in dlf_global_module.f90
        write(suffix,'(i10)') glob%iam
        open(unit=stdout,file='output.proc'//trim(adjustl(suffix)))
     else
        open(unit=stdout,file='/dev/null')
     end if
  endif

  if (glob%nprocs > 1) then
     ! write some info on the parallelization
     write(stdout,'(1x,a,i10,a)')"I have rank ",glob%iam," in mpi_comm_world"
     write(stdout,'(1x,a,i10)')"Total number of processors = ",glob%nprocs
     if (keep_alloutput) then
        write(stdout,'(1x,a)')"Keeping output from all processors"
     else
        write(stdout,'(1x,a)')"Not keeping output from processors /= 0"
     end if
  end if

end subroutine dlf_output


! **********************************************************************
! **********************************************************************
! The following routine either writes random numbers to a file, or reads
! them. This is to have equal starting conditions for different compilers
subroutine read_rand(arr)
  use dlf_parameter_module, only: rk
  use dlf_global, only : glob
  real(rk)  :: arr(:)
  integer, parameter :: si1=12
  integer, parameter :: si2=3000
  logical,parameter :: readf=.true.
  real(rk) :: ar(si2)
  integer :: l(1),length
  l=ubound(arr)
  length=l(1)
  if(readf) then
    if(length<=si1) then
      open(unit=201,file="random1.bin",form="unformatted")
    else if(length<=si2) then
      open(unit=201,file="random2.bin",form="unformatted")
    else
      call dlf_mpi_finalize() ! only necessary for a parallel build;
                              ! can be present for a serial build
      stop "Too many coordinates to be read from random.bin file"
    end if
    read(201) ar(1:length)
    close(201)
    arr=ar(1:length)
  else
    if (glob%iam == 0) then
       call random_number(ar)
       open(unit=201,file="random1.bin",form="unformatted")
       write(201) ar(1:si1)
       close(201)
       open(unit=201,file="random2.bin",form="unformatted")
       write(201) ar(1:si2)
       close(201)
    end if
  end if
end subroutine read_rand
