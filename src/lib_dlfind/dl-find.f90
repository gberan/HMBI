! **********************************************************************
! **                                                                  **
! **                      DL-FIND main module                         **
! **                                                                  **
! **                    Johannes Kaestner 2006                        **
! **                                                                  **
! **                                                                  **
! **********************************************************************
!!****J* DL-FIND/main
!!
!! NAME
!! DL-FIND
!!
!! FUNCTION
!! Main unit of the optimiser DL-FIND
!!
!!
!! DATA
!! $Date: 2010-05-07 19:30:43 $
!! $Revision: 1.1 $
!! $Author: gberan $
!! $URL: http://ccpforge.cse.rl.ac.uk/svn/dl-find/branches/release_chemsh3.3/dl-find.f90 $
!! $Id: dl-find.f90,v 1.1 2010-05-07 19:30:43 gberan Exp $
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

!! Main layout:
!! external program calls dl_find
!!
!! dl-find
!!   |
!!   dlf_read_in
!!   |
!!   dlf_task (calls dlf_run once or several times)
!!       |
!!       dlf_run 
!!       (restart is at the moment handeled there. This does not work really IMPROVE)
!!          |
!!          dlf_init
!!          main optimisation cycle
!!   |<--   dlf_destroy
!!   |
!!   shut down
!!   |
!! return to calling program
!!
subroutine dl_find(nvarin,nvarin2,nspec,master&
#ifdef GAMESS
    ,core&
#endif
    )
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stdout
  use dlf_stat, only: stat
  use dlf_allocate, only: allocate_report,allocate,deallocate
  use dlf_store, only: store_delete_all
  implicit none
  integer   ,intent(in)    :: nvarin ! number of variables to read in
                                     !  3*nat
  integer   ,intent(in)    :: nvarin2! number of variables to read in
                                     !  in the second array (coords2)
  integer   ,intent(in)    :: nspec  ! number of values in the integer
                                     !  array spec
  integer   ,intent(in)    :: master ! 1 if this task is the master of
                                     ! a parallel run, 0 otherwise
#ifdef GAMESS
  real(rk) :: core(*) ! GAMESS memory, not used in DL-FIND
#endif
  ! 
! **********************************************************************

  ! read input parameters, set defaults
  call dlf_read_in(nvarin,nvarin2,nspec,master)

  ! task manager, main optimisation cycle
  call dlf_task( &
#ifdef GAMESS
    core&
#endif
    )

  ! shut down finally
  call dlf_deallocate_glob()

  ! deallocate arrays in formstep_set_tsmode
  call dlf_formstep_set_tsmode(1,-2,1.d0)

  ! delete dlf_store
  call store_delete_all

  call clock_stop("TOTAL")
  call time_report

  call allocate_report

end subroutine dl_find

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* main/dlf_read_in
!!
!! FUNCTION
!!
!! * Init clock
!! * Get input parameters from calling code
!! * set defaults
!! * call dlf_allocate_glob to allocate arrays globally to the optimiser
!!
!! SYNOPSIS
subroutine dlf_read_in(nvarin,nvarin2,nspec,master)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,pi,stdout,printl,printf
  use dlf_stat, only: stat
  use dlf_allocate, only: allocate,deallocate
  use dlf_store, only: store_initialise
  implicit none
  integer   ,intent(in)    :: nvarin ! number of variables to read in
                                     !  3*nat
  integer   ,intent(in)    :: nvarin2! number of variables to read in
                                     !  in the second array (coords2)
  integer   ,intent(in)    :: nspec  ! number of values in the integer
                                     !  array spec
  integer   ,intent(in)    :: master ! 1 if this task is the master of
                                     ! a parallel run, 0 otherwise
  real(rk),allocatable :: tmpcoords(:),tmpcoords2(:)
  integer, allocatable :: spec(:)
  integer              :: ivar,nat,nframe,nmass,nweight,nz,tsrel
  integer              :: massweight,ierr
  integer              :: tdlf_farm
  integer              :: n_po_scaling
! **********************************************************************

  pi=4.D0*datan(1.D0)

  call time_init
  call clock_start("TOTAL")

  if(nvarin<=0) call dlf_fail("A positive number of variables is needed")

  call allocate(tmpcoords,nvarin)
  call allocate(tmpcoords2,max(nvarin2,1))
  call allocate(spec,max(nspec,1))

  ! ====================================================================
  ! READ EXTERNAL PARAMETERS, SET DEFAULTS
  ! ====================================================================

  ! set parameters to something useless to recognise user input
  call dlf_default_init(nspec,spec)
  nz=0

  ! ====================================================================
  ! The input arrays and their meaning:
  ! tmpcoords(nvarin) : one set of coordinates (main point)
  ! tmpcoords2(nvarin2) : a multi-purpose real array
  !    nframe*nat*3    entries of coordinates of nframe structures
  !    nweight         entries of weights (nat or 0)
  !    nmass           entries of atomic masses (nat or 0)
  !    n_po_scaling    entries of radii scaling factors in the parallel 
  !                    optimization (0 [meaning all radii set to the base 
  !                    value], or a pre-known nivar)
  !    i.e. nvarin2= nframe*nat*3 + nweight + nmass + n_po_scaling
  ! spec(nspec) : integer specification array
  !    nat     entries of freezing/residue number
  !    nz      entries of nuclear charges (same order as coords)
  !    5*ncons entries of constraints (typ, atom1,atom2, atom3, atom4)
  !    2*nconn entries of connections (atom1 atom2)
  !    i.e. nspec= nat + nz + 5*ncons + 2*nconn
  ! ====================================================================
  ! get input parameters
  ivar=1
  massweight=0
  tdlf_farm=1 ! set default value
  n_po_scaling=0 ! set default value
  call dlf_get_params(nvarin,max(nvarin2,1),max(nspec,1), &
      tmpcoords,tmpcoords2,spec, ierr, &
      glob%tolerance,printl,glob%maxcycle,glob%maxene,&
      ivar,glob%icoord,glob%iopt,glob%iline,glob%maxstep, &
      glob%scalestep,glob%lbfgs_mem,glob%nimage,glob%nebk, &
      glob%dump,glob%restart,nz,glob%ncons,glob%nconn,&
      glob%update,glob%maxupd,glob%delta,glob%soft,glob%inithessian, &
      glob%carthessian,tsrel,glob%maxrot,glob%tolrot,nframe,nmass,nweight,&
      glob%timestep,glob%fric0,glob%fricfac,glob%fricp, &
      glob%imultistate, glob%state_i, glob%state_j, &
      glob%pf_c1, glob%pf_c2, glob%gp_c3, glob%gp_c4, glob%ln_t1, glob%ln_t2, &
      printf, glob%tolerance_e, glob%distort, massweight, glob%minstep, &
      glob%maxdump, glob%task, glob%temperature, &
      glob%po_pop_size, glob%po_radius_base, glob%po_contraction, &
      glob%po_tol_r_base, glob%po_tolerance_g, glob%po_distribution, &
      glob%po_maxcycle, glob%po_init_pop_size, glob%po_reset, &
      glob%po_mutation_rate, glob%po_death_rate, glob%po_scalefac, &
      glob%po_nsave,glob%ntasks,tdlf_farm,n_po_scaling)

  ! write(*,*) 'printl = ',printl

  if(ierr/=0) call dlf_fail("Failed to read parameters")

  if (glob%ntasks <= 0) then
    write(stdout,'("glob%ntasks = ",i6)') glob%ntasks
    call dlf_fail("Number of task farms must be positive")
  end if

  ! call this subroutine even if glob%ntasks == 1 to sort out the 
  ! writing to files from each processor and the communicators
  call dlf_make_taskfarm(tdlf_farm)

  ! get logical variables (communication with c-code requires integers)
  glob%tatoms=(ivar==1)
  glob%tsrelative=(tsrel==1)
  glob%massweight=(massweight==1)

  ! Do we need to calculate the interstate coupling gradient?
  if (glob%imultistate > 1) glob%needcoupling = 1

  ! set parameters that have not been set by the user to the default
  ! this routine defines the default values
  call dlf_default_set(nvarin)

  ! check consistency of the array sizes:
  nat=nvarin/3
  if(nspec/=nat+nz+5*glob%ncons+2*glob%nconn) then
    write(stdout,'("nspec ",i6)') nspec
    write(stdout,'("nat   ",i6)') nat
    write(stdout,'("nz    ",i6)') nz
    write(stdout,'("ncons ",i6)') glob%ncons
    write(stdout,'("nconn ",i6)') glob%nconn
    write(stdout,'("nspec should be: nat + nz + 5*ncons + 2*nconn")')
    call dlf_fail("Inconsistent size of array spec - interface error")
  end if

  if(nvarin2/=nframe*nat*3+nweight+nmass+n_po_scaling &
      .or. (nweight/=0.and.nweight/=nat) &
      .or. (nmass/=0.and.nmass/=nat) &
      .or. (n_po_scaling < 0) ) then
    write(stdout,'("nvarin2      ",i6)') nvarin2
    write(stdout,'("nframe       ",i6)') nframe
    write(stdout,'("nat          ",i6)') nat
    write(stdout,'("nweight      ",i6)') nweight
    write(stdout,'("nmass        ",i6)') nmass
    write(stdout,'("n_po_scaling ",i6)') n_po_scaling
    write(stdout,'("varin2 should be: nframe*nat*3 + nweight + nmass + n_po_scaling")')
    write(stdout,'("nweight should be either 0 or nat")')
    write(stdout,'("nmass should be either 0 or nat")')
    call dlf_fail("Inconsistent size of array coords2 - interface error")
  end if

  ! check consistency of multistate calculations
  if (glob%imultistate > 0) call dlf_conint_check_consistency

  ! initialise printout
  if(master==0) then
    ! this is a slave, do not print anything!
    printl=-2
    printf=-2
  end if

  if(printl>=2) call dlf_printheader

  ! ====================================================================
  ! INITIALISE VARIOUS INSTANCES
  ! ====================================================================

  ! allocate storage
  call dlf_allocate_glob(nvarin,nvarin2,nspec,tmpcoords,tmpcoords2,spec,&
      nz,nframe,nmass,nweight,n_po_scaling)
  call deallocate(tmpcoords)
  call deallocate(tmpcoords2)
  call deallocate(spec)

  ! initialise (reset) all statistics counters
  stat%sene=0
  call dlf_stat_reset

  ! initialise dlf_store
  call store_initialise

end subroutine dlf_read_in
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* main/dlf_run
!!
!! FUNCTION
!!
!! * Initialise coordinates, formstep, linesearch
!! * do main optimisation cycle
!! * Destroy coordinates, formstep, linesearch
!!
!! SYNOPSIS
subroutine dlf_run( &
#ifdef GAMESS
    core&
#endif
    )
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stdout,printl,printf
  use dlf_stat, only: stat
  use dlf_allocate, only: allocate,deallocate
  implicit none
#ifdef GAMESS
  real(rk) :: core(*) ! GAMESS memory, not used in DL-FIND
#endif
  integer  :: iimage,ivar,status
  real(rk) :: svar
  logical  :: tconv,trestarted,trerun_energy
  logical  :: needhessian ! do we need a Hessian?
  logical  :: testconv ! is convergence checked in coords_xtoi ?
  logical  :: trestarted_report
! **********************************************************************
  
  iimage=1
  glob%toldenergy_conv=.false.

  ! if thermal analysis (hessian only) is called, make sure that mass-weighted
  ! cartesians are used!
  if(glob%iopt==11) then
    glob%massweight=.true.
    glob%icoord=0
    glob%distort=0.D0
  end if

  ! initialise coordinate transform, allocate memory for it
  call dlf_coords_init
  
  ! initialise search algorithm
  call dlf_formstep_init(needhessian)

  ! initialise line search
  call linesearch_init

  ! prepare for main optimisation cycle

  glob%taccepted=.true.
  tconv=.false.
  trerun_energy=.false.

  ! read checkpoint file 
  if(glob%restart==1) then
    call clock_start("CHECKPOINT")
    call dlf_checkpoint_read(status,trestarted)
    call clock_stop("CHECKPOINT")
    if(.not.trestarted) call dlf_fail("Restart attempt failed")
    trestarted_report=.true.
  else
    trestarted=.false.
    trestarted_report=.false.
  end if

  ! report at the start
  call dlf_report(trestarted_report)
  

  ! open trajectory file
  if(glob%tatoms.or.glob%nimage.gt.1) then
    call clock_start("XYZ")
    if(printf>=3) then
      if(stat%sene>0) THEN
        if (glob%iam == 0) then
           open(unit=30,file="path.xyz",position="APPEND")
           open(unit=31,file="path_active.xyz",position="APPEND")
        end if
      ELSE
        if (glob%iam == 0) then
           open(unit=30,file="path.xyz")
           open(unit=31,file="path_active.xyz")
        end if
      end if
    end if
    call clock_stop("XYZ")
  else
    if (glob%iam == 0) then
       open(unit=30,file="path.inc")
       !open(unit=30,file="path.inc",position="append")
    end if
  end if

  ! ====================================================================
  ! MAIN OPTIMISATION CYCLE
  ! ====================================================================
  do ! exit conditions implemented via exit statements

    if (glob%iopt/10 == 5) then ! parallel optimisation
       call dlf_parallel_opt(trestarted_report, tconv &
#ifdef GAMESS
       ,core&
#endif
       )
       exit ! the MAIN OPTIMISATION CYCLE
    end if

    if(trestarted) goto 1000

    if(.not.trerun_energy) stat%ccycle= stat%ccycle+1
    stat%sene=stat%sene+1
    if(stat%ccycle > glob%maxcycle) then
      stat%ccycle= stat%ccycle-1
      stat%sene=stat%sene-1
      if(printl>0) write(stdout,"(&
          &'Stopping: maximum number of cycles reached')")
      exit
    end if

    ! ==================================================================
    ! EVALUATE THE ENERGY
    ! ==================================================================

    ! get out of main cycle if energy evaluated more often than glob%maxene
    if(stat%sene > glob%maxene) then
      stat%ccycle= stat%ccycle-1
      stat%sene=stat%sene-1
      if(printl>=2) write(stdout,"(&
          &'Stopping: maximum number of energy evaluations reached')")
      exit
    end if

    if(printl>=6) write(stdout,"('Calculating the energy',i5)") stat%sene

    call clock_start("EANDG")

    if (glob%imultistate == 0) then
       ! An ordinary single state gradient
       !write(*,*) 'GJB DLF: getting the energy & gradient'
       !write(*,*) 'GJB Before the call: xcoords(1,1) = ', glob%xcoords(1,1)
       call dlf_get_gradient(glob%nvar,glob%xcoords,glob%energy, &
            glob%xgradient,iimage,&
#ifdef GAMESS
            core,&
#endif
            status)
       !write(*,*) 'GJB xcoords(1,1) after get_gradient= ', glob%xcoords(1,1)
       !write(*,*) 'GJB gradient = ',glob%xgradient
    else
       ! Multiple state gradient calculation
       if (printl >= 6) write(stdout, '(a)') "Calculating multistate energies"
       call dlf_get_multistate_gradients(glob%nvar,glob%xcoords,glob%msenergy, &
            glob%msgradient,glob%mscoupling,glob%needcoupling,iimage,status)
       ! Form the objective function and gradient from the individual
       ! state gradients
       if (printl >= 4) then
          write(stdout, '(a, f20.10)') "Lower state energy: ", glob%msenergy(1)
          write(stdout, '(a, f20.10)') "Upper state energy: ", glob%msenergy(2)
          write(stdout, '(a, f20.10)') "Energy difference:  ", &
               abs(glob%msenergy(1) - glob%msenergy(2))
       endif
       if (printl >= 6) write(stdout, '(a)') "Forming objective function"
       call dlf_make_conint_gradient
    endif

    !write(*,*) 'Just got HMBI energy and gradient'
    !write(*,*) 'energy = ',glob%energy

    call clock_stop("EANDG")

    ! check of NaN in the energy (comparison of NaN with any number 
    ! returns .false. , pgf90 does not understand isnan() )
    if( abs(glob%energy) > huge(1.D0) ) then
      status=1
    else
      if (.not. abs(glob%energy) < huge(1.D0) ) status=1
    end if

    if(status/=0) then
      call dlf_report(trestarted_report)
      call dlf_fail("Energy evaluation failed")
    end if

    if(iimage==1) then
      if(printl>=2) write(stdout,'(1x,a,es16.9)') &
          "Energy calculation finished, energy: ", &
          glob%energy
    else
      if(printl>=2) write(stdout,'(1x,a,i4,a,es16.9)') &
          "Energy calculation of image ",iimage,&
          " finished, energy: ",glob%energy
    end if

    ! send coordinates to the calling program
    if(printf>=1) then
      !write(*,*) 'GJB xcoords(1) = ', glob%xcoords
      call dlf_put_coords(glob%nvar,1,glob%energy,glob%xcoords,glob%iam)
    end if

    ! write restart information
    if(stat%sene<=glob%maxdump) then
      if(glob%dump>0) then
        if(mod(stat%sene,glob%dump)==0) then
          call clock_start("CHECKPOINT")
          if(printl>=6) write(stdout,"('Writing restart information')")
          call dlf_checkpoint_write(status)
          call clock_stop("CHECKPOINT")
        end if
      end if
    end if
    ! come here if checkpoint file successfully read
1000 trestarted=.false.

    ! ==================================================================
    ! TRANSFORM CARTESIANS TO INTERNALS (COORDS AND GRAD)
    ! ==================================================================
    call clock_start("COORDS")
    call dlf_coords_xtoi(trerun_energy,testconv,iimage)
    call clock_stop("COORDS")
 
    ! ==================================================================
    ! Multistate gradient: post-transformation
    ! ==================================================================
    if (glob%imultistate == 3) then
       call clock_start("EANDG")
       call dlf_make_ln_gradient_posttrans
       call clock_stop("EANDG")
    endif

    ! ==================================================================
    ! Calculate the Hessian if necessary
    ! ==================================================================
    ! glob%ihessian is allocated and deallocated in formstep_init and
    !  formstep_destroy
    if(needhessian) then
      ! Make sure the Hessian exists. This can imply energy cycles (trerun_energy)
      ! or a call to dlf_get_hessian
      call dlf_makehessian(trerun_energy,tconv &
#ifdef GAMESS
          ,core&
#endif
          )
      if (tconv) exit
    end if

    ! this may be parallelised in the future ...
    ! ==================================================================
    if(trerun_energy) cycle ! main optimisation cycle
    ! ==================================================================

    ! exit if only the hessian and thermal analysis should be calculated
    if(glob%iopt == 11) then
      call dlf_thermal
      exit ! no optimisation cycles
    end if

    ! write trajectory
    if(printf>=3 .and. glob%iam == 0) then
      call write_xyz(30,glob%nat,glob%znuc,glob%xcoords)
      call write_xyz_active(31,glob%nat,glob%znuc,glob%spec,glob%xcoords)
    end if

    ! if trust-radius, test for step acceptance. 
    ! If rejected, do not form a new step and keep old energy.
    if (glob%iline==1) then
      call test_acceptance
    end if
    if (glob%iline==2) then
      call test_acceptance_g
    end if
    if (glob%iline==3) then
      call linesearch
    end if

    ! ==================================================================
    ! TEST FOR CONVERGENCE
    ! ==================================================================
    if(glob%taccepted) then
      stat%caccepted=stat%caccepted+1
      tconv=.false.
      if(.not.testconv) call convergence_test(stat%ccycle,.true.,tconv)
      if(tconv) exit
    end if

    ! ==================================================================
    ! FORM AN OPTIMISATION STEP
    ! ==================================================================
    if(glob%taccepted) then
      call clock_start("FORMSTEP")
      call dlf_formstep
      call clock_stop("FORMSTEP")
    end if

    ! ==================================================================
    ! DO A LINE SEARCH OR A TRUST RADIUS STEP
    ! ==================================================================

    if(glob%taccepted) then
      
      ! scale the step 
      call dlf_scalestep

      call clock_start("COORDS")

      ! check the step in special cases
      ! Set step of frozen NEB images to zero
      if(glob%icoord/100==1) call dlf_neb_checkstep

      ! Check step in case of dimer
      if(glob%icoord/100==2) call dlf_dimer_checkstep
      
      call clock_stop("COORDS")

      ! do the step
      glob%icoords(:)=glob%icoords(:) + glob%step(:)

    end if

    if(glob%taccepted) then
      ! store old values
      glob%oldenergy=glob%energy
      glob%toldenergy=.true.
    end if

    ! ==================================================================
    ! TRANSFORM INTERNAL COORDINATES TO CARTESIANS (COORDS ONLY)
    ! ==================================================================

    call clock_start("COORDS")
    call dlf_coords_itox(iimage)
    call clock_stop("COORDS")

  end do ! main simulation cycle

  ! Job finished, prepare for shutdown
  if(glob%iopt /= 11) then
    if(tconv) then
      if(printl>=2)write(stdout,'(1x,a)') "converged"
    else
      if(printl>=2)write(stdout,'(1x,a)') "NOT CONVERGED"
    end if
  end if

  ! close trajectory files
  if(printf>=3 .and. glob%iam == 0) then
    close(30)
    close(31)
  end if
  ! again some rubbish for printing ...
  if(glob%icoord/100==1.and..not.glob%tatoms) then
    if (glob%iam == 0) then
       open(unit=30,file="path.inc")
       do ivar=1,glob%nimage
         write(30,"('sphere{<',f12.4,',',f12.4,',',f12.4,'> 0.04}')") &
             glob%icoords(ivar*2-1),1.D0,glob%icoords(ivar*2)
         if(ivar>1) then
           write(30,"('cylinder{<',f12.4,',',f12.4,',',f12.4,'>,"//&
               &"<',f12.4,',',f12.4,',',f12.4,'> 0.01} ')") &
               glob%icoords(ivar*2-1),1.D0,glob%icoords(ivar*2), &
               glob%icoords((ivar-1)*2-1),1.D0,glob%icoords((ivar-1)*2)
         end if
       end do

       close(30)
    end if
  end if

  ! ====================================================================
  ! CLOSE DOWN
  ! ====================================================================

  ! delete memory for line search 
  if (glob%iline==1.or.glob%iline==2.or.glob%iline==3) then
    call linesearch_destroy
  end if

  ! delete memory for search algorithm
  call dlf_formstep_destroy

  ! delete memory for internal coordinates
  call dlf_coords_destroy

  !write report on this optimisation
  call dlf_report(trestarted_report)

end subroutine dlf_run
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* main/dlf_allocate_glob
!!
!! FUNCTION
!!
!! Allocate the arrays globally available to the optimiser,
!! set a few starting values.
!!
!! SYNOPSIS
subroutine dlf_allocate_glob(nvarin,nvarin2,nvarspec, &
    tmpcoords,tmpcoords2,spec,nz,nframe,nmass,nweight,n_po_scaling)
!! SOURCE
  use dlf_parameter_module, only: rk,ik
  use dlf_global, only: glob,stdout,printl
  use dlf_allocate, only: allocate
  implicit none
  integer, intent(in) :: nvarin
  integer, intent(in) :: nvarin2
  integer, intent(in) :: nvarspec
  real(rk),intent(in) :: tmpcoords(nvarin)
  real(rk),intent(in) :: tmpcoords2(nvarin2)
  integer, intent(in) :: spec(nvarspec)
  integer, intent(in) :: nz
  integer, intent(in) :: nframe
  integer, intent(in) :: nmass
  integer, intent(in) :: nweight
  integer, intent(in) :: n_po_scaling
  integer :: nat,ivar
! **********************************************************************
  if(glob%tinit) return ! this instance has been initialised
  glob%tcoords2=(nframe>0)
  if(glob%tatoms) then
    ! input contains atoms
    if(mod(nvarin,3)/=0) call dlf_fail("nvarin has to be 3*nat")
    nat=nvarin/3
    glob%nvar=nvarin
    if(printl>4) write(stdout,'(a,i5,a)') &
      "Input contains ",nat," atoms"
    glob%nat=nat

    call allocate( glob%xcoords,3,nat)

    if(glob%tcoords2) then
      call allocate( glob%xcoords2,3,nat,nframe)
      glob%xcoords2 = reshape(tmpcoords2(1:nat*3*nframe),(/3,nat,nframe/))
    end if

    call allocate( glob%xgradient,3,nat)
    call allocate( glob%spec,nat)
    !nuclear charges
    call allocate( glob%znuc,nat)
    if(nz==nat) then
      glob%znuc(:)=spec(nat+1:nat+nz)
    else
      glob%znuc(:)=1
    end if
    
    glob%xcoords = reshape(tmpcoords,(/3,nat/))
    glob%spec(:) = spec(1:nat)

    if(glob%ncons>0) then
      call allocate( glob%icons,5,glob%ncons)
      ivar=nat+nz
      glob%icons=reshape(spec(ivar+1:ivar+5*glob%ncons),(/5,glob%ncons/))
    else
      call allocate( glob%icons,5,1)
      glob%icons(:,:)=0
    end if

    ! user input connections - additional to the hdlc-primitive created
    if(glob%nconn>0) then
      call allocate( glob%iconn,2,glob%nconn)
      ivar=nat+nz+glob%ncons
      glob%iconn=reshape(spec(ivar+1:ivar+2*glob%nconn),(/2,glob%nconn/))
    else
      call allocate( glob%iconn,2,1)
      glob%iconn(:,:)=0
    end if

    call allocate( glob%weight,nat)
    if(nweight>0) then
      glob%weight=tmpcoords2(nat*3*nframe+1:nat*3*nframe+nat)
    else
      glob%weight(:)=1.D0
    end if
    call allocate( glob%mass,nat)
    if(nmass>0) then
      glob%mass=tmpcoords2(nat*3*nframe+nweight+1:nat*3*nframe+nweight+nat)
    else
      glob%mass(:)=1.D0
    end if

    if (glob%iopt/10==5) then ! parallel optimizers
       call allocate(glob%po_radius_scaling,max(1,n_po_scaling))
       if (n_po_scaling == 0) then 
          ! this is a shorthand for setting all elements of the glob%po_radius(:) and 
          ! glob%po_tolerance_r(:) arrays to their respective base values.
          glob%po_radius_scaling(:) = 1.0D0 ! there should only be one component anyway
       else
          ! n_po_scaling < 0 will be caught in dlf_read_in
          glob%po_radius_scaling(:) = tmpcoords2(nat*3*nframe + nweight + nmass + 1 : &
                                      nat*3*nframe + nweight + nmass + n_po_scaling)
       end if 
    end if

    ! multistate arrays for conical intersection search
    if (glob%imultistate > 0) then
       call allocate(glob%msenergy, 2)
       call allocate(glob%msgradient, 3, nat, 2)
       if (glob%needcoupling == 1) then
          call allocate(glob%mscoupling, 3, nat)
       else
          ! dummy array just to avoid problems in argument list
          ! of dlf_get_multistate_gradients - I assume you can't
          ! send unallocated arrays to a C subroutine
          call allocate(glob%mscoupling, 1, 1)
       endif
    endif

  else
    ! input contains sequential variables (as the 2D-potentials)
    glob%nat=-1
    glob%nvar=nvarin
    if(printl>4) write(stdout,'(a,i5,a)') &
      "Input contains ",nvarin," degrees of freedom (no atoms)"
    call allocate( glob%xcoords,1,nvarin)
    if(glob%tcoords2) then
      call allocate( glob%xcoords2,1,nvarin,nframe)
      glob%xcoords2 = reshape(tmpcoords2,(/1,nvarin,nframe/))
    end if
    call allocate( glob%xgradient,1,nvarin)
    ! NOT TO BE USED
    call allocate( glob%spec,1)
    call allocate( glob%icons,5,1)
    glob%icons(:,:)=0
    call allocate( glob%znuc,1)
    call allocate( glob%iconn,2,1)
    glob%iconn(:,:)=0
    call allocate( glob%weight,1)
    call allocate( glob%mass,1)

    glob%xcoords = reshape(tmpcoords,(/1,nvarin/))

  end if
  glob%toldenergy=.false.
  glob%tinit=.true.
end subroutine dlf_allocate_glob
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* main/dlf_deallocate_glob
!!
!! FUNCTION
!!
!! Deallocate global arrays
!!
!! SYNOPSIS
subroutine dlf_deallocate_glob
!! SOURCE
  use dlf_parameter_module, only: rk,ik
  use dlf_global, only: glob
  use dlf_allocate, only: deallocate
  implicit none
! **********************************************************************
  if(.not.glob%tinit) return ! this instance has not been initialised
  call deallocate( glob%xcoords )
  if(glob%tcoords2) then
    call deallocate( glob%xcoords2 )
  end if
  call deallocate( glob%xgradient )

  ! arrays only used for atoms
  call deallocate( glob%spec  )
  call deallocate( glob%icons )
  call deallocate( glob%znuc  )
  call deallocate( glob%iconn )

  call deallocate( glob%mass )
  call deallocate( glob%weight )

  glob%tinit=.false.

  ! conical intersection search
    if (glob%imultistate > 0) then
       call deallocate(glob%msenergy)
       call deallocate(glob%msgradient)
       call deallocate(glob%mscoupling)
    endif

  ! parallel optimisation algorithms
  if (glob%iopt/10==5) call deallocate(glob%po_radius_scaling)

end subroutine dlf_deallocate_glob
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* main/dlf_fail
!!
!! FUNCTION
!!
!! Shut down DL-FIND after a severe error. Calls dlf_error, which has
!! to be provided by the calling code.
!!
!! SYNOPSIS
subroutine dlf_fail(msg)
!! SOURCE
  use dlf_global, only: stdout, stderr
  implicit none
  character(*),intent(in) :: msg
  call flush(stdout)
  call flush(stderr)
  write(stderr,"(/,a,/,a,/)") "DL-FIND ERROR:",msg
  write(stdout,"(/,a,/,a,/)") "DL-FIND ERROR:",msg
  call flush(stdout)
  call flush(stderr)
  call dlf_error()
  ! this should not be reached, as dlf_error should not return
  stop
end subroutine dlf_fail
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* main/dlf_report
!!
!! FUNCTION
!!
!! Report information about the run to stdout. Called at the beginning 
!! and at the end of a calculations.
!!
!! SYNOPSIS
subroutine dlf_report(trestarted)
!! SOURCE
  use dlf_global, only: glob,stdout,printl
  use dlf_stat, only: stat
  implicit none
  logical ,intent(in) :: trestarted
  integer             :: ia3(3)
  integer             :: kk

! **********************************************************************


  ! Can't return as dlf_mpi_counters requires communication
  ! if(printl<=0) return
  if (printl <= 0) goto 111 ! to skip the report printing up to the 
                            ! call to dlf_mpi_counters

  write(stdout,'(/,a)') "DL-FIND Report:"
  write(stdout,'(a)')   "==============="

  ! ====================================================================
  ! Multistate calculations
  ! ====================================================================
  if (glob%imultistate == 1) then
     write(stdout,'(a)') &
          "Conical intersection algorithm: penalty function"
     write(stdout,1000) "Penalty function constant c1", &
          glob%pf_c1
     write(stdout,1000) "Penalty function constant c2", &
          glob%pf_c2
  end if
  if (glob%imultistate == 2) then
     write(stdout,'(a)') &
       "Conical intersection algorithm: gradient projection"
     write(stdout,1000) "Gradient projection constant c3", &
          glob%gp_c3
     write(stdout,1000) "Gradient projection constant c4", &
          glob%gp_c4
  end if
  if (glob%imultistate == 3) then
     write(stdout,'(a)') &
       "Conical intersection algorithm: Lagrange-Newton"
     write(stdout,1000) "Lagrange-Newton threshold t1", &
          glob%ln_t1
     write(stdout,1000) "Lagrange-Newton threshold t2", &
          glob%ln_t2
  end if
  if (glob%imultistate > 0) then
     write(stdout,2000) "Lower electronic state", &
          glob%state_i
     write(stdout,2000) "Upper electronic state", &
          glob%state_j
  end if

  ! ====================================================================
  ! Optimiser
  ! ====================================================================
  if(glob%iopt==0) write(stdout,'(a)') &
      "Optimisation algorithm: Steepest descent"
  if(glob%iopt==1) write(stdout,'(a)') &
      "Optimisation algorithm: Conjugate gradient"
  if(glob%iopt==2) write(stdout,'(a)') &
      "Optimisation algorithm: Conjugate gradient"
  if(glob%iopt==3) then
    write(stdout,'(a)') "Optimisation algorithm: L-BFGS"
    write(stdout,2000) "Number of steps in L-BFGS memory", &
         glob%lbfgs_mem
  end if
  if (glob%iopt==51) write(stdout,'(a)') &
      "Optimisation algorithm: Stochastic search"
  if (glob%iopt==52) write(stdout,'(a)') &
      "Optimisation algorithm: Genetic algorithm"

  if(glob%iopt==10) write(stdout,'(a)') &
      "Optimisation algorithm: P-RFO"
  if(glob%iopt==11) then 
    write(stdout,'(a)') &
        "Optimisation algorithm: no optimisation, just calculate the"
    write(stdout,'(a)') "    Hessian and thermal corrections"
  end if
  if(glob%iopt==20) write(stdout,'(a)') &
      "Optimisation algorithm: Newton-Raphson"
  if(glob%iopt==30) write(stdout,'(a)') &
      "Optimisation algorithm: Damped dynamics"
  if(glob%iopt==40) write(stdout,'(a)') &
      "Optimisation algorithm: Lagrange-Newton"

  ! Hessian options
  if(glob%iopt==10.or.glob%iopt==20.or.glob%iopt==40) then
    if (glob%inithessian == 0) write(stdout,'(a)') &
         "Initial Hessian from external program"
    if (glob%inithessian == 1) write(stdout,'(a)') &
         "Initial Hessian by one-point finite difference"
    if (glob%inithessian == 2) write(stdout,'(a)') &
         "Initial Hessian by two-point finite difference"
    if (glob%inithessian == 3) write(stdout,'(a)') &
         "Initial Hessian by diagonal one-point finite difference"
    if (glob%inithessian == 4) write(stdout,'(a)') &
         "Initial Hessian is identity matrix"
    if(glob%update==0) write(stdout,'(a)') &
      "No Hessian updates"
    if(glob%update==1) write(stdout,'(a)') &
      "Hessian update mechanism: Powell"
    if(glob%update==2) write(stdout,'(a)') &
      "Hessian update mechanism: Bofill"
    if(glob%update==3) write(stdout,'(a)') &
      "Hessian update mechanism: BFGS"
    write(stdout,2000) "Maximum Number of Hessian updates before recalc.", &
         glob%maxupd
    write(stdout,1000) "Finite difference for Hessian calculation", &
         glob%delta
    if(glob%soft>0.D0) then
      write(stdout,1000) "Eigenmodes below this value are considered soft", &
          glob%soft
    else
      write(stdout,'(a)') "No eigenmodes are considered soft"
    end if
    write(stdout, 1000) "Minimum step size for Hessian update", &
         glob%minstep
     
  end if
  if(glob%iopt==11) then
    write(stdout,1000) "Finite difference for Hessian calculation", &
         glob%delta
  end if

  ! Damped dynamics options
  if(glob%iopt==30) then
    write(stdout,1000) "Time step (a.u.)", &
          glob%timestep
    write(stdout,1000) "Start friction", &
          glob%fric0
    write(stdout,1000) "Friction decreasing factor if energy decreases", &
          glob%fricfac
    write(stdout,1000) "Friction to apply if energy increases", &
          glob%fricp
  end if

  ! Parallel optimisation options
  if(glob%iopt/10==5) then
    write(stdout,2000) "Working population size", &
          glob%po_pop_size
    write(stdout,1000) "Base sample radius", &
          glob%po_radius_base
    write(stdout,1000) "Tolerance on max component of mod g", &
          glob%po_tolerance_g
    write(stdout,2000) "Maximum number of cycles", &
          glob%po_maxcycle
    do kk = 1, SIZE(glob%po_radius_scaling,1)
       write(stdout,1000) "Scaling factor for radius component", &
       glob%po_radius_scaling(kk)
    end do
    if (glob%iopt==51) then
       write(stdout,'(a)') "Stochastic-search-specific options:"
       write(stdout,1000) "Radius contraction factor", &
          glob%po_contraction
       write(stdout,1000) "Base tolerance on radius", &
          glob%po_tol_r_base
       if (glob%po_distribution == 1) write(stdout,'(a)') &
          "Search strategy: uniform"
       if (glob%po_distribution == 2) write(stdout,'(a)') &
          "Search strategy: force_direction_bias"
       if (glob%po_distribution == 3) then 
          write(stdout,'(a)') "Search strategy: force_bias"
          write(stdout,1000) "Scaling factor for absolute gradient vector", &
             glob%po_scalefac
       end if
    else if (glob%iopt==52) then
       write(stdout,'(a)') "Genetic-algorithm-specific options:"
       write(stdout,2000) "Initial population size", &
          glob%po_init_pop_size
       write(stdout,2000) "Number of cycles before resetting population", &
          glob%po_reset
       write(stdout,1000) "Mutation rate", &
          glob%po_mutation_rate
       write(stdout,1000) "Death rate", &
          glob%po_death_rate
       write(stdout,2000) "Number of low-energy minima to store", &
          glob%po_nsave
    end if
  end if

  write(stdout,*)

  ! ====================================================================
  ! Line search / trust radius
  ! ====================================================================
  if (glob%iopt/10 /= 5) then ! we're not running a parallel optimisation 
     if(glob%iline==0) write(stdout,'(a)') &
         "Step length: simple scaled"
     if(glob%iline==1) write(stdout,'(a)') &
         "Trust radius based on energy"
     if(glob%iline==2) write(stdout,'(a)') &
         "Trust radius based on the gradient"
     if(glob%iline==3) write(stdout,'(a)') &
         "Line search"
     write(stdout,1000) "Maximum step length", glob%maxstep
     if(glob%iline==0 .or. glob%iopt == 0) write(stdout,1000) &
          "Scaling step by", glob%scalestep

     write(stdout,*)
  end if 

  ! ====================================================================
  ! Coordinate system
  ! ====================================================================

  ! direct coordinate system
  if(mod(glob%icoord,10)==0) then
    if(glob%massweight) then
      write(stdout,'(a)') "Coordinate system: Mass-weighted Cartesian coordinates"
    else
      write(stdout,'(a)') "Coordinate system: Cartesian coordinates"
    end if
  end if
  if(mod(glob%icoord,10)==1) write(stdout,'(a)') &
      "Coordinate system: Hybrid delocalised internal coordinates (HDLC)"
  if(mod(glob%icoord,10)==2) write(stdout,'(a)') &
      "Coordinate system: Hybrid delocalised total connection scheme (HDLC-TC)"
  if(mod(glob%icoord,10)==3) write(stdout,'(a)') &
      "Coordinate system: Delocalised internal coordinates (DLC)"
  if(mod(glob%icoord,10)==4) write(stdout,'(a)') &
      "Coordinate system: Delocalised total connection scheme (DLC-TC)"

  ! multi-image approaches
  if(glob%icoord/10==10) write(stdout,'(a)') &
        "Nudged elastic band with minimised start and endpoint"
  if(glob%icoord/10==11) write(stdout,'(a)') &
        "Nudged elastic band with start and endpoint minimised perpendicular to path"
  if(glob%icoord/10==12) write(stdout,'(a)') &
        "Nudged elastic band with frozen start and endpoint"
  if(glob%icoord>=100.and.glob%icoord<200) then
    write(stdout,2000) "Number of images",glob%nimage
    write(stdout,1000) "NEB spring constant",glob%nebk
  end if
  ! Dimer method
  if(glob%icoord/100==2) then
    if(glob%icoord/10==20) write(stdout,'(a)') &
        "Dimer method, rotation by optimiser"
    if(glob%icoord/10==21) write(stdout,'(a)') &
        "Dimer method, rotation by line search without extrapolation"
    if(glob%icoord/10==22) write(stdout,'(a)') &
        "Dimer method, rotation by line search with extrapolation"
    write(stdout,1000) "Dimer distance (mid- to endpoint)",glob%delta
    write(stdout,2000) "Maximum number of dimer rotations",glob%maxrot
    write(stdout,1000) "Tolerance in rotation",glob%tolrot,"degrees"
    if(.not.glob%tcoords2) write(stdout,'(a)') &
        "Initial dimer direction randomised"

  end if

  write(stdout,*)

  ! ====================================================================
  ! System size
  ! ====================================================================

  if(glob%tatoms) &
      write(stdout,2000) "Number of atoms",glob%nat 
  if(glob%tcoords2) then
    ia3=shape(glob%xcoords2)
  else
    ia3=(/1,1,1/)
  end if
  ! input geometries are coords and coords2
  write(stdout,2000) "Number of input geometries",ia3(3)+1
  if(glob%tcoords2.and.abs(glob%distort)>0.D0) then
    write(stdout,1000) "Distorting start coordinates along coords2 by",&
        glob%distort
  end if

  write(stdout,2000) "Variables to be optimised",glob%nivar
  if(glob%dump>0) then
    write(stdout,2000) "Restart information is written every", &
        glob%dump,"steps"
  else
    write(stdout,'(a)') "No restart information is written"
  end if
  if(trestarted) then
    write(stdout,'(a)') "This run has been restarted from files."
  else
    write(stdout,'(a)') "This run has not been restarted."
  end if

  if(stat%ccycle > 0) then
    write(stdout,2000) "Number of energy evaluations on this processor",stat%sene
    write(stdout,2000) "Number of steps",stat%ccycle
    write(stdout,2000) "Number of accepted steps / line searches",stat%caccepted
  else
    write(stdout,2000) "Maximum number of steps",glob%maxcycle
    write(stdout,2000) "Maximum number of energy evaluations",glob%maxene
  end if

111 continue ! skip to here if printl is less than or equal to zero

  if (glob%ntasks > 1) call dlf_mpi_counters()

  if (printl > 0) write(stdout,*)
  call flush(stdout)

  ! real number
  1000 format (t1,'................................................', &
           t1,a,' ',t50,es10.3,1x,a)
  ! integer
  2000 format (t1,'................................................', &
           t1,a,' ',t50,i10,1x,a)
end subroutine dlf_report
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* main/dlf_printheader
!!
!! FUNCTION
!!
!! Print header to stdout
!!
!! SYNOPSIS
subroutine dlf_printheader
!! SOURCE
  use dlf_global, only: stdout
  implicit none
  character(50)   :: svn_v
  call dlf_svnversion(svn_v)
  if(svn_v=="") then
    ! svn not available, use the revision number of this file
    svn_v="$Revision: 1.1 $"
    svn_v= svn_v(11:)
  end if

  write(stdout,'(a)') &
      "***********************************************************************"      
  write(stdout,'(a)') &
      "**                                                                   **"      
  write(stdout,'(a)') &
      "**                              DL-FIND                              **"      
  write(stdout,'(a)') &
      "**                       Geometry Optimisation                       **"      
  write(stdout,'(a)') &
      "**               written by Johannes Kaestner, in 2006               **"      
  write(stdout,'(a)') &
      "**               Copyright:  STFC Daresbury Laboratory               **"      
!                    g  f  e  d  c  b  a  C  a  b  c  d  e  f  g
!  write(stdout,'(a)') &
!      "**                          $Revision: 1.1 $                         **"      
  write(stdout,'("**",27x,"Revision: ",a30,"**")') svn_v
  write(stdout,'(a)') &
      "**                                                                   **"      
  write(stdout,'(a)') &
      "***********************************************************************"      
  call flush(stdout)
end subroutine dlf_printheader
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* main/dlf_default_init
!!
!! FUNCTION
!!
!! set input parameters to something useless to recognise user input
!! a few parameters are set here, see comments in the code below
!!
!! SYNOPSIS
subroutine dlf_default_init(nspec,spec)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,printl,printf
  implicit none
  integer,intent(in) :: nspec
  integer,intent(out):: spec(nspec)
! **********************************************************************
  printl=-1
  printf=-1
  glob%maxcycle=-1
  glob%maxene=-1
  glob%tolerance=-1.D0
  glob%tolerance_e=-1.D0  
  glob%toldenergy=.false.
  glob%tinit=.false.
  glob%tatoms=.false.
  glob%iopt=-1
  glob%iline=-1
  glob%maxstep=-1.D0
  glob%scalestep=-1.D0
  glob%maxdump=-1
  
  glob%lbfgs_mem=-1
  glob%update=-1
  glob%maxupd=-1
  glob%havehessian=.false.
  glob%delta=-1.D0
  glob%soft=1.D20
  glob%inithessian = -1
  glob%carthessian=-1
  glob%tsrelative=.false.
  glob%minstep=-1.D0
  
  glob%icoord=-1
  glob%nimage=-1
  glob%nebk=-1.D0
  glob%maxrot=-1
  glob%tolrot=-1.D20
  spec(:)=0 ! for spec, the default is set here (as it may contain
            ! positive and negative values)
  glob%ncons=-1
  glob%nconn=-1
  glob%dump=-1
  glob%restart=-1

  glob%timestep=-1.D0
  glob%fric0=-1.D0
  glob%fricfac=-1.D0
  glob%fricp=-1.D0

  glob%imultistate = 0
  glob%needcoupling = 0
  glob%state_i = 1
  glob%state_j = 2
  glob%pf_c1 = 5.0d0
  glob%pf_c2 = 5.0d0
  glob%gp_c3 = 1.0d0
  glob%gp_c4 = 0.9d0
  glob%ln_t1 = 1.0d-4
  glob%ln_t2 = 1.0d0

  glob%distort = 0.D0 ! default given here

  glob%task = -1
  
  glob%temperature = -1.D0

  glob%po_pop_size= -1
  glob%po_radius_base= -1.0D0
  glob%po_contraction= -1.0D0
  glob%po_tol_r_base= -1.0D0
  glob%po_tolerance_g= -1.0D0
  glob%po_distribution= -1
  glob%po_maxcycle= -1
  glob%po_init_pop_size= -1
  glob%po_reset= -1
  glob%po_mutation_rate= -1.0D0
  glob%po_death_rate= -1.0D0
  glob%po_scalefac= -1.0D0
  glob%po_nsave= -1

  glob%ntasks= -1 ! why???  no real need for this not to be 1...

end subroutine dlf_default_init
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* main/dlf_default_set
!!
!! FUNCTION
!!
!! set parameters that have not been set by the user to the default
!! this routine defines the default values
!!
!! SYNOPSIS
subroutine dlf_default_set(nvarin)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,printl,printf,stdout
  implicit none
  integer,   intent(in)  :: nvarin
! **********************************************************************
  if(printl < 0) printl=2
  if(printf<0) printf=2
  if(glob%maxcycle < 0) glob%maxcycle=100
  if(glob%maxene < 0) glob%maxene=100000

  if(glob%tolerance < 0.D0) glob%tolerance=4.5D-4
  if(glob%tolerance_e < 0.D0) glob%tolerance_e= glob%tolerance/ 450.D0
  if(glob%iopt < 0) glob%iopt=3 ! L-BFGS

  if(glob%iline < 0) then
    if(glob%iopt < 3 .or. glob%imultistate == 2) then
      glob%iline=2
    else
      glob%iline=0
      if(glob%iopt==3.and.glob%icoord<10) glob%iline=1
    end if
  end if
  if(glob%maxstep < 0.D0) glob%maxstep=0.5D0
  if(glob%scalestep < 0.D0) then
    !if(glob%iopt==3) then
      glob%scalestep=1.D0
    !else
    !  glob%scalestep=0.2D0
    !end if
  end if
  if(glob%lbfgs_mem<0) then
    glob%lbfgs_mem=max(nvarin,5)
    glob%lbfgs_mem=min(nvarin,50)
  end if
  if(glob%update<0) glob%update=2
  if(glob%maxupd<0) glob%maxupd=50
  if(glob%delta<0.D0) glob%delta=0.01D0
  if(glob%inithessian == -1) glob%inithessian = 0
  if(glob%carthessian==-1) glob%carthessian=0
  if(glob%minstep < 0.D0) glob%minstep = 1.D-5
  if(glob%maxdump < 0 ) glob%maxdump=100000

  if(glob%nimage<0) then
    if(glob%icoord>=100.and.glob%icoord<200) then
      !NEB
      glob%nimage=10
    else
      glob%nimage=1
    end if
  end if
  if(glob%icoord<0) then
    if(glob%nimage==1) then
      glob%icoord=0
    else
      glob%icoord=110 ! NEB with endpoints perpendicular to tau
    end if
  end if
  if(abs(glob%soft)>1.D19) then
    if ((glob%icoord >= 3 .and. glob%icoord <= 4) .or. &
        (glob%icoord >= 13 .and. glob%icoord <= 14)) then
      ! in case of internals only, there is nothing like soft modes
      glob%soft=-1.D0
    else
      glob%soft=5.D-3
    end if
  end if
  if(glob%nebk<0.D0) glob%nebk=0.01D0
  if(glob%maxrot<0) glob%maxrot=10
  if(glob%tolrot<-1.D19) glob%tolrot=5.D0

  ! for "spec" the default is set in dlf_default_init!
  if(glob%ncons<0) glob%ncons=0
  if(glob%nconn<0) glob%nconn=0
  if(glob%dump<0) glob%dump=0
  if(glob%restart<0) glob%restart=0
  ! glob%weight and mass defaults are set in dlf_allocate_glob

  if(glob%timestep <=0.D0) glob%timestep= 1.D0
  if(glob%fric0 <=0.D0) glob%fric0= 0.3D0
  if(glob%fricfac <=0.D0) glob%fricfac=0.95D0
  if(glob%fricp <=0.D0) glob%fricp=0.3D0

  ! Conical intersection search parameters
  ! Note that <= 0 is allowed for ln_t1/ln_t2
  if (glob%pf_c1 .le. 0.D0) glob%pf_c1 = 5.0d0
  if (glob%pf_c2 .le. 0.D0) glob%pf_c2 = 5.0d0
  if (glob%gp_c3 .le. 0.D0) glob%gp_c3 = 1.0d0
  if (glob%gp_c4 .le. 0.D0) glob%gp_c4 = 0.9d0

  ! the default for glob%distort is in dlf_default_init

  if (glob%task < 0) glob%task = 0

  if(glob%temperature < 0.D0) glob%temperature=300.D0

  if (glob%po_pop_size < 1) glob%po_pop_size = 25
  if (glob%po_radius_base <= 0.0D0) glob%po_radius_base = 1.0D0
  if (glob%po_contraction <= 0.0D0) glob%po_contraction = 0.9D0
  if (glob%po_tol_r_base <= 0.0D0) glob%po_tol_r_base = 1.0D-8
  if (glob%po_tolerance_g <= 0.0D0) glob%po_tolerance_g = 1.0D-3
  if (glob%po_distribution < 0) glob%po_distribution = 3 ! force_bias
  if (glob%po_maxcycle <= 0) glob%po_maxcycle = 10000
  if (glob%po_init_pop_size < 1) glob%po_init_pop_size = 2*glob%po_pop_size
  if (glob%po_reset < 1) glob%po_reset = 500
  if (glob%po_mutation_rate < 0.0D0) glob%po_mutation_rate = 0.15D0
  if (glob%po_death_rate < 0.0D0) glob%po_death_rate = 0.5D0
  if (glob%po_scalefac <= 0.0D0) glob%po_scalefac = 10.0D0
  if (glob%po_nsave < 0) glob%po_nsave = 10

  if (glob%ntasks <= 0) glob%ntasks = 1

  ! Change any parameters that are incompatible with parallel optimisation
  if (glob%iopt/10==5) then
     !if (glob%icoord/=0) then ! change this in future!!!
     !   write(stdout,'(1x,a)') "Warning: parallel optimisation incompatible &
     !   &with icoord /= 0 currently; switching to Cartesians"
     !   glob%icoord=0
     !end if
     if (glob%iline/=0) then
        write(stdout,'(1x,a)') "Warning: parallel optimisation incompatible &
        &with iline /= 0; setting to zero" 
        glob%iline=0
     end if
  end if

end subroutine dlf_default_set
!!****
