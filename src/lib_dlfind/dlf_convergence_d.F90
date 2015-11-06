! **********************************************************************
! **                     Test for convergence                         **
! **********************************************************************
!!****h* DL-FIND/convergence
!!
!! NAME
!! convergence
!!
!! FUNCTION
!! Test if the optimisation is converged
!!
!!
!! DATA
!! $Date: 2010-05-07 19:30:43 $
!! $Rev: 307 $
!! $Author: gberan $
!! $URL: http://ccpforge.cse.rl.ac.uk/svn/dl-find/branches/release_chemsh3.3/dlf_convergence.f90 $
!! $Id: dlf_convergence.f90,v 1.1 2010-05-07 19:30:43 gberan Exp $
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
module dlf_convergence
  use dlf_parameter_module, only: rk
  real(rk),save     :: tole    ! Energy tolerance
  real(rk),save     :: tolg    ! Maximum gradient tolerance
  real(rk),save     :: tolrmsg ! RMS gradient tolerance
  real(rk),save     :: tols    ! Maximum step tolerance
  real(rk),save     :: tolrmss ! RMS step tolerance
  ! for external providing:
  logical ,save     :: texternal=.false. ! energy, grad, and step set externally
  character(30),save:: message="" 
  real(rk),save     :: vale    ! Energy value
  real(rk),save     :: valg    ! Maximum gradient value
  real(rk),save     :: valrmsg ! RMS gradient value
  real(rk),save     :: vals    ! Maximum step value
  real(rk),save     :: valrmss ! RMS step value
end module dlf_convergence

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* convergence/convergence_test
!!
!! FUNCTION
!!
!! test for convergence and print convergence information
!! gradient and step information is obtained from convergence_set_info
!!
!! If convergence_set_info has not been called previously, it is called
!! from here with the glob%... arrays as input
!!
!! INPUTS
!!
!! glob%tolerance, glob%tolerance_e
!!
!! SYNOPSIS
subroutine convergence_test(icycle,tene,tconv)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stdout,printl
  use dlf_convergence
  implicit none
  integer,intent(in)  :: icycle ! just for printing
  logical,intent(in)  :: tene   ! store current energy as old energy?
  logical,intent(out) :: tconv  ! converged?
  !
  logical  :: le,lrmss,ls,lrmsg,lg
  real(rk) :: svar
! **********************************************************************

    ! This is initialisation in reality
  if(glob%tolerance<0.D0) call dlf_fail("Convergence tolerance < 0")
  tolg=glob%tolerance

  tole=glob%tolerance_e
  tolrmsg= tolg / 1.5D0
  tols=    tolg * 4.D0
  tolrmss= tolg * 8.D0/3.D0

  ! get values from glob%... arrays
  if(.not.texternal) call convergence_set_info("", &
      glob%nivar,glob%energy,glob%igradient,glob%step)
  texternal=.false.

  if(printl>0) write(stdout,'(3a,i5)') &
      "Testing convergence ",trim(message)," in cycle",icycle


  if(glob%toldenergy_conv) then

    ! energy convergence
    svar=abs(vale-glob%oldenergy_conv)
    le=(svar < tole)
    if(printl>0) then
      if(le) then
        write(stdout,1000) "Energy",svar,tole,"yes"
      else 
        write(stdout,1000) "Energy",svar,tole,"no"
      end if
    end if

    ! Maximum step convergence
    ls=(vals < tols)
    if(printl>0) then
      if(ls) then
        write(stdout,1000) "Max step",vals,tols,"yes"
      else
        write(stdout,1000) "Max step",vals,tols,"no"
      end if
    end if
    
    ! RMS step convergence
    lrmss=(valrmss < tolrmss)
    if(printl>0) then
      if(lrmss) then
        write(stdout,1000) "RMS step",valrmss,tolrmss,"yes"
      else
        write(stdout,1000) "RMS step",valrmss,tolrmss,"no"
      end if
    end if
    
  else
    if(tene) then
      le=.true.
      ls=.true.
      lrmss=.true.
      glob%toldenergy_conv=.true.
    end if
  end if
  if(tene) glob%oldenergy_conv=vale

  ! Maximum gradient convergence
  lg=(valg < tolg)
  if(printl>0) then
    if(lg) then
      write(stdout,1000) "Max grad",valg,tolg,"yes"
    else
      write(stdout,1000) "Max grad",valg,tolg,"no"
    end if
  end if

  ! RMS gradient convergence
  lrmsg=(valrmsg < tolrmsg)
  if(printl>0) then
    if(lrmsg) then
      write(stdout,1000) "RMS grad",valrmsg,tolrmsg,"yes"
    else
      write(stdout,1000) "RMS grad",valrmsg,tolrmsg,"no"
    end if
  end if

  ! overall convergence
  !tconv=(le.and.lg.and.lrmsg)
  tconv=(le.and.lg.and.lrmsg.and.ls.and.lrmss)
  if(printl>0.and.tconv) write(stdout,"('Converged!')")

  ! send convergence information to task module
  call dlf_task_set_l("CONVERGED",tconv)

! formats
1000 format (a10,2x,es10.4," Target: ",es10.4," converged? ",a4)
end subroutine convergence_test
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* convergence/convergence_get
!!
!! FUNCTION
!!
!! Provides the currently set tolerance to other modules
!!
!! SYNOPSIS
subroutine convergence_get(name,val)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob
  use dlf_convergence
  implicit none
  character(*),intent(in)  :: name
  real(rk)    ,intent(out) :: val
! **********************************************************************
  if(name=="TOLG") then
    val=glob%tolerance
  else
    call dlf_fail("Wrong name in conv_get")
  end if
end subroutine convergence_get
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* convergence/convergence_set_info
!!
!! FUNCTION
!!
!! used to provide complete convergence info externally (from other
!! modules), rather than using the glob%... arrays
!!
!! SYNOPSIS
subroutine convergence_set_info(msg,nvar,energy,gradient,step)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_convergence
  implicit none
  character(*),intent(in)  :: msg
  integer     ,intent(in)  :: nvar
  real(rk)    ,intent(in)  :: energy
  real(rk)    ,intent(in)  :: gradient(nvar)
  real(rk)    ,intent(in)  :: step(nvar)
! **********************************************************************
  message=msg
  texternal=.true.
  vale=energy
  valg=maxval(abs(gradient(:)))
  valrmsg=sqrt(sum(gradient(:)**2)/dble(nvar))
  vals=maxval(abs(step(:)))
  valrmss=sqrt(sum(step(:)**2)/dble(nvar))
end subroutine convergence_set_info
!!****
    
