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
module dlf_parameter_module
  integer, parameter :: rk = kind(1.d0) ! read kind
  integer, parameter :: ik = kind(1)    ! integer kind (for memory tracking)
end module dlf_parameter_module

module dlf_stat
  implicit none
  type stat_type
    integer    :: sene=0      ! number of total energy evaluations (on the 
                              ! current processor)
    integer    :: ccycle=0    ! number of cycles/steps in the current run
    integer    :: caccepted=0 ! number of accepted steps in the current run
  end type stat_type
  type(stat_type),save :: stat
end module dlf_stat

subroutine dlf_stat_reset
  use dlf_stat, only: stat
  implicit none
  stat%ccycle=0
  stat%caccepted=0 
end subroutine dlf_stat_reset
