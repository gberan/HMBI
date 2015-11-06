!!****h* utilities/linalg
!!
!! FUNCTION
!! Linear algebra utilities
!!
!! These subroutines provide wrappers for lapack and blas routines
!! so that the latter are only called once in the code. This 
!! should facilitate porting.
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
!!****f* linalg/dlf_matrix_multiply
!!
!! FUNCTION
!!
!! C= alpha * AxB + beta * C
!!
!! SYNOPSIS
subroutine dlf_matrix_multiply(M,N,K,alpha,A,B,beta,C)
  use dlf_parameter_module, only: rk
  implicit none
  integer  ,intent(in)    :: M,N,K
  real(rk) ,intent(in)    :: alpha,beta
  real(rk) ,intent(in)    :: A(M,K) ! K is the common index
  real(rk) ,intent(in)    :: B(K,N)
  real(rk) ,intent(inout) :: C(M,N)
!! SOURCE
! **********************************************************************
  CALL dgemm('N','N',M,N,K,alpha, A , M, B, K, beta, C, M)
end subroutine dlf_matrix_multiply
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* linalg/dlf_matrix_invert
!!
!! FUNCTION
!!
!! A <- A^-1
!!
!! SYNOPSIS
subroutine dlf_matrix_invert(N,tdet,a,det)
  use dlf_parameter_module, only: rk
  use dlfhdlc_matrixlib, only: array_invert 
  implicit none
  integer  ,intent(in)    :: N
  logical  ,intent(in)    :: tdet
  real(rk) ,intent(inout) :: A(N,N)
  real(rk) ,intent(out)   :: det
!! SOURCE
  integer :: idum
! **********************************************************************
  idum = array_invert(A,det,tdet,N)
end subroutine dlf_matrix_invert
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* linalg/dlf_matrix_diagonalise
!!
!! FUNCTION
!!
!! Calculate eigenvalues and eigenvectors of a real, symmetric matrix A
!!
!! The eigenvector to the eigenvalue evals(ivar) is evecs(:,ivar)
!!
!! SYNOPSIS
subroutine dlf_matrix_diagonalise(N,a,evals,evecs)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlfhdlc_matrixlib, only: array_diagonalise
  implicit none
  integer  ,intent(in)    :: N
  real(rk) ,intent(in)    :: A(N,N)
  real(rk) ,intent(out)   :: evals(N)
  real(rk) ,intent(out)   :: evecs(N,N)
  integer :: idum
! **********************************************************************
  idum = array_diagonalise(a,evecs,evals,n,n,n,.true.)
end subroutine dlf_matrix_diagonalise
!!****
