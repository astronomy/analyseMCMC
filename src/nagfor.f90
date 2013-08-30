!> \file nagfor.f90  Provide some redirection/dummy routines to compile the code with NAG Fortran and g95

! 
! LICENCE:
! 
! Copyright (c) 2007-2013  Marc van der Sluys
!  
! This file is part of the AnalyseMCMC package.
!  
! This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
! by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License along with this code (LICENCE).  If not, see 
! <http://www.gnu.org/licenses/>.
! 




!***********************************************************************************************************************************
!> \brief  Dummy version of the system() intrinsic function

function system(str)
   implicit none
   character, intent(in) :: str*(*)
   integer :: system
   character :: dummystr*(99)
   
   dummystr = trim(str)
   system = 0
   
end function system
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Dummy version of the sleep() intrinsic routine

subroutine sleep(nr)
   implicit none
   integer, intent(in) :: nr
   integer :: dummyint
   
   dummyint = nr
   
end subroutine sleep
!***********************************************************************************************************************************

