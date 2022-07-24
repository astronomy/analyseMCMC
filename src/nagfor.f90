!> \file nagfor.f90  Provide some redirection/dummy routines to compile the code with NAG Fortran and g95

! 
! LICENCE:
! 
! Copyright (c) 2007-2022  Marc van der Sluys
!  
! This file is part of the AnalyseMCMC package, see http://analysemcmc.sf.net and https://github.com/Astronomy/AnalyseMCMC.
!  
! This is free software: you can redistribute it and/or modify it under the terms of the European Union
! Public Licence 1.2 (EUPL 1.2).
! 
! This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the EU Public License for more details.
! 
! You should have received a copy of the European Union Public License along with this code.  If not, see
! <https://www.eupl.eu/1.2/en/>.
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

