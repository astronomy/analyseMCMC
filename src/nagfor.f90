!> \file nagfor.f90  Provide some redirection/dummy routines to compile the code with NAG Fortran and g95




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

