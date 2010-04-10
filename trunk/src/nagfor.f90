!> \file nagfor.f90
!! \brief Provide some redirection/dummy routines for NAG Fortran
!<



function system(str1)
   implicit none
   character, intent(in) :: str1*(*)
   integer :: system
   character :: dummystr*99
   dummystr = trim(str1)
   system = 0
end function system

subroutine sleep(nr)
   implicit none
   integer, intent(in) :: nr
   integer :: dummyint
   dummyint = nr
end subroutine sleep

