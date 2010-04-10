!> \file nagfor.f90
!! \brief Provide some redirection/dummy routines for NAG Fortran
!<



subroutine getenv(str1,str2)
   implicit none
   character, intent(in) :: str1*(*)
   character, intent(out) :: str2*(*)
   call get_environment_variable(str1,str2)
end subroutine getenv

function iargc()
   implicit none
   integer :: iargc
   iargc = command_argument_count()
end function iargc

subroutine getarg(nr,str)
   implicit none
   integer, intent(in) :: nr
   character, intent(out) :: str*(*)
   call get_command_argument(nr,str)
end subroutine getarg




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

