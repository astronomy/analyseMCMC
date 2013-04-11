!> \file analyseMCMC_plot.f90  Plotting routines for analyseMCMC

! 
! LICENCE:
! 
! Copyright 2007-2013 Marc van der Sluys
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
!> \brief Initialise PGPlot for use with AnalyseMCMC
!!
!! \param colour   Use colour:  0:no, 1-yes
!! \param file     File type: 0-screen, 1-png, 2-eps, 3-pdf
!! \param whiteBG  Use a white background:  0:no, 1-yes

subroutine pginitl(colour,file,whiteBG)
  use aM_constants, only: use_PLplot
  implicit none
  integer, intent(in) :: colour,file,whiteBG
  integer :: i
  
  if(whiteBG.ge.1) then
     call pgscr(0,1.,1.,1.)                ! Background colour always white (also on screen, bitmap)
     call pgscr(1,0.,0.,0.)                ! Default foreground colour always black
     
     if(file.le.1) then                    ! png: create white background
        call pgsvp(-100.,100.,-100.,100.)
        call pgswin(0.,1.,0.,1.)
        call pgsci(0)
        call pgrect(-1.,2.,-1.,2.)
        
        if(use_PLplot) then
           call pgsvp(0.08,0.92,0.12,0.94) ! Default viewport size (?)
        else
           call pgsvp(0.08,0.95,0.06,0.87) ! Default viewport size (?)
        end if
     end if
     call pgsci(1)
  end if
  
  if(colour.eq.0) then
     do i=0,99
        call pgscr(i,0.,0.,0.)
     end do
     call pgscr(0,1.,1.,1.)      ! White
     call pgscr(14,0.3,0.3,0.3)  ! Dark grey
     call pgscr(15,0.6,0.6,0.6)  ! Light grey
  else
     call pgscr(2,1.,0.1,0.1)    ! Make default red lighter
     call pgscr(3,0.,0.5,0.)     ! Make default green darker
     call pgscr(4,0.,0.,0.8)     ! Make default blue darker
     call pgscr(7,0.9,0.9,0.)    ! Make default yellow darker
     call pgscr(10,0.5,0.3,0.)   ! 10: brown
     call pgscr(11,0.6,0.,0.)    ! 11: dark red
  end if
  
end subroutine pginitl
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief Set the colour to ci, but use a darker shade if the background is black or a lighter shade if it is white
!! 
!! \param ci       Colour index
!! \param file     Output file type
!! \param whiteBG  Have a white background (1) or nor (0)

subroutine pgscidark(ci,file,whiteBG)  
  implicit none
  integer, intent(in) :: ci,file,whiteBG
  integer :: ci1
  real :: r,g,b,weight
  
  call pgqcr(ci,r,g,b)
  
  ci1 = 99
  ! Use half the RGB value to create a darker shade:
  call pgscr(ci1,r*0.5,g*0.5,b*0.5)
  
  ! Use the weighted mean of the RGB value and 1. to create a lighter shade:
  weight = 3.
  if(file.ge.2.or.whiteBG.ge.1) call pgscr(ci1,(r+weight)/(weight+1.),(g+weight)/(weight+1.),(b+weight)/(weight+1.)) 
  call pgsci(ci1)
  
end subroutine pgscidark
!***********************************************************************************************************************************



