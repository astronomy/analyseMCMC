!> \file analyseMCMC_functions.f90  General routines and functions for analyseMCMC

! 
! LICENCE:
! 
! Copyright (c) 2007-2022  Marc van der Sluys, Vivien Raymond, Ben Farr, Chris Chambers
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
!> \brief  Compute right ascension (in radians) from longitude (radians) and GPS time (seconds)
!!
!! \param lon     Longitude
!! \param GPSsec  GPS time in seconds
!!
!! - Declination == latitude for equatorial coordinates

function lon2ra(lon, GPSsec)
  use SUFR_kinds, only: double
  use SUFR_constants, only: pi2
  
  implicit none
  real(double), intent(in) :: lon,GPSsec
  real(double) :: lon2ra,gmst
  
  lon2ra = mod(lon + gmst(GPSsec) + 10*pi2,pi2)
  
end function lon2ra
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute longitude (in radians) from right ascension (radians) and GPS time (seconds)
!!
!! \param ra      Right Ascension
!! \param GPSsec  GPS time in seconds
!!
!! - Declination == latitude for equatorial coordinates.

function ra2lon(ra, GPSsec)
  use SUFR_kinds, only: double
  use SUFR_constants, only: pi2
  
  implicit none
  real(double), intent(in) :: ra,GPSsec
  real(double) :: ra2lon,gmst
  
  ra2lon = mod(ra - gmst(GPSsec) + 10*pi2,pi2)
  
end function ra2lon
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the 'Greenwich Mean Sidereal Time' (in radians) from GPS time (in seconds)
!!
!! \param GPSsec  GPS time in seconds
!!
!! \see K.R. Lang (1999), p.80sqq.

function gmst(GPSsec)
  use SUFR_kinds, only: double
  use SUFR_constants, only: stdErr
  use SUFR_constants, only: pi2  ! 2*pi
  use SUFR_system, only: warn
  
  implicit none
  real(double), intent(in) :: GPSsec
  real(double) :: gmst,seconds,days,centuries,secCurrentDay
  real(double) :: gps0,leapseconds
  
  gps0 = 630720013.d0  ! GPS time at 1/1/2000 at midnight
  leapseconds = 32.d0  ! At Jan 1st 2000
  if(GPSsec.gt.820108813.d0) leapseconds = leapseconds + 1.d0  ! leapsecond after 1/1/2006
  if(GPSsec.gt.914803214.d0) leapseconds = leapseconds + 1.d0  ! leapsecond after 1/1/2009
  if(GPSsec.lt.630720013.d0) call warn('GMSTs before 01/01/2000 are inaccurate!', stdErr)
  
  ! Time since 1/1/2000 midnight:
  seconds       = (GPSsec - gps0) + (leapseconds - 32.d0)
  days          = floor(seconds/86400.d0) - 0.5d0
  secCurrentDay = mod(seconds, 86400.d0)
  centuries     = days/36525.d0
  gmst = 24110.54841d0 + (centuries*(8640184.812866d0 + centuries*(0.093104d0 + centuries*6.2d-6)))
  gmst = gmst + secCurrentDay * 1.002737909350795d0   ! UTC day is 1.002 * MST day
  gmst = mod(gmst/86400.d0,1.d0)
  gmst = gmst * pi2
  
end function gmst
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief
!!
!! Uses lubksb,ludcmp

subroutine savgol(c, np,nl,nr,ld,m)
  implicit none
  integer, intent(in) :: np,nl,nr,ld,m
  real, intent(out) :: c(np)
  
  integer, parameter :: mmax=6
  
  integer :: imj,ipj,j,k,kk,mm,indx(mmax+1)
  real :: d,fac,sum,a(mmax+1,mmax+1),b(mmax+1)
  
  !if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.mmax.or.nl+nr.lt.m) pause 'bad args in savgol'
  if(np.lt.nl+nr+1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.mmax.or.nl+nr.lt.m) write(0,'(A)')' Bad args in savgol'
  
  do ipj=0,2*m
     sum = 0.
     if(ipj.eq.0)sum = 1.
     do k=1,nr
        sum = sum+float(k)**ipj
     end do
     do k=1,nl
        sum = sum+float(-k)**ipj
     end do
     mm = min(ipj,2*m-ipj)
     do imj=-mm,mm,2
        a(1+(ipj+imj)/2,1+(ipj-imj)/2) = sum
     end do
  end do
  
  call ludcmp(a,m+1,mmax+1,indx,d)
  
  do j=1,m+1
     b(j) = 0.
  end do
  
  b(ld+1) = 1.
  
  call lubksb(a,m+1,mmax+1,indx,b)
  
  do kk=1,np
     c(kk) = 0.
  end do
  
  do k=-nl,nr
     sum = b(1)
     fac = 1.
     do mm=1,m
        fac = fac*real(k)
        sum = sum + b(mm+1)*fac
     end do
     kk = mod(np-k,np) + 1
     c(kk) = sum
  end do
  
end subroutine savgol
!***********************************************************************************************************************************


!***********************************************************************************************************************************
subroutine lubksb(a,n,np,indx,b)
  use SUFR_numerics, only: sne
  
  implicit none
  integer, intent(in) :: n,np,indx(n)
  real, intent(in) :: a(np,np)
  real, intent(out) :: b(n)
  
  integer :: i,ii,j,ll
  real :: sum
  
  ii = 0
  do i=1,n
     ll = indx(i)
     sum = b(ll)
     b(ll) = b(i)
     if (ii.ne.0) then
        do j=ii,i-1
           sum = sum-a(i,j)*b(j)
        end do
     else if (sne(sum,0.)) then
        ii = i
     end if
     b(i) = sum
  end do
  
  do i=n,1,-1
     sum = b(i)
     do j=i+1,n
        sum = sum - a(i,j)*b(j)
     end do
     b(i) = sum/a(i,i)
  end do
  
end subroutine lubksb
!***********************************************************************************************************************************



!***********************************************************************************************************************************
subroutine ludcmp(a,n,np,indx,d)
  use SUFR_numerics, only: seq
  implicit none
  integer, intent(in) :: n,np
  integer, intent(out) :: indx(n)
  real, intent(inout) :: a(np,np)
  real, intent(out) :: d
  
  integer, parameter :: nmax=500
  real, parameter :: tiny=1.0e-20
  integer :: i,imax,j,k
  real :: aamax,dum,sum,vv(nmax)
  
  imax = 0
  
  d = 1.
  do i=1,n
     aamax = 0.
     do j=1,n
        if (abs(a(i,j)).gt.aamax) aamax = abs(a(i,j))
     end do
     !if (seq(aamax,0.)) pause 'singular matrix in ludcmp'
     if(seq(aamax,0.)) write(0,'(A)')' Singular matrix in ludcmp'
     vv(i) = 1./aamax
  end do
  
  do j=1,n
     do i=1,j-1
        sum = a(i,j)
        do k=1,i-1
           sum = sum-a(i,k)*a(k,j)
        end do
        a(i,j) = sum
     end do
     
     aamax = 0.
     
     do i=j,n
        sum = a(i,j)
        do k=1,j-1
           sum = sum-a(i,k)*a(k,j)
        end do
        a(i,j) = sum
        dum = vv(i)*abs(sum)
        if (dum.ge.aamax) then
           imax = i
           aamax = dum
        end if
     end do
     
     if (j.ne.imax) then
        do k=1,n
           dum = a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
        end do
        d = -d
        vv(imax) = vv(j)
     end if
     
     indx(j) = imax
     if(seq(a(j,j),0.)) a(j,j) = tiny
     
     if(j.ne.n) then
        dum = 1./a(j,j)
        do i=j+1,n
           a(i,j) = a(i,j)*dum
        end do
     end if
  end do
  
end subroutine ludcmp
!***********************************************************************************************************************************






!***********************************************************************************************************************************
!> \brief  Returns angle in radians between 0 and 2pi (double precision)
!!
!! \param x  Angle (rad)

function drev2pi(x)
  use SUFR_kinds, only: double
  use SUFR_constants, only: pi
  
  implicit none
  real(double), intent(in) :: x
  real(double) :: drev2pi
  
  drev2pi = x - floor(x/(2*pi))*2*pi
  
end function drev2pi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns periodic value x between 0 and per
!!
!! \param x    Input value
!! \param per  Period of cycle

function revper(x,per)
  implicit none
  real, intent(in) :: x,per
  real :: revper
  
  revper = x - real(floor(x/per)) * per
  
end function revper
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in radians between -pi and pi
!!
!! \param x  Angle (rad)

function revpipi(x)
  use SUFR_constants, only: rpi,rpi2
  
  implicit none
  real, intent(in) :: x
  real :: revpipi
  
  revpipi = x - real(floor(x/rpi2)) * rpi2
  if(revpipi.gt.rpi) revpipi = revpipi - rpi2
  
end function revpipi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in degrees between 0 and 360
!!
!! \param x  Angle (deg)

function rev360(x)
  implicit none
  real, intent(in) :: x
  real :: rev360
  
  rev360 = x - real(floor(x/360.)) * 360.
  
end function rev360
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in degrees between 0 and 180
!!
!! \param x  Angle (deg)

function rev180(x)
  implicit none
  real, intent(in) :: x
  real :: rev180
  
  rev180 = x - real(floor(x/180.)) * 180.
  
end function rev180
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in hours between 0 and 24
!!
!! \param x  Angle (hours)

function rev24(x)
  implicit none
  real, intent(in) :: x
  real :: rev24
  
  rev24 = x - real(floor(x/24.)) * 24.
  
end function rev24
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in radians between 0 and 2pi
!!
!! \param x  Angle (rad)

function rev2pi(x)
  implicit none
  real, intent(in) :: x
  real :: rev2pi,pi
  
  pi = 4*atan(1.)
  rev2pi = x - real(floor(x/(2.0*pi))) * 2.0*pi
  
end function rev2pi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in radians between 0 and pi - double
!!
!! \param x  Angle (rad)

function drevpi(x)
  use SUFR_kinds, only: double
  use SUFR_constants, only: pi
  
  implicit none
  real(double), intent(in) :: x
  real(double) :: drevpi
  
  drevpi = x - real(floor(x/pi)) * pi
  
end function drevpi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Returns angle in radians between 0 and pi - real
!!
!! \param x  Angle (rad)

function rrevpi(x)
  use SUFR_constants, only: rpi
  implicit none
  real, intent(in) :: x
  real :: rrevpi
  
  rrevpi = x - real(floor(x/rpi)) * rpi
  
end function rrevpi
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Print angle as mm:ss.s string, input in hours
!!
!! a1  Angle (hours)

function tms(a1)
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: a1
  
  real(double) :: a,s
  integer :: m
  character :: tms*(8),mm*(2),ss*(4)
  
  a = a1
  m = int((a)*60.d0)
  s = (a-m/60.d0)*3600.d0
  
  write(mm,'(i2.2)') m
  write(ss,'(f4.1)') s
  if(nint(s*10).lt.100) write(ss,'(a1,f3.1)') '0',s
  write(tms,'(a2,a1,a4,a1)') mm,'m',ss,'s'
  
end function tms
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Determine the operating system type: 1-Linux, 2-MacOSX

function getos()
  use SUFR_constants, only: stdErr, homedir
  use SUFR_system, only: warn
  
  implicit none
  integer :: status,system,getos
  character :: ostype*(25),filename*(99)
  
  filename = trim(homedir)//'/.analysemcmc.uname.temp'
  status = system('uname &> '//trim(filename))  ! This should return "Linux" or "Darwin"
  open(unit=16,file=trim(filename), status='old', form='formatted',iostat=status)
  if(status.ne.0) then  ! Something went wrong - guess Linux
     call warn('getOS(): cannot determine OS - guessing Linux...', stdErr)
     getos = 1
     return
  end if
  read(16,'(A)', iostat=status) ostype
  close(16, status='delete')
  
  !write(stdOut,*) ostype
  getos = 1  ! Linux
  if(index(trim(ostype),'Darwin').ne.0) getos = 2  ! MacOSX
  
end function getos
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Get time stamp in seconds since 1970-01-01 00:00:00 UTC, mod countmax

function timestamp()
  use SUFR_kinds, only: double
  
  implicit none
  real(double) :: timestamp
  integer :: count,countmax,countrate
  
  call system_clock(count,countrate,countmax)
  timestamp = dble(mod(count+countmax,countmax))/dble(countrate)
  
end function timestamp
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Transforms longitude l, latitude b and radius r into a vector with length r.  Use r=1 for a unit vector
!!
!! \param l     Longitude (rad)
!! \param b     Latitude (rad)
!! \param r     Radius
!! \retval vec  3D vector with the same units as r

subroutine lbr2vec(l,b,r,vec)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: l,b,r
  real(double), intent(out) :: vec(3)
  
  real(double) :: sinb,cosb
  
  sinb = sin(b)
  cosb = sqrt(1.d0-sinb*sinb)
  vec(1) = cos(l) * cosb         ! 'Greenwich'
  vec(2) = sin(l) * cosb         ! 'Ganges'
  vec(3) = sinb;                 ! 'North Pole'
  vec = vec*r
  
end subroutine lbr2vec
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Compute the length of a 3D cartesian vector
!!
!! \param vec  3D vector

function veclen(vec)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: vec(3)
  real(double) :: veclen
  
  veclen = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
  
end function veclen
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Create a unit vector from a 3D cartesian vector
!!
!! \param vec  3D vector (I/O)

subroutine normvec(vec)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(inout) :: vec(3)
  real(double) :: veclen
  
  vec = vec/veclen(vec)
  
end subroutine normvec
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Convert chirp mass and q to m1 and m2 - double precision
!!
!! \param  mc   Chirp mass (Mo)
!! \param  q    Mass ratio q
!! 
!! \retval m1   M1 (Mo)
!! \retval m2   M2 (Mo)

subroutine mc_q_2_m1_m2(mc,q, m1,m2)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: mc,q
  real(double), intent(out) :: m1,m2
  real(double) :: factor

  factor = mc*(1.d0 + q)**(0.2d0)
  m1 = factor*q**(-0.6d0)
  m2 = factor*q**(0.4d0)
  
end subroutine mc_q_2_m1_m2
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert chirp mass and eta to m1 and m2 - single precision
!! 
!! \param  mcr   Chirp mass (Mo)
!! \param  qr    Mass ratio q
!! 
!! \retval m1r   M1 (Mo)
!! \retval m2r   M2 (Mo)

subroutine mc_q_2_m1_m2r(mcr,qr,m1r,m2r)
  use SUFR_kinds, only: double
  implicit none
  real, intent(in) :: mcr,qr
  real, intent(out) :: m1r,m2r

  real(double) :: mc,q,m1,m2

  mc = dble(mcr)
  q = dble(qr)
  call mc_q_2_m1_m2(mc,q, m1,m2)
  m1r = real(m1)
  m2r = real(m2)

end subroutine mc_q_2_m1_m2r
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert symmetric mass ratio eta to asymmetric mass ratio q (0 - 1)
!!
!! \param  eta   Symmetric mass ratio (0 - 0.25)

function eta2q(eta)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: eta
  real(double) :: eta2q, var

  var = sqrt(1.d0 - 4 * min(max(eta,0.d0),0.25d0))
  eta2q = (1.d0-var) / (1.d0+var)
  
end function eta2q
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert symmetric mass ratio eta to asymmetric mass ratio q (0 - 1)
!!
!! \param  eta   Symmetric mass ratio (0 - 0.25)

function eta2qr(eta)
  implicit none
  real, intent(in) :: eta
  real :: eta2qr, var

  var = sqrt(1. - 4 * min(max(eta,0.),0.25))
  eta2qr = (1.-var) / (1.+var)
  
end function eta2qr
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert asymmetric mass ratio q (0 - 1) to  symmetric mass ratio eta (0 - 0.25)
!!
!! \param  q   Asymmetric mass ratio (0 - 1)

function q2eta(q)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: q
  real(double) :: q2eta

  q2eta = q / (1.d0+q)**2
  
end function q2eta
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert asymmetric mass ratio q (0 - 1) to  symmetric mass ratio eta (0 - 0.25)
!!
!! \param  q   Asymmetric mass ratio (0 - 1)

function q2etar(q)
  implicit none
  real, intent(in) :: q
  real :: q2etar

  q2etar = q / (1.+q)**2
  
end function q2etar
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert chirp mass and eta to m1 and m2 - double precision
!!
!! \param  mc   Chirp mass (Mo)
!! \param  eta  Eta
!!
!! \retval m1   M1 (Mo)
!! \retval m2   M2 (Mo)

subroutine mc_eta_2_m1_m2(mc,eta, m1,m2)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: mc,eta
  real(double), intent(out) :: m1,m2
  real(double) :: dvar,mtot
  
  mtot = mc*eta**(-0.6d0)
  if(eta.le.0.25d0) then
     dvar = sqrt(1.d0-4*eta)
     m1 = mtot/2.d0 * (1.0 + dvar);
     m2 = mtot/2.d0 * (1.0 - dvar);
  else                                 ! Allow 0.25<eta<0.50
     dvar = sqrt(4*eta-1.d0)
     m1 = mtot/2.d0 * (1.0 - dvar);
     m2 = mtot/2.d0 * (1.0 + dvar);
  end if
  
end subroutine mc_eta_2_m1_m2
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert chirp mass and eta to m1 and m2 - single precision
!!
!! \param  mcr   Chirp mass (Mo)
!! \param  etar  Eta
!!
!! \retval m1r   M1 (Mo)
!! \retval m2r   M2 (Mo)

subroutine mc_eta_2_m1_m2r(mcr,etar,m1r,m2r)
  use SUFR_kinds, only: double
  implicit none
  real, intent(in) :: mcr,etar
  real, intent(out) :: m1r,m2r
  
  real(double) :: mc,eta,m1,m2
  
  mc = dble(mcr)
  eta = dble(etar)
  call mc_eta_2_m1_m2(mc,eta,m1,m2)
  m1r = real(m1)
  m2r = real(m2)
  
end subroutine mc_eta_2_m1_m2r
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert M1,M2 to Mchirp
!!
!! \param  m1   M1 (Mo)
!! \param  m2   M2 (Mo)
!!
!! \retval mc   Chirp mass (Mo)

function m1m2_2_mc(m1,m2)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: m1,m2
  real(double) :: m1m2_2_mc, mtot, eta
  
  mtot = m1+m2
  eta  = m1*m2/(mtot*mtot)
  
  m1m2_2_mc = mtot * eta**0.6d0
  
end function m1m2_2_mc
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert M1,M2 to Mchirp - single precision
!!
!! \param  m1   M1 (Mo)
!! \param  m2   M2 (Mo)
!!
!! \retval mc   Chirp mass (Mo)

function m1m2_2_mcr(m1,m2)
  implicit none
  real, intent(in) :: m1,m2
  real :: m1m2_2_mcr, mtot, eta
  
  mtot = m1+m2
  eta  = m1*m2/(mtot*mtot)
  
  m1m2_2_mcr = mtot * eta**0.6
  
end function m1m2_2_mcr
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert M1,M2 to Mchirp, eta  (double precision)
!!
!! \param  m1   M1 (Mo)
!! \param  m2   M2 (Mo)
!!
!! \retval mc   Chirp mass (Mo)
!! \retval eta  Eta

subroutine m1_m2_2_mc_eta(m1,m2, mc,eta)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: m1,m2
  real(double), intent(out) :: mc,eta
  real(double) :: mtot
  
  mtot = m1+m2
  eta = m1*m2/(mtot*mtot)
  mc = mtot*eta**0.6d0
  
end subroutine m1_m2_2_mc_eta
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert M1,M2 to Mchirp, eta  (single precision)
!!
!! \param  m1r   M1 (Mo)
!! \param  m2r   M2 (Mo)
!!
!! \retval mcr   Chirp mass (Mo)
!! \retval etar  Eta

subroutine m1_m2_2_mc_etar(m1r,m2r, mcr,etar)
  use SUFR_kinds, only: double
  implicit none
  real, intent(in) :: m1r,m2r
  real, intent(out) :: mcr,etar
  real(double) :: m1,m2,mc,eta
  
  m1 = dble(m1r)
  m2 = dble(m2r)
  call m1_m2_2_mc_eta(m1,m2,mc,eta)
  mcr = real(mc)
  etar = real(eta)
  
end subroutine m1_m2_2_mc_etar
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert longitude, latitude (rad) to a 3D unit vector
!!
!! \param l     Longitude, in [0,2pi[
!! \param b     Latitude, in [-pi,pi]
!!
!! \retval vec  3D unit vector

subroutine ang2vec(l,b, vec)
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: l,b
  real(double), intent(out) :: vec(3)
  real(double) :: cosb
  
  cosb = cos(b)
  vec(1) = cos(l) * cosb
  vec(2) = sin(l) * cosb
  vec(3) = sin(b)
  
end subroutine  ang2vec
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Convert a 3D vector to longitude, latitude (rad)
!!
!! \param  vec  3D vector
!!
!! \retval l    Longitude, in [0,2pi[
!! \retval b    Latitude, in [-pi,pi]

subroutine vec2ang(vec, l,b)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: vec(3)
  real(double), intent(out) :: l,b
  real(double) :: vec1(3)
  
  vec1 = vec
  call normvec(vec1) !Make sure vec1 is normalised
  l = atan2(vec1(2),vec1(1))
  b = asin(vec1(3))
  
end subroutine  vec2ang
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the dot product of two 3D cartesian vectors
!!
!! \param vec1  3D vector 1
!! \param vec2  3D vector 2

function dotproduct(vec1,vec2)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: vec1(3),vec2(3)
  real(double) :: dotproduct
  
  dotproduct = vec1(1)*vec2(1) + vec1(2)*vec2(2) + vec1(3)*vec2(3)
  
end function dotproduct
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the cross (outer) product of two cartesian vectors
!!
!! \param vec1  3D vector 1
!! \param vec2  3D vector 2
!!
!! \retval crpr  Cross/outer product (vec1 x vec2)

subroutine crossproduct(vec1,vec2, crpr)
  use SUFR_kinds, only: double
  implicit none
  real(double), intent(in) :: vec1(3),vec2(3)
  real(double), intent(out) :: crpr(3)
  
  crpr(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
  crpr(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
  crpr(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
  
end subroutine crossproduct
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief Compute the polarisation angle of a source with position unit vector p and orientation normal vector o
!!
!! \param p  3D position unit vector
!! \param o  3D orientation unit vector
!!
!! \see Apostolatos et al. 1994, Eq.5

function polangle(p,o)  
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: p(3),o(3)
  real(double) :: polangle
  real(double) :: z(3),denom,ocz(3),numer,dotproduct!,datan2
  
  z = (/0.d0,0.d0,1.d0/)                                     ! Vertical unit vector
  denom = dotproduct(o,z) - dotproduct(o,p)*dotproduct(z,p)  ! Denominator
  call crossproduct(o,z,ocz)
  numer = dotproduct(p,ocz)                                  ! Numerator
  
  polangle = atan(denom/(numer+1.d-30))                      ! Take into account the degeneracy in psi, hence no atan2
  
end function polangle
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the position angle of a source with position unit vector p and orientation unit vector o
!!
!! \param p  3D position unit vector
!! \param o  3D orientation unit vector

function posangle(p,o)
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: p(3),o(3)
  real(double) :: posangle
  real(double) :: x1(3),o1(3),z(3),z1(3),dotproduct
  
  call crossproduct(p,o,x1)
  call crossproduct(x1,p,o1) !o1: projection of o in the plane of the sky
  
  z = (/0.d0,0.d0,1.d0/) !Vertical unit vector
  call crossproduct(p,z,x1)
  call crossproduct(x1,p,z1) !z1: projection of z in the plane of the sky
  
  call normvec(o1)
  call normvec(z1)
  
  posangle = acos(dotproduct(o1,z1))
  
end function posangle
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Compute the inclination and polarisation angle for a source with position (pl,pb) and orientation (ol,ob)
!! 
!! \param  pl   Position: longtitude, in [0,2pi[  (rad)
!! \param  pb   Position: latitude, in [0,pi]  (rad)
!! \param  ol   Orientation: longtitude, in [0,2pi[  (rad)
!! \param  ob   Orientation: latitude, in [0,pi]  (rad)  
!!
!! \retval i    Inclination angle (rad)
!! \retval psi  Polarisation angle (rad)
!! 
!! \note 
!! - all variables are angles (no cos, sin)
!! - pb,ob used to be in ([-pi/2,pi/2]) now [0,pi], conf John V. & Christian R.

subroutine compute_incli_polang(pl,pb,ol,ob, i,psi) 
  use SUFR_kinds, only: double
  
  implicit none
  
  real(double), intent(in) :: pl,pb,ol,ob
  real(double), intent(out) :: i,psi
  real(double) :: p(3),o(3),dotproduct,polangle,drevpi
  
  call ang2vec(pl,pb,p)       ! Position unit vector
  call ang2vec(ol,ob,o)       ! Orientation unit vector
  
  ! Definition 1:
  ! Compute inclination angle: <0: points towards us, >0 points away from us:
  !i = pi2 - acos(dotproduct(p,o))
  ! Compute polarisation angle [-pi/2,pi/2]:
  !psi = polangle(p,o)
  
  ! Definition 2:
  ! Compute inclination angle: 0: points exactly away from us, 180 points exactly towards us, 90: in the plane of the sky:
  i = acos(dotproduct(p,o))  
  ! Compute polarisation angle [0,pi]:
  psi = drevpi(polangle(p,o))
  
end subroutine compute_incli_polang
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the inclination and polarisation angle for a source with position (pl,pb) and orientation (ol,ob) - single prec.
!! 
!! \param  plr   Position: longtitude, in [0,2pi[  (rad)
!! \param  pbr   Position: latitude, in [0,pi]  (rad)
!! \param  olr   Orientation: longtitude, in [0,2pi[  (rad)
!! \param  obr   Orientation: latitude, in [0,pi]  (rad)  
!! \retval ir    Inclination angle (rad)
!! \retval psir  Polarisation angle (rad)
!! 
!! \note
!! - all variables are angles (no cos, sin)
!! - single-precision wrapper for compute_incli_polang()

subroutine compute_incli_polangr(plr,pbr,olr,obr, ir,psir)
  use SUFR_kinds, only: double
  
  implicit none
  real, intent(in) :: plr,pbr,olr,obr
  real, intent(out) :: ir,psir
  real(double) :: pl,pb,ol,ob,i,psi
  
  pl = dble(plr)
  pb = dble(pbr)
  ol = dble(olr)
  ob = dble(obr)
  call compute_incli_polang(pl,pb,ol,ob, i,psi)
  ir = real(i)
  psir = real(psi)
  
end subroutine compute_incli_polangr
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the inclination and position angle for a source with position (pl,pb) and orientation (ol,ob)
!!
!! \param  pl   Position: longtitude, in [0,2pi[  (rad)
!! \param  pb   Position: latitude, in  (rad)
!! \param  ol   Orientation: longtitude, in [0,2pi[  (rad)
!! \param  ob   Orientation: latitude, in  (rad)  
!! \retval i    Inclination angle (rad)
!! \retval pa   Polarisation angle (rad)
!!
!! \note  Position angle, not polarisation angle!

subroutine compute_incli_posang(pl,pb,ol,ob, i,pa) 
  use SUFR_kinds, only: double
  
  implicit none
  real(double), intent(in) :: pl,pb,ol,ob
  real(double), intent(out) :: i,pa
  real(double) :: p(3),o(3),dotproduct,posangle
  
  call ang2vec(pl,pb,p)       ! Position unit vector
  call ang2vec(ol,ob,o)       ! Orientation unit vector
  
  ! Definition 1:
  ! Compute inclination angle: <0: points towards us, >0 points away from us
  !i = pi2 - acos(dotproduct(p,o))
  ! Compute position angle:
  !pa = drevpi(posangle(p,o))
  
  ! Definition 2:
  ! Compute inclination angle: 0: points exactly away from us, 180 points exactly towards us, 90: in the plane of the sky:
  i = acos(dotproduct(p,o))  
  ! Compute position angle:
  pa = posangle(p,o)
  
end subroutine compute_incli_posang
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Determine the sky position pointed at by the vector that connects two detectors
!!
!! \param d1  ID Detector 1
!! \param d2  ID Detector 2
!! \param jd  Julian day
!!
!! \todo  Finish and use

subroutine detectorvector(d1,d2,jd)
  use SUFR_kinds, only: double
  implicit none
  integer, intent(in) :: d1,d2
  real(double), intent(inout) :: jd  ! should become (in) once used
  real(double) :: detcoords(3,2),vec1(3),vec2(3),dvec(3),l,b
  
  jd = 0 !get rid of warnings
  detcoords(1,:) = (/-119.41,46.45/)  !H1; l,b
  detcoords(2,:) = (/-90.77,30.56/)   !L1
  detcoords(3,:) = (/10.50,43.63/)    !V
  
  call ang2vec(detcoords(d1,1),detcoords(d1,2),vec1)
  call ang2vec(detcoords(d2,1),detcoords(d2,2),vec2)
  
  dvec = vec2 - vec1
  
  call vec2ang(dvec,l,b)  !Searched point is in zenith/nadir for an observer on this location on the globe
  
end subroutine detectorvector
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Find files in the current working directory
!!
!! \param match    Search string to match
!! \param nff      Maximum number of files to return
!! \param all      All files?  0-select manually from list, 1-always return all files in list
!!
!! \retval fnames  Array that contains the files found; make sure it has the same length as the array in the calling programme
!! \retval nf      The actual number of files returned in fnames ( = min(number found, nff))

subroutine findFiles(match,nff,all, fnames,nf)  
  use SUFR_constants, only: stdErr, homedir
  use SUFR_system, only: quit_program_error, warn
  
  implicit none
  character, intent(in) :: match*(*)
  integer, intent(in) :: nff,all
  character, intent(out) :: fnames(nff)*(99)
  integer, intent(out) :: nf
  
  integer :: i,j,k,fnum,status,system,io
  character :: names(nff)*(99),tempfile*(99)
  
  if(len_trim(homedir).ge.99) call quit_program_error('FindFiles: variable homedir not defined (forgot to call setconstants?)',1)
  
  tempfile = trim(homedir)//'/.findFile.tmp'
  ! Shell command to list all the files with the search string and pipe them to a temporary file:
  status = system('ls '//trim(match)//' 1> '//trim(tempfile)//' 2> /dev/null') 
  status = status  ! Remove 'set but never used' warning
  
  do i=1,nff
     names(i)=''
  end do
  
  k=0
  open(10,file=trim(tempfile), status='old', form='formatted',iostat=io)  ! Read the temp file and delete it when closing
  if(io.ne.0) then
     call warn('findFiles(): cannot list files in current directory...', stdErr)
     fnames = ''
     nf = 0
     return
  end if
  rewind(10)
  do i=1,nff
     read(10,'(A99)',end=100) names(i)
     k=k+1
  end do
100 continue
  close(10, status='delete')
  fnames(1) = names(1)
  nf = 1
  j = 0
  
  if(k.gt.1) then
     if(all.eq.0) then  ! Select files manually
        write(6,'(A)')'  Files found:'  ! Don't use stdOut here!
        do i=1,k
           write(6,'(I5,A3,A)')i,':  ',trim(names(i))
        end do
        write(6,*)
        write(6,'(A,I3)')'  Enter the number of the file you want to select: 1 -',k
        write(6,'(A,I3,A1)')'    (max',nff,')'
        write(6,'(A)')'      or:   0 - to select all files in the list'
        write(6,'(A)')'           -1 - when done'
        do j=1,nff
           read*,fnum
           if(fnum.lt.0) then
              nf = j-1
              return
           end if
           if(fnum.eq.0) then
              nf = min(k,nff)
              fnames(1:nf) = names(1:nf)
              return
           end if !if(fnum.eq.0)
           fnames(j) = names(fnum)
           nf = j
        end do !j 
     else  ! Select all files (all=1)
        nf = min(k,nff)
        fnames(1:nf) = names(1:nf)
        return
     end if
  end if
  
  if(k.eq.0) then
     fnames(1)='                                                                                                   '
     !write(stdErr,'(A)')'  No file found in this directory'
     nf = 0
  end if
  
end subroutine findFiles
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Compute the size needed for PGPlot to get the desired bitmap size in pixels
!!
!! \param bmpXSz   Desired x-size in pixels
!! \param bmpYSz   Desired y-size in pixels
!! \param scFac    Scale factor; produce a larger bitmap, the shrink to get smoother graphics
!!
!! \retval bmpsz   Size of the bitmap (x)
!! \retval bmprat  Aspect ration of the bitmap

subroutine compBitmapSize(bmpXSz,bmpYSz, scFac, bmpsz,bmprat)
  use aM_constants, only: use_PLplot
  implicit none
  integer, intent(in) :: bmpXSz,bmpYSz
  real, intent(in) :: scFac
  real, intent(out) :: bmpsz,bmprat
  
  if(use_PLplot) then
     bmpsz = real(bmpXSz)/300.               ! PLplot (300 dpi?), no resizing afterwards
  else
     bmpsz = real(bmpXSz-1)/85. * scFac      ! PGPlot: Make png larger, so that convert interpolates and makes the plot smoother
  end if
  bmprat = real(bmpYSz-1)/real(bmpXSz-1)
  
end subroutine compBitmapSize
!***********************************************************************************************************************************





!***********************************************************************************************************************************
!> \brief  Print a single output line to specify when and were AnalyseMCMC was run
!!
!! \param op  Output unit

subroutine print_rundata(op)
  use SUFR_constants, only: workdir,hostname,username,currenttimezonestr,currenttimestr,currentdatestr
  use analysemcmc_settings, only: htmlOutput
  use general_data, only: infiles,nchains0
  
  implicit none
  integer, intent(in) :: op
  
  if(htmlOutput.ge.1) then
     write(op,'(A)', advance="no")'  Analysed'
  else
     write(op,'(A)', advance="no")'  Analysing'
  end if
  if(nchains0.eq.1) then
     write(op,'(A)', advance="no")' 1 chain from'
  else
     write(op,'(I3,A)', advance="no") nchains0,' chains from'
  end if
  if(index(infiles(1),'SPINspiral.output').ne.0 .or. index(infiles(1),'mcmc.output').ne.0) then
     write(op,'(A)', advance="no")' SPINspiral'
  else
     write(op,'(A)', advance="no")' LALInference'
  end if
  write(op,'(A)', advance='no')',  in '//trim(username)//'@'//trim(hostname)//':'//trim(workdir)//'/'
  write(op,'(A)')',  on '//trim(currentdatestr)//', '//trim(currenttimestr)//' ('//trim(currenttimezonestr)//').'
  
end subroutine print_rundata
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!>  Setup a index.html file for output

subroutine create_html_index_file()
  use SUFR_constants, only: stdOut
  use SUFR_system, only: quit_program_error
  
  use analysemcmc_settings, only: prRunInfo,prChainInfo,prConv,prStat,prIval,prCorr
  use analysemcmc_settings, only: plot, plLogL,plChain,plParL,plJump,plAcorr,plRhat, plPDF1D,plPDF2D
  use aM_constants, only: stdOutFile
  use mcmcrun_data, only: t0
  
  implicit none
  integer :: io
  
  ! Some initial output was already written to (tmp)stdOut.  Close that file, reopen it and overwrite it:
  close(stdOut)
  open(unit=stdOut,action='write',form='formatted',status='replace',file=trim(stdOutFile), iostat=io)
  if(io.ne.0) call quit_program_error('Error reopening output file '//trim(stdOutFile),1)
  
  write(stdOut,'(A)') '<!DOCTYPE html>'
  write(stdOut,'(A)') '<html>'
  write(stdOut,'(A)') '  <head>'
  write(stdOut,'(4x,A,I11,A)') '<title>AnalyseMCMC: event',nint(t0),'</title>'
  write(stdOut,'(A)') '  </head>'
  write(stdOut,'(A)') '  <body>'
  
  write(stdOut,'(4x,A)') '<a name="top"></a>'
  write(stdOut,'(4x,A)') '<font size="2">'
  
  write(stdOut,'(6x,A)') '<b>Jump to:</b> &nbsp;'
  if(prRunInfo.gt.0)              write(stdOut,'(6x,A)') '<a href="#runinfo">Run info</a>'
  if(prChainInfo.gt.0)            write(stdOut,'(6x,A)') '&nbsp; &ndash; &nbsp; <a href="#chaininfo">Chain info</a>'
  if(prConv.ge.2)                 write(stdOut,'(6x,A)') '&nbsp; &ndash; &nbsp; <a href="#mixing">Mixing</a>'
  if(prStat.gt.0)                 write(stdOut,'(6x,A)') '&nbsp; &ndash; &nbsp; <a href="#stats">Main statistics</a>'
  if(prIval.eq.1.or.prIval.eq.3)  write(stdOut,'(6x,A)') '&nbsp; &ndash; &nbsp; <a href="#prob">Probability intervals</a>'
  if(prIval.ge.2)                 write(stdOut,'(6x,A)') '&nbsp; &ndash; &nbsp; <a href="#statsprob">Stats &amp; prob.ivals </a>'
  if(prCorr.gt.0)                 write(stdOut,'(6x,A)') '&nbsp; &ndash; &nbsp; <a href="#corr">Correlations</a>'
  if(plot.gt.0)                   write(stdOut,'(6x,A)') '&nbsp; &ndash; &nbsp; <a href="#plots">Plots</a>:'
  
  if(plLogL.gt.0)   write(stdOut,'(6x,A)') '&ndash; <a href="#postchains">Post.chains</a>'
  if(plChain.gt.0)  write(stdOut,'(6x,A)') '&ndash; <a href="#parchains">Par.chains</a>'
  if(plParL.gt.0)   write(stdOut,'(6x,A)') '&ndash; <a href="#par-l">Par.-L</a>'
  if(plJump.gt.0)   write(stdOut,'(6x,A)') '&ndash; <a href="#jumps">Jumps</a>'
  if(plAcorr.gt.0)  write(stdOut,'(6x,A)') '&ndash; <a href="#acorrs">A.corrs</a>'
  if(plRhat.gt.0)   write(stdOut,'(6x,A)') '&ndash; <a href="#rhat">R-hat</a>'
  if(plPDF1D.gt.0)  write(stdOut,'(6x,A)') '&ndash; <a href="#1dpdfs">1D PDFs</a>'
  if(plPDF2D.gt.0)  write(stdOut,'(6x,A)') '&ndash; <a href="#2dpdfs">2D PDFs</a>'
  
  write(stdOut,'(4x,A)') '</font>'
  write(stdOut,'(4x,A)') '<br>'
  
  write(stdOut,'(4x,A,I11,A)') '<h1><a href="http://analysemcmc.sourceforge.net/" target="_blank" title="AnalyseMCMC homepage">'// &
       'AnalyseMCMC</a>: event',nint(t0),'</h1>'
  
  write(stdOut,'(A)') '<pre>'
  
end subroutine create_html_index_file
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!>  Setup the 2dpdf.html file for output

subroutine create_html_2dpdf_file(op)
  use SUFR_system, only: quit_program_error
  use mcmcrun_data, only: t0
  
  implicit none
  integer, intent(in) :: op
  integer :: io
  character :: filename*(99)
  
  write(filename,'(A)') 'html/2dpdfs.html'
  open(unit=op, action='write', form='formatted', status='replace', file=trim(filename), iostat=io)
  if(io.ne.0) call quit_program_error('Error opening output file '//trim(filename),1)
  
  write(op,'(A)') '<!DOCTYPE html>'
  write(op,'(A)') '<html>'
  write(op,'(A)') '  <head>'
  write(op,'(4x,A,I11,A)') '<title>AnalyseMCMC: event',nint(t0),' - 2D PDF matrix</title>'
  write(op,'(A)') '  </head>'
  write(op,'(A)') '  <body>'
  
  write(op,'(4x,A)') '<a name="top"></a>'
  write(op,'(4x,A)') '<font size="2">'
  write(op,'(6x,A)') '<a href="index.html#2dpdfs" title="Go back to the main page">Main page</a>'
  write(op,'(4x,A)') '</font>'
  write(op,'(4x,A)') '<br>'
  
  write(op,'(4x,A,I11,A)') '<h1><a href="http://analysemcmc.sourceforge.net/" target="_blank" title="AnalyseMCMC homepage">'// &
       'AnalyseMCMC</a>: event',nint(t0),'</h1>'
  write(op,'(4x,A)') '<h2>2D PDF matrix</h2>'
  
end subroutine create_html_2dpdf_file
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Report an undefined parameter that the code tries to use
!!
!! \param parName  Name of the undefined parameter
!! \param parID    ID of the undefined parameter
!! \param routine  Name of the caller routine

subroutine report_undefined_parameter(parName, parID, routine)
  use SUFR_constants, only: stdOut,stdErr
  use analysemcmc_settings, only: prProgress
  implicit none
  character, intent(in) :: parName*(*), routine*(*)
  integer, intent(in) :: parID
  
  select case(prProgress)
  case(0)
  case(1)
     write(stdOut,'(A,I0,A)', advance='no') ' !!par. "'//trim(parName)//'" (',parID,') undefined!! '
  case default  ! 2,3: verbose and debug output
     write(stdErr,'(/,A,I0,A)')'  * Warning:  '//trim(routine)//'():  parameter "'//trim(parName)// &
          '" with parID ',parID,' is not defined, check plPars() in the input file.  Skipping...'
  end select
  
end subroutine report_undefined_parameter
!***********************************************************************************************************************************
