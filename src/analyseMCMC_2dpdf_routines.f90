!> \file  analyseMCMC_2dpdf_routines.f90  Routines and functions to help produce 2D marginalised PDFs

! 
! LICENCE:
! 
! Copyright 2007-2012 Marc van der Sluys
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
!> \brief  Bin data and 'normalise' 2D PDF
!!

subroutine bin_and_normalise_2D_data(ic,p1,p2, xmin,xmax, ymin,ymax, z,tr, sky_position,binary_orientation)
  use SUFR_constants, only: stdOut, pi,rpi
  use SUFR_statistics, only: bin_data_2d
  use SUFR_text, only: replace_substring
  
  use analysemcmc_settings, only: normPDF2D, maxChs, Nbin2Dx,Nbin2Dy, prProgress, Nival,prIval,ivals
  use general_data, only: selDat, n, maxIter, startval, pgUnits
  use stats_data, only: probArea,probAreas, injectionranges2d
  use mcmcrun_data, only: parID
  
  implicit none
  integer, intent(in) :: ic, p1,p2
  real, intent(inout) :: xmin,xmax, ymin,ymax, z(Nbin2Dx+1,Nbin2Dy+1), tr(6)
  logical, intent(in) :: sky_position,binary_orientation
  
  integer :: i, injectionrange2d
  real :: xx(maxChs*maxIter), yy(maxChs*maxIter), zz(maxChs*maxIter)
  character :: areaunit*(19)
  
  
  xx(1:n(ic)) = selDat(ic,p1,1:n(ic))  ! Parameter 1
  yy(1:n(ic)) = selDat(ic,p2,1:n(ic))  ! Parameter 2
  zz(1:n(ic)) = selDat(ic,1,1:n(ic))   ! Likelihood
  
  
  if(normPDF2D.le.2.or.normPDF2D.eq.4) then
     
     
     ! Bin data in 2D:
     call bin_data_2d( xx(1:n(ic)), yy(1:n(ic)), 0, Nbin2Dx,Nbin2Dy, xmin,xmax, ymin,ymax, z, tr )
     
     
     !Test:
     !call check_binned_data(Nbin2Dx,Nbin2Dy,z)
     
     !do Nbin2Dx = 10,200,10
     !   Nbin2Dy = Nbin2Dx
     !   xmin1 = xmin
     !   xmax1 = xmax
     !   ymin1 = ymin
     !   ymax1 = ymax
     !   
     !   !Bin data:  compute bin number rather than find it, ~10x faster:
     !   call bin_data_2d(xx(1:n(ic)),yy(1:n(ic)),0,Nbin2Dx,Nbin2Dy,xmin1,xmax1,ymin1,ymax1,z,tr)
     !   
     !   !Test!
     !   call check_binned_data(Nbin2Dx,Nbin2Dy,z)
     !   
     !end do
     !stop
     
     
     
     if(normPDF2D.eq.1) z = max(0.,log10(z + 1.e-30))
     if(normPDF2D.eq.2) z = max(0.,sqrt(z + 1.e-30))
     
     if(normPDF2D.eq.4) then
        
        ! Get 2D probability ranges; identify to which range each bin belongs:
        if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  identifying 2D ranges...'
        call identify_2d_ranges(p1,p2,Nival,Nbin2Dx+1,Nbin2Dy+1,z,tr)
        
        ! Compute 2D probability areas; sum the areas of all bins:
        if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  computing 2D areas...'
        call calc_2d_areas(p1,p2,Nival,Nbin2Dx+1,Nbin2Dy+1,z,tr,probArea)
        injectionranges2d(p1,p2) = injectionrange2d(z,Nbin2Dx+1,Nbin2Dy+1,startval(1,p1,1),startval(1,p2,1),tr)
        
        do i=1,Nival
           if(prIval.ge.1.and.prProgress.ge.2 .and. (sky_position .or. binary_orientation)) then  
              ! For sky position and orientation only:
              if(i.eq.1) write(stdOut,'(/,1x,A10,A13,3A23)') 'Nr.','Ival frac.','Area (sq.deg) ', &
                   'Circ. area rad. (deg) ','Fraction of sky '
              write(stdOut,'(I10,F13.2,3(2x,F21.5))') i,ivals(i),probArea(i),sqrt(probArea(i)/pi)*2, &
                   probArea(i)*(pi/180.)**2/(4*pi)  ! 4pi*(180/pi)^2 = 41252.961 sq. degrees in a sphere
           else
              areaunit = trim(pgUnits(parID(p1)))//' '//trim(pgUnits(parID(p2)))
              if(trim(pgUnits(parID(p1))) .eq. trim(pgUnits(parID(p2)))) areaunit = trim(pgUnits(parID(p1)))//'^2'  ! mm->m^2
              call replace_substring(areaunit, '\(2218)', 'deg')    ! degrees
              call replace_substring(areaunit, '\d\(2281)\u', 'o')  ! Mo
              call replace_substring(areaunit, '\dh\u', 'hr')       ! hr
              areaunit = ' '//trim(areaunit)  ! Add space between value and unit
              
              if(i.eq.1) write(stdOut,'(/,1x,A10,A13,A23)') 'Nr.','Ival frac.','Area'
              write(stdOut,'(I10,F13.2,2x,1p,G21.3,1x,A)') i,ivals(i),probArea(i),trim(areaunit)
           end if
           probAreas(p1,p2,i,1) = probArea(i)*(rpi/180.)**2/(4*rpi)  ! Fraction of the sky
           probAreas(p1,p2,i,2) = sqrt(probArea(i)/rpi)*2            ! Equivalent diameter
           probAreas(p1,p2,i,3) = probArea(i)                        ! Square degrees
        end do
     end if
  end if  ! if(normPDF2D.le.2.or.normPDF2D.eq.4)
  
  
  ! 'Bin' the data by weighing by likelihood value:
  if(normPDF2D.eq.3) then
     if(prProgress.ge.3) write(stdOut,'(A)',advance="no")'  binning 2D data...'
     ! Measure amount of likelihood in each bin:
     call bin_data_2d_a( n(ic), xx(1:n(ic)), yy(1:n(ic)), zz(1:n(ic)), Nbin2Dx,Nbin2Dy, xmin,xmax,ymin,ymax, z, tr )
  end if
  
  
  ! Normalise the binned data:
  z = z/(maxval(z)+1.e-30)
  
end subroutine bin_and_normalise_2D_data
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief 'Bin data' in 2 dimensions  -  Measure the amount of likelihood in each bin
!! 
!! \param ndat   Input data: array size
!! \param xdat   Input data: x values
!! \param ydat   Input data: y values
!! \param zdat   Input data: likelihood values
!!
!! \param nxbin  Desired number of bins in the x direction
!! \param nybin  Desired number of bins in the y direction
!!
!! \param xmin1  Lower limit for the binning range in the x direction - autodetermine if xmin1=xmax1
!! \param xmax1  Upper limit for the binning range in the x direction - autodetermine if xmin1=xmax1
!! \param ymin1  Lower limit for the binning range in the y direction - autodetermine if ymin1=ymax1
!! \param ymax1  Upper limit for the binning range in the y direction - autodetermine if ymin1=ymax1
!!
!! \retval zz    'Binned' data set z(nxbin,nybin) (real)
!! \retval tr    Transformation elements for pgplot tr(6) (real)
!!
!! \todo  Should z be replaced?

subroutine bin_data_2d_a(ndat, xdat,ydat, zdat, nxbin,nybin, xmin1,xmax1,ymin1,ymax1, zz, tr)
  implicit none
  integer, intent(in) :: ndat, nxbin,nybin
  integer :: i,bx,by
  real, intent(in) :: xdat(ndat),ydat(ndat),zdat(ndat)
  real, intent(inout) :: xmin1,xmax1,ymin1,ymax1
  real, intent(out) :: zz(nxbin+1,nybin+1), tr(6)
  
  real :: xbin(nxbin+1),ybin(nybin+1), zztot, xmin,xmax,ymin,ymax, dx,dy ,zmin
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  zmin = minval(zdat)
  
  if(abs((xmin-xmax)/(xmax+1.e-30)).lt.1.e-20) then  ! Autodetermine
     xmin = minval(xdat(1:ndat))
     xmax = maxval(xdat(1:ndat))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs((ymin-ymax)/(ymax+1.e-30)).lt.1.e-20) then  ! Autodetermine
     ymin = minval(ydat(1:ndat))
     ymax = maxval(ydat(1:ndat))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  do bx=1,nxbin+1
     !xbin(bx) = xmin + (real(bx)-0.5)*dx  ! x is the centre of the bin
     xbin(bx) = xmin + real(bx-1)*dx           ! x is the left of the bin
  end do
  do by=1,nybin+1
     !ybin(by) = ymin + (real(by)-0.5)*dy  ! y is the centre of the bin
     ybin(by) = ymin + real(by-1)*dy           ! y is the left of the bin
  end do
  
  !write(stdOut,'(50F5.2)'),xdat(1:50)
  !write(stdOut,'(50F5.2)'),ydat(1:50)
  !write(stdOut,'(20F8.5)'),xbin
  !write(stdOut,'(20F8.5)'),ybin
  
  zz = 0.
  zztot = 0.
  do bx=1,nxbin
     do by=1,nybin
        zz(bx,by) = 0.
        do i=1,ndat
           !if(xdat(i).ge.xbin(bx).and.xdat(i).lt.xbin(bx+1) .and. ydat(i).ge.ybin(by).and.ydat(i).lt.ybin(by+1))  &
           !zz(bx,by) = zz(bx,by) + 1.
           if(xdat(i).ge.xbin(bx).and.xdat(i).lt.xbin(bx+1) .and. ydat(i).ge.ybin(by).and.ydat(i).lt.ybin(by+1))  &
                zz(bx,by) = zz(bx,by) + exp(zdat(i) - zmin)
           !write(stdOut,'(2I4,8F10.5)')bx,by,xdat(i),xbin(bx),xbin(bx+1),ydat(i),ybin(by),ybin(by+1),zz(bx,by),zdat(i)
        end do
        zztot = zztot + zz(bx,by) 
        !write(stdOut,'(2I4,5x,4F6.3,5x,10I8)')bx,by,xbin(bx),xbin(bx+1),ybin(by),ybin(by+1),nint(zz(bx,by))
     end do
     !write(stdOut,'(I4,5x,2F6.3,5x,10I8)')bx,xbin(bx),xbin(bx+1),nint(zz(bx,1:nybin))
  end do
  
  if(abs((xmin1-xmax1)/(xmax1+1.e-30)).lt.1.e-20) then  ! Autodetermine
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs((ymin1-ymax1)/(ymax1+1.e-30)).lt.1.e-20) then  ! Autodetermine
     ymin1 = ymin
     ymax1 = ymax
  end if
  
  ! Determine transformation elements for pgplot (pggray, pgcont):
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
end subroutine bin_data_2d_a
!***********************************************************************************************************************************




!***********************************************************************************************************************************
!> \brief  Get the 2d probability intervals; z lies between 1 (in 100% range) and ni (in lowest-% range, e.g. 90%)
!!
!! \param p1  Parameter ID 1
!! \param p2  Parameter ID 2
!! \param ni  Number of probability intervals
!! \param nx  Number of bins in the x direction
!! \param ny  Number of bins in the y direction
!! \param z   Binned data (nx,ny) -> probability-interval data
!! \param tr  Transformation elements used by PGPlot

subroutine identify_2d_ranges(p1,p2,ni,nx,ny,z,tr)
  use SUFR_constants, only: stdOut
  use SUFR_constants, only: rd2r
  use SUFR_sorting, only: sorted_index_list
  
  use analysemcmc_settings, only: changeVar,ivals,prProgress
  use mcmcrun_data, only: parID
  
  implicit none
  integer, intent(in) :: p1,p2,ni,nx,ny
  real, intent(inout) :: z(nx,ny)
  real, intent(in) :: tr(6)
  
  integer :: nn,indx(nx*ny),i,b,ib,full(ni),ix,iy, zxy,dix,diy,same,same0,iso, nonzero
  real :: x1(nx*ny),x2(nx*ny),tot,np,y
  
  
  ! Weigh number of points in each bin by bin size for position/orientation plots:
  do iy = 1,ny
     if(changeVar.ge.1) then
        if((parID(p1).eq.31.and.parID(p2).eq.32) .or. (parID(p1).eq.52.and.parID(p2).eq.51)) then  
           ! Then: RA-Dec or (phi/theta_Jo)/(psi/i) plot, convert lon -> lon * 15 * cos(lat)
           y = tr(4) + tr(6)*real(iy)
           if(parID(p1).eq.31) then
              if(abs(y).le.90.) then
                 z(1:nx,iy) = z(1:nx,iy)/(cos(y*rd2r)+1.e-30)
              else  !This can happen when the PDF lies close to the pole
                 z(1:nx,iy) = 0.
              end if
           else if(parID(p1).eq.52) then
              if(y.ge.0..and.y.lt.180.) then
                 z(1:nx,iy) = z(1:nx,iy)/(abs(sin(y*rd2r))+1.e-30)
              else  !This can happen when the PDF lies close to the pole
                 z(1:nx,iy) = 0.
                 !write(stdErr,'(//,A,//)')'  *** identify_2d_ranges:  sin(y)<0. '// &
                 !' Please check whether the if(y.ge.0..and.y.lt.180.) statement works properly ***'
              end if
           end if
        end if
     end if  ! if(changeVar.ge.1)
  end do  ! iy
  
  
  nn = nx*ny
  x1 = reshape(z,(/nn/))                 ! x1 is an 1D array with the same data as the 2D array z
  call sorted_index_list(dble(-x1(1:nn)), indx(1:nn))  ! -x1: sort the 1D array to descending value
  
  np = sum(z)
  tot = 0.
  full = 0
  do b=1,nn  ! Loop over bins in 1D array
     ib = indx(b)
     x2(ib) = 0.
     if(x1(ib).eq.0.) cycle
     tot = tot + x1(ib)
     do i=ni,1,-1  ! Loop over intervals
        if(tot.le.np*ivals(i)) then
           x2(ib) = real(ni-i+1)  ! e.g. x2(b) = ni if within 68%, ni-1 if within 95%, etc, and 1 if within 99.7%
        else
           if(prProgress.ge.3.and.full(i).eq.0) then  ! Report the number of points in the lastly selected bin
              if(i.eq.1) write(stdOut,'(A)',advance="no") 'Last bin:'
              !write(stdOut,'(F6.3,I5)',advance="no") ivals(i),nint(x1(ib))
              write(stdOut,'(I5)',advance="no") nint(x1(ib))
              full(i) = 1
           end if
        end if
        !write(stdOut,'(2I4, F6.2, 3F20.5)')b,i, ivals(i), np,tot,np*ivals(i)
     end do
  end do
  
  z = reshape(x2, (/nx,ny/))  ! z lies between 1 and ni
  
  ! Print the matrix:
  !do b=ny,1,-1
  !   print*,nint(z(:,b))
  !end do
  
  ! Count isolated pixels (surrounded by pixels that are all of a different value than it):
  if(prProgress.ge.3) then
     write(*,*)
     do same0=1,8
        iso = 0
        nonzero = 0
        do iy=ny-1,2,-1
           do ix=2,nx-1
              
              zxy = nint(z(ix,iy))
              if(zxy.eq.0) cycle  ! Consider coloured pixels only
              
              same = 0
              do diy = -1,1
                 do dix = -1,1
                    if(dix.eq.0.and.diy.eq.0) cycle
                    !if( zxy.eq.nint(z(ix+dix,iy+diy)) ) same = same+1
                    if( nint(z(ix+dix,iy+diy)) .ne. 0 ) same = same+1  ! white pixels only
                 end do
              end do
              
              if(same.lt.same0) iso = iso+1
              
              if(zxy.gt.0) nonzero = nonzero + 1
           end do
        end do
        
        write(*,'(A,I2,A,I6,A1,F7.2,A1)') '  Isolated pixels (fewer than',same0,' identical pixels around it):',iso,',', &
             real(iso)/real(nn)*100,'%'
     end do
     write(*,'(A,I6,A1,F7.2,A1)') '  Non-zero pixels :',nonzero,',', &
          real(nonzero)/real(nn)*100,'%'
     
  end if  ! if(prProgress.ge.3)
  
end subroutine identify_2d_ranges
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Compute 2D probability areas
!!
!! \param p1  Parameter ID 1
!! \param p2  Parameter ID 2
!! \param ni  Number of probability intervals
!! \param nx  Number of bins in the x direction
!! \param ny  Number of bins in the y direction
!! \param z   Binned data (nx,ny) -> probability-interval data
!! \param tr  Transformation elements used by PGPlot
!! \retval area  Probability areas

subroutine calc_2d_areas(p1,p2,ni,nx,ny,z,tr,area)
  use SUFR_constants, only: rd2r
  use analysemcmc_settings, only: changeVar
  use mcmcrun_data, only: parID
  
  implicit none
  integer, intent(in) :: p1,p2,ni,nx,ny
  real, intent(in) :: z(nx,ny),tr(6)
  real, intent(out) :: area(ni)
  
  integer :: ix,iy,i,i1,iv
  real :: y,dx,dy
  
  
  area = 0.
  
  do ix = 1,nx
     do iy = 1,ny
        dx = tr(2)
        dy = tr(6)
        if(changeVar.ge.1) then
           if((parID(p1).eq.31.and.parID(p2).eq.32) .or. (parID(p1).eq.52.and.parID(p2).eq.51)) then  
              !Then: RA-Dec or (phi/theta_Jo)/(psi/i) plot, convert lon -> lon * 15 * cos(lat)
              y = tr(4) + tr(6)*real(iy)
              if(parID(p1).eq.31) then
                 dx = dx*cos(y*rd2r)
              else if(parID(p1).eq.52) then
                 dx = dx*abs(sin(y*rd2r))  !Necessary for i-psi plot?
              end if
              if(parID(p1).eq.31) dx = dx*15
           end if
        end if
        iv = nint(z(ix,iy))
        do i=1,ni
           if(iv.ge.i) then
              i1 = ni-i+1
              area(i1) = area(i1) + dx*dy
           end if
        end do !i
        
     end do !iy
  end do !ix
  
end subroutine calc_2d_areas
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Get the smallest probability area in which the injection values lie
!!
!! \param z   Binned data (nx,ny) -> probability-interval data
!! \param nx  Number of bins in the x direction
!! \param ny  Number of bins in the y direction
!! \param injectionx  Injection value of parameter 1
!! \param injectiony  Injection value of parameter 2
!! \param tr  Transformation elements used by PGPlot

function injectionrange2d(z,nx,ny,injectionx,injectiony,tr)
  implicit none
  integer, intent(in) :: nx,ny
  real, intent(in) :: injectionx,injectiony,z(nx,ny),tr(6)
  integer :: ix,iy,injectionrange2d
  
  !x = tr(1) + tr(2)*ix + tr(3)*iy
  !y = tr(4) + tr(5)*ix + tr(6)*iy
  ix = floor((injectionx - tr(1))/tr(2))
  iy = floor((injectiony - tr(4))/tr(6))
  if(ix.lt.1.or.ix.gt.nx.or.iy.lt.1.or.iy.gt.ny) then
     injectionrange2d = 0
  else
     injectionrange2d = nint(z(ix,iy))
  end if
  
end function injectionrange2d
!***********************************************************************************************************************************








!***********************************************************************************************************************************
!> \brief Find the fraction of non-zero bins for which more than frac*8 of their neighbours are smaller
!!
!! \param nxbin  Desired number of bins in the x direction
!! \param nybin  Desired number of bins in the y direction
!! \param z     Binned data set z(nxbin,nybin) (real)

subroutine check_binned_data(nxbin,nybin,z)
  
  implicit none
  integer, intent(in) :: nxbin,nybin
  real, intent(in) :: z(nxbin+1,nybin+1)
  
  integer :: n,bx,by,bbx,bby,count,nsb
  real :: frac
  
  frac = 1.0
  
  n = 0
  count = 0
  do bx = 2,nxbin-1
     do by = 2,nybin-1
        
        if(z(bx,by).lt.1.e-30) cycle  ! Don't do empty bins
        
        nsb = 0
        !print*,bx,by,z(bx,by)
        do bbx = -1,1
           do bby = -1,1
              if(bbx.eq.0.and.bby.eq.0) cycle
              if(z(bx+bbx,by+bby).lt.z(bx,by)) nsb = nsb + 1  ! If the neigbouring bin is smaller, count it
              !if(z(bx+bbx,by+bby).lt.1.e-30) nsb = nsb + 1   ! If the neigbouring bin is empty, count it
           end do
        end do
        
        n = n+1
        if(nsb.ge.nint(frac*real(8))) count = count+1
        
     end do
  end do
  
  write(6,'(//,A,F10.2,2I6,3I8,F10.4,//)') 'Check_binned_data:',frac,nxbin,nybin,(nxbin-2)*(nybin-2),n,count,real(count)/real(n)
  
end subroutine check_binned_data
!***********************************************************************************************************************************




!***********************************************************************************************************************************

!> \brief  Prepare binning for a cute sky map in 2D PDF
!!
!! \param xmin  Minimum value of x (RA) range
!! \param xmax  Maximum value of x (RA) range
!! \param ymin  Minimum value of y (Dec) range
!! \param ymax  Maximum value of y (Dec) range

subroutine prepare_skymap_binning(xmin,xmax, ymin,ymax)
  use SUFR_constants, only: stdOut
  use analysemcmc_settings, only: prProgress
  
  implicit none
  real, intent(inout) :: xmin,xmax, ymin,ymax
  real :: rat, avg, dx,dy
  
  rat = 0.5  ! scrRat
  
  dx = xmax - xmin
  dy = ymax - ymin
  
  if(abs(dx)*15.lt.dy/rat) then  ! Expand x
     dx = dy/(15*rat)
     avg = (xmin+xmax)*0.5
     xmin = avg - 0.5*dx
     xmax = avg + 0.5*dx
     if(prProgress.ge.3) write(stdOut,'(2(A,F6.1),A)',advance='no')'  Changing RA binning range to ',xmin,' - ',xmax,' h.'
  end if
  
  if(abs(dx)*15.gt.dy/rat) then  ! Expand y
     dy = abs(dx)*rat*15
     avg = (ymin+ymax)*0.5
     ymin = avg - 0.5*dy
     ymax = avg + 0.5*dy
     if(prProgress.ge.3) write(stdOut,'(2(A,F6.1),A)',advance='no')'  Changing Dec. binning range to ',ymin,' - ',ymax,' deg.'
  end if
  
end subroutine prepare_skymap_binning
!***********************************************************************************************************************************




!***********************************************************************************************************************************
!> \brief  Define the colours for the 2D probability areas
!!
!! \retval clr1  First colour of colour-index range for pgimag
!! \retval clr2  Second colour of colour-index range for pgimag

subroutine set_2D_probability_colours(clr1,clr2)
  use analysemcmc_settings, only: colour, Nival, file
  
  implicit none
  integer, intent(out) :: clr1,clr2
  
  
  if(colour.eq.0) then
     call pgscr(30,1.,1.,1.)  ! BG colour
     if(Nival.eq.2) then
        call pgscr(31,0.5,0.5,0.5)  ! Grey
        call pgscr(32,0.,0.,0.)     ! Black
     end if
     if(Nival.eq.3) then
        call pgscr(31,0.7,0.7,0.7)  ! Light grey
        call pgscr(32,0.4,0.4,0.4)  ! Dark grey
        call pgscr(33,0.0,0.0,0.0)  ! Black
     end if
     if(Nival.eq.4) then
        call pgscr(31,0.75,0.75,0.75)  ! Light grey
        call pgscr(32,0.50,0.50,0.50)  ! 
        call pgscr(33,0.25,0.25,0.25)  ! Dark grey
        call pgscr(34,0.00,0.00,0.00)  ! Black
     end if
     if(Nival.eq.5) then
        call pgscr(31,0.8,0.8,0.8)  ! Light grey
        call pgscr(32,0.6,0.6,0.6)  ! 
        call pgscr(33,0.4,0.4,0.4)  ! 
        call pgscr(34,0.2,0.2,0.2)  ! Dark grey
        call pgscr(35,0.0,0.0,0.0)  ! Black
     end if
  end if
  
  if(colour.ge.1) then
     call pgscr(30,1.,1.,1.)  ! BG colour
     if(Nival.eq.2) then
        call pgscr(31,1.,1.,0.)  ! Yellow
        if(file.ge.2) call pgscr(31,0.8,0.7,0.)  ! Dark yellow
        call pgscr(32,1.,0.,0.)  ! Red
     end if
     if(Nival.eq.3) then
        call pgscr(31,0.,0.,1.)  ! Blue
        call pgscr(32,1.,1.,0.)  ! Yellow
        if(file.ge.2) call pgscr(32,0.8,0.7,0.)  ! Dark yellow
        call pgscr(33,1.,0.,0.)  ! Red
     end if
     if(Nival.eq.4) then
        call pgscr(31,0.,0.,1.)  ! Blue
        call pgscr(32,0.,1.,0.)  ! Green
        call pgscr(33,1.,1.,0.)  ! Yellow
        if(file.ge.2) call pgscr(33,0.8,0.7,0.)  ! Dark yellow
        call pgscr(34,1.,0.,0.)  ! Red
     end if
     if(Nival.eq.5) then
        call pgscr(31,0.,0.,1.)  ! Blue
        call pgscr(32,0.,1.,0.)  ! Green
        call pgscr(33,1.,1.,0.)  ! Yellow
        if(file.ge.2) call pgscr(33,0.8,0.7,0.)  ! Dark yellow
        call pgscr(34,1.,0.5,0.)  ! Orange
        call pgscr(35,1.,0.,0.)  ! Red
     end if
  end if
  
  clr1 = 30
  clr2 = 30+Nival
  call pgscir(clr1,clr2)  ! Set colour-index range for pgimag

end subroutine set_2D_probability_colours
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Plot contours in a 2D PDF
!!
!! \param z            2D binned data
!! \param tr           Transformation elements used by PGPlot
!! \param project_map  Use map projection?
!! \param lw           Default line width

subroutine plot_2D_contours(z, tr, project_map, lw)
  use analysemcmc_settings, only: normPDF2D, plotSky, Nival, Nbin2Dx,Nbin2Dy
  
  implicit none
  real, intent(in) :: z(Nbin2Dx+1,Nbin2Dy+1), tr(6)
  logical, intent(in) :: project_map
  integer, intent(in) :: lw
  
  integer :: i, Ncont
  real :: cont(11)
  
  
  if(normPDF2D.lt.4) then
     Ncont = 11
     do i=1,Ncont
        cont(i) = 0.01 + 2*real(i-1)/real(Ncont-1)
        if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) cont(i) = 1.-cont(i)
     end do
     Ncont = min(4,Ncont)  ! Only use the first 4
  else if(normPDF2D.eq.4) then
     Ncont = Nival
     do i=1,Ncont
        cont(i) = max(1. - real(i-1)/real(Ncont-1),0.001)
        !if(project_map) cont(i) = 1.-cont(i)
     end do
  end if
  
  call pgsls(1)
  if((.not.project_map .or. plotSky.ne.1.or.plotSky.ne.3) .and. normPDF2D.ne.4) then  ! First in bg colour
     call pgslw(2*lw)
     call pgsci(0)
     call pgcont(z,Nbin2Dx+1,Nbin2Dy+1,1,Nbin2Dx+1,1,Nbin2Dy+1,cont(1:Ncont),Ncont,tr)
  end if
  
  call pgslw(lw)
  call pgsci(1)
  if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) call pgsci(7)
  
  call pgcont(z,Nbin2Dx+1,Nbin2Dy+1,1,Nbin2Dx+1,1,Nbin2Dy+1,cont(1:Ncont),Ncont,tr)  ! Plot contours
  
end subroutine plot_2D_contours
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief Plot injection value, median, ranges, etc. in 2D PDF
!!
!! \param ic           Chain ID
!! \param p1           ID of parameter 1
!! \param p2           ID of parameter 2
!!
!! \param xmin         Lower limit of horizontal plot range
!! \param xmax         Upper limit of horizontal plot range
!! \param ymin         Lower limit of vertical plot range
!! \param ymax         Upper limit of vertical plot range
!! \param dx           Width of horizontal plot range
!! \param dy           Width of vertical plot range
!!
!! \param sch          Default (character) scaling
!! \param lw           Default line width
!! \param project_map  Use map projection?


subroutine plot_values_in_2D_PDF(ic, p1,p2, xmin,xmax, ymin,ymax, dx,dy, sch,lw, project_map)
  use analysemcmc_settings, only: plotSky, plLmax,plInject,plRange,plMedian, mergeChains, ivals, normPDF2D, map_projection
  use general_data, only: allDat,startval, c0, ranges,stats, icloglmax,iloglmax, wrap,shifts,shIvals,raCentre
  use mcmcrun_data, only: parID
  
  implicit none
  integer, intent(in) :: ic, p1,p2, lw
  real, intent(in) :: xmin,xmax, ymin,ymax, dx,dy, sch
  logical, intent(in) :: project_map
  
  real :: plx,ply
  character :: delta*(19)
  
  
  if(.not.project_map.or.plotSky.eq.1) then
     
     call pgsci(1)
     
     
     ! Plot max likelihood in 2D PDF:
     if(plLmax.ge.1) then
        call pgsci(1)
        call pgsls(5)
        
        plx = allDat(icloglmax,p1,iloglmax)
        if(wrap(ic,p1).ne.0) plx = mod(plx + shifts(ic,p1), shIvals(ic,p1)) - shifts(ic,p1)
        call pgline(2,(/plx,plx/),(/ymin,ymax/))  ! Max logL
        
        ply = allDat(icloglmax,p2,iloglmax)
        if(wrap(ic,p2).ne.0) ply = mod(ply + shifts(ic,p2), shIvals(ic,p2)) - shifts(ic,p2)
        call pgline(2,(/xmin,xmax/),(/ply,ply/))  ! Max logL
        
        call pgpoint(1,plx,ply,18)
     end if
     
     
     if(project_map .and. (plotSky.eq.1.or.plotSky.eq.3)) call pgsci(0)
     call pgsls(2)
     
     
     ! Plot injection value in 2D PDF:
     if((plInject.eq.1.or.plInject.eq.3).and.(.not.project_map) .or. &
          ((plInject.eq.2.or.plInject.eq.4) .and.  &
          (parID(p1).eq.61.and.parID(p2).eq.62 .or. parID(p1).eq.63.and.parID(p2).eq.64 .or. &
          parID(p2).eq.61.and.parID(p2).eq.67  .or. parID(p2).eq.61.and.parID(p2).eq.68)) ) then
        
        ! CHECK The units of the injection values haven't changed (e.g. from rad to deg) for ic>1 
        ! (but they have for the starting values, why?)
        if(mergeChains.ne.1.or.ic.le.1) then
           
           call pgsls(3)  ! Dash-dotted line for injection value
           call pgsci(1)
           
           ! x:
           plx = startval(ic,p1,1)
           if(wrap(ic,p1).ne.0) plx = mod(plx + shifts(ic,p1), shIvals(ic,p1)) - shifts(ic,p1)
           call pgline(2,(/plx,plx/),(/ymin,ymax/))  ! Injection value
           
           ! y:
           ply = startval(ic,p2,1)
           if(wrap(ic,p2).ne.0) ply = mod(ply + shifts(ic,p2), shIvals(ic,p2)) - shifts(ic,p2)
           call pgline(2,(/xmin,xmax/),(/ply,ply/))  ! Injection value
           
           call pgpoint(1,plx,ply,18)
        end if
     end if  !If plotting injection values in 2D plot
     
     call pgsci(1)
     call pgsls(4)
     
     
     ! Plot starting values in 2D PDF:
     !call pgline(2,(/startval(ic,p1,2),startval(ic,p1,2)/),(/ymin,ymax/))
     !call pgline(2,(/xmin,xmax/),(/startval(ic,p2,2),startval(ic,p2,2)/))
     
     call pgsci(2)
     
     
     ! Plot probability ranges in 2D PDF:
     if(plRange.eq.2.or.plRange.eq.3.or.plRange.eq.5.or.plRange.eq.6) then
        write(delta,'(A,I3.3,A)')'\(2030)\d',nint(ivals(c0)*100),'%\u'
        if(nint(ivals(c0)*100).lt.100) write(delta,'(A,I2.2,A)')'\(2030)\d',nint(ivals(c0)*100),'%\u'
        
        call pgsls(1)
        call pgsch(sch*0.6)
        call pgsah(1,45.,0.1)
        
        call pgarro( ranges(ic,c0,p1,3), ymin+dy*0.017*sch, ranges(ic,c0,p1,1), ymin+dy*0.017*sch)
        call pgarro( ranges(ic,c0,p1,3), ymin+dy*0.017*sch, ranges(ic,c0,p1,2), ymin+dy*0.017*sch)
        call pgptxt( ranges(ic,c0,p1,3), ymin+dy*0.033*sch, 0., 0.5, trim(delta) )
        
        call pgarro( xmin+dx*0.023*sch, ranges(ic,c0,p2,3), xmin+dx*0.023*sch, ranges(ic,c0,p2,1) )
        call pgarro( xmin+dx*0.023*sch, ranges(ic,c0,p2,3), xmin+dx*0.023*sch, ranges(ic,c0,p2,2) )
        call pgptxt( xmin+dx*0.01*sch, ranges(ic,c0,p2,3), 90., 0.5, trim(delta) )
     end if
     
     call pgsch(sch)
     call pgsls(2)
     
     
     ! Plot medians in 2D PDF:
     if(plMedian.eq.2.or.plMedian.eq.3.or.plMedian.eq.5.or.plMedian.eq.6) then
        call pgline(2,(/stats(ic,p1,1),stats(ic,p1,1)/),(/ymin,ymax/))
        call pgline(2,(/xmin,xmax/),(/stats(ic,p2,1),stats(ic,p2,1)/))
        call pgpoint(1,stats(ic,p1,1),stats(ic,p2,1),18)
     end if
     
     call pgsls(1)
     
  end if  ! if(.not.project_map.or.plotSky.eq.1)
  
  
  
  ! Plot big symbol at injection position in sky map:
  if(project_map .and. (plInject.eq.1.or.plInject.eq.3)) then
     call pgsch(sch*1.5)               ! Use 1.5 for plsym=8, 2 for plsym=18
     call pgslw(lw*2)
     call pgsci(9)
     if(normPDF2D.eq.4) call pgsci(1)  ! Black
     
     plx = startval(ic,p1,1)
     ply = startval(ic,p2,1)
     
     if(plotSky.eq.2.or.plotSky.eq.4) call project_skymap(plx,ply,raCentre,map_projection)
     call pgpoint(1,plx,ply,8)
     
     call pgsch(sch)
     call pgslw(lw)
     call pgsci(1)
  end if
  
  
end subroutine plot_values_in_2D_PDF
!***********************************************************************************************************************************



