!> \file  analyseMCMC_2dpdf_routines.f90  Routines and functions to help produce 2D marginalised PDFs





!***********************************************************************************************************************************
!> \brief  Bin data 2D by counting the number of points in each bin - old routine
!! 
!! \param n      Number of data points
!! \param x      Data to be binned (n points), parameter 1
!! \param y      Data to be binned (n points), parameter 2
!! \param norm   Normalise histogram (1) or not (0)
!! \param nxbin  Desired number of bins, parameter 1
!! \param nybin  Desired number of bins, parameter 2
!!
!! \param xmin1  Minimum value of the binning range, parameter 1.  Set xmin=xmax to auto-determine (I/O)
!! \param xmax1  Maximum value of the binning range, parameter 1.  Set xmin=xmax to auto-determine (I/O)
!! \param ymin1  Minimum value of the binning range, parameter 2.  Set ymin=ymax to auto-determine (I/O)
!! \param ymax1  Maximum value of the binning range, parameter 2.  Set ymin=ymax to auto-determine (I/O)
!!
!! \retval z     Binned data, height of the bins
!! \retval tr    Transformation elements for pgplot; determines ranges of the output matrix z

subroutine bindata2dold(n,x,y,norm,nxbin,nybin, xmin1,xmax1,ymin1,ymax1, z,tr)
  
  implicit none
  integer, intent(in) :: n, nxbin,nybin, norm
  real, intent(in) :: x(n),y(n)
  real, intent(inout) :: xmin1,xmax1,ymin1,ymax1
  real, intent(out) :: z(nxbin+1,nybin+1),tr(6)
  
  integer :: i,bx,by
  real :: xbin(nxbin+1),ybin(nybin+1)
  real :: xmin,xmax,ymin,ymax,dx,dy
  
  !write(stdOut,'(A4,5I8)')'n:',norm,nxbin,nybin
  !write(stdOut,'(A4,2F8.3)')'x:',xmin1,xmax1
  !write(stdOut,'(A4,2F8.3)')'y:',ymin1,ymax1
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  
  if(abs((xmin-xmax)/(xmax+1.e-30)).lt.1.e-20) then  ! Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs((ymin-ymax)/(ymax+1.e-30)).lt.1.e-20) then  ! Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  do bx=1,nxbin+1
     !xbin(bx) = xmin + (real(bx)-0.5)*dx  ! x is the centre of the bin
     xbin(bx) = xmin + (bx-1)*dx           ! x is the left of the bin
  end do
  do by=1,nybin+1
     !ybin(by) = ymin + (real(by)-0.5)*dy  ! y is the centre of the bin
     ybin(by) = ymin + (by-1)*dy           ! y is the left of the bin
  end do
  
  !write(stdOut,'(50F5.2)'),x(1:50)
  !write(stdOut,'(50F5.2)'),y(1:50)
  !write(stdOut,'(20F8.5)'),xbin
  !write(stdOut,'(20F8.5)'),ybin
  
  z = 0.
  !ztot = 0.
  do i=1,n
     bxl: do bx=1,nxbin
        do by=1,nybin
           if(x(i).ge.xbin(bx)) then
              if(x(i).lt.xbin(bx+1)) then
                 if(y(i).ge.ybin(by)) then
                    if(y(i).lt.ybin(by+1)) then
                       z(bx,by) = z(bx,by) + 1.
                       exit bxl    ! Exit bx loop; if point i fits this bin, don't try other bins. Speeds things up ~2x
                    end if
                 end if
              end if
           end if
           
        end do !by
     end do bxl !bx
  end do !i
  !if(norm.eq.1) z = z/(ztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs((xmin1-xmax1)/(xmax1+1.e-30)).lt.1.e-20) then  ! Autodetermine
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs((ymin1-ymax1)/(ymax1+1.e-30)).lt.1.e-20) then  ! Autodetermine
     ymin1 = ymin
     ymax1 = ymax
  end if
  
  ! Determine transformation elements for pgplot (pggray, pgcont, pgimag):
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
end subroutine bindata2dold
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief Bin data in 2 dimensions - computing the bin number rather than searching for it is ~10x faster
!! 
!! \param n  Number of input data points (integer)
!! \param x  Input data: x values (real)
!! \param y  Input data: y values (real)
!! \param norm  Normalise the bins (1) or not (0) (integer)
!! \param nxbin  Desired number of bins in the x direction
!! \param nybin  Desired number of bins in the y direction
!! \param xmin1  Lower limit for the binning range in the x direction - autodetermine if xmin1=xmax1 (real)
!! \param xmax1  Upper limit for the binning range in the x direction - autodetermine if xmin1=xmax1 (real)
!! \param ymin1  Lower limit for the binning range in the y direction - autodetermine if ymin1=ymax1 (real)
!! \param ymax1  Upper limit for the binning range in the y direction - autodetermine if ymin1=ymax1 (real)
!! \retval z     Binned data set z(nxbin,nybin) (real)
!! \retval tr    Transformation elements for pgplot tr(6) (real)

subroutine bindata2d(n,x,y,norm,nxbin,nybin,xmin1,xmax1,ymin1,ymax1,z,tr)
  
  implicit none
  integer, intent(in) :: n, nxbin,nybin, norm
  real, intent(in) :: x(n),y(n)
  real, intent(inout) :: xmin1,xmax1,ymin1,ymax1
  real, intent(out) :: z(nxbin+1,nybin+1),tr(6)
  integer :: i,bx,by
  real :: xmin,xmax,ymin,ymax,dx,dy
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  
  if(abs((xmin-xmax)/(xmax+1.e-30)).lt.1.e-20) then !Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs((ymin-ymax)/(ymax+1.e-30)).lt.1.e-20) then !Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  
  
  
  ! Determine transformation elements for pgplot (pggray, pgcont, pgimag):
  tr(1) = xmin - dx/2.
  tr(2) = dx
  tr(3) = 0.
  tr(4) = ymin - dy/2.
  tr(5) = 0.
  tr(6) = dy
  
  z = 0.
  do i=1,n
     bx = floor((x(i) - xmin)/dx) + 1 
     by = floor((y(i) - ymin)/dy) + 1
     !if(bx.lt.1.or.bx.gt.nxbin.or.by.lt.1.or.by.gt.nybin) then
     !   if(bx.eq.0.or.bx.eq.nxbin+1) bx = max(min(bx,nxbin),1)  !Treat an error of 1 x bin as round-off
     !   if(by.eq.0.or.by.eq.nybin+1) by = max(min(by,nybin),1)  !Treat an error of 1 y bin as round-off
     !   
     !   if(bx.lt.0.or.bx.gt.nxbin+1) then
     !      !write(stdErr,'(A,I7,A2,F8.3,A,I4,A,I4,A1)') &
     !'  Bindata2d:  error for X data point',i,' (',x(i),').  I found bin',bx,', but it should lie between 1 and',nxbin,'.'
     !   else if(by.lt.0.or.by.gt.nybin+1) then
     !      !write(stdErr,'(A,I7,A2,F8.3,A,I4,A,I4,A1)') &
     !'  Bindata2d:  error for Y data point',i,' (',y(i),').  I found bin',by,', but it should lie between 1 and',nybin,'.'
     !   else
     !      z(bx,by) = z(bx,by) + 1.
     !   end if
     !else
     !   z(bx,by) = z(bx,by) + 1.
     !end if
     if(bx.ge.1.and.bx.le.nxbin.and.by.ge.1.and.by.le.nybin) z(bx,by) = z(bx,by) + 1.  ! Don't treat 1-bin errors as round-off
  end do
  
  !if(norm.eq.1) z = z/(ztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs((xmin1-xmax1)/(xmax1+1.e-30)).lt.1.e-20) then  ! Autodetermine
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs((ymin1-ymax1)/(ymax1+1.e-30)).lt.1.e-20) then  ! Autodetermine
     ymin1 = ymin
     ymax1 = ymax
  end if
  
end subroutine bindata2d
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief 'Bin data' in 2 dimensions  -  Measure the amount of likelihood in each bin
!! 
!! \param n      Number of input data points (integer)
!! \param x      Input data: x values (real)
!! \param y      Input data: y values (real)
!! \param z      Input data: likelihood values
!!
!! \param norm   Normalise the bins (1) or not (0) (integer)
!! \param nxbin  Desired number of bins in the x direction
!! \param nybin  Desired number of bins in the y direction
!!
!! \param xmin1  Lower limit for the binning range in the x direction - autodetermine if xmin1=xmax1 (real)
!! \param xmax1  Upper limit for the binning range in the x direction - autodetermine if xmin1=xmax1 (real)
!! \param ymin1  Lower limit for the binning range in the y direction - autodetermine if ymin1=ymax1 (real)
!! \param ymax1  Upper limit for the binning range in the y direction - autodetermine if ymin1=ymax1 (real)
!!
!! \retval zz    'Binned' data set z(nxbin,nybin) (real)
!! \retval tr    Transformation elements for pgplot tr(6) (real)
!!
!! \todo  Should z be replaced?

subroutine bindata2da(n, x,y,z, norm, nxbin,nybin, xmin1,xmax1,ymin1,ymax1, zz, tr)
  
  implicit none
  integer, intent(in) :: n, norm, nxbin,nybin
  integer :: i,bx,by
  real, intent(in) :: x(n),y(n)
  real, intent(inout) :: z(n), xmin1,xmax1,ymin1,ymax1
  real, intent(out) :: zz(nxbin+1,nybin+1), tr(6)
  
  real :: xbin(nxbin+1),ybin(nybin+1), zztot, xmin,xmax,ymin,ymax, dx,dy ,zmin
  
  !write(stdOut,'(A4,5I8)')'n:',norm,nxbin,nybin
  !write(stdOut,'(A4,2F8.3)')'x:',xmin1,xmax1
  !write(stdOut,'(A4,2F8.3)')'y:',ymin1,ymax1
  
  xmin = xmin1
  xmax = xmax1
  ymin = ymin1
  ymax = ymax1
  zmin = minval(z)
  
  if(abs((xmin-xmax)/(xmax+1.e-30)).lt.1.e-20) then  ! Autodetermine
     xmin = minval(x(1:n))
     xmax = maxval(x(1:n))
  end if
  dx = abs(xmax - xmin)/real(nxbin)
  if(abs((ymin-ymax)/(ymax+1.e-30)).lt.1.e-20) then  ! Autodetermine
     ymin = minval(y(1:n))
     ymax = maxval(y(1:n))
  end if
  dy = abs(ymax - ymin)/real(nybin)
  do bx=1,nxbin+1
     !xbin(bx) = xmin + (real(bx)-0.5)*dx  ! x is the centre of the bin
     xbin(bx) = xmin + (bx-1)*dx           ! x is the left of the bin
  end do
  do by=1,nybin+1
     !ybin(by) = ymin + (real(by)-0.5)*dy  ! y is the centre of the bin
     ybin(by) = ymin + (by-1)*dy           ! y is the left of the bin
  end do
  
  !write(stdOut,'(50F5.2)'),x(1:50)
  !write(stdOut,'(50F5.2)'),y(1:50)
  !write(stdOut,'(20F8.5)'),xbin
  !write(stdOut,'(20F8.5)'),ybin
  
  zz = 0.
  zztot = 0.
  do bx=1,nxbin
     do by=1,nybin
        zz(bx,by) = 0.
        do i=1,n
           !if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1))  &
           !zz(bx,by) = zz(bx,by) + 1.
           if(x(i).ge.xbin(bx).and.x(i).lt.xbin(bx+1) .and. y(i).ge.ybin(by).and.y(i).lt.ybin(by+1))  &
                zz(bx,by) = zz(bx,by) + exp(z(i) - zmin)
           !write(stdOut,'(2I4,8F10.5)')bx,by,x(i),xbin(bx),xbin(bx+1),y(i),ybin(by),ybin(by+1),zz(bx,by),z(i)
        end do
        zztot = zztot + zz(bx,by) 
        !write(stdOut,'(2I4,5x,4F6.3,5x,10I8)')bx,by,xbin(bx),xbin(bx+1),ybin(by),ybin(by+1),nint(zz(bx,by))
     end do
     !write(stdOut,'(I4,5x,2F6.3,5x,10I8)')bx,xbin(bx),xbin(bx+1),nint(zz(bx,1:nybin))
  end do
  !if(norm.eq.1) z = z/(zztot+1.e-30)
  if(norm.eq.1) z = z/maxval(z+1.e-30)
  
  if(abs((xmin1-xmax1)/(xmax1+1.e-30)).lt.1.e-20) then  !Autodetermine
     xmin1 = xmin
     xmax1 = xmax
  end if
  if(abs((ymin1-ymax1)/(ymax1+1.e-30)).lt.1.e-20) then  !Autodetermine
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
  
end subroutine bindata2da
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
  use constants, only: rd2r
  use analysemcmc_settings, only: changeVar,ivals,prProgress
  use mcmcrun_data, only: parID
  
  implicit none
  integer, intent(in) :: p1,p2,ni,nx,ny
  real, intent(inout) :: z(nx,ny)
  real, intent(in) :: tr(6)
  
  integer :: nn,indx(nx*ny),i,b,ib,full(ni),iy
  real :: x1(nx*ny),x2(nx*ny),tot,np,y
  
  
  ! Weigh number of points in each bin by bin size for position/orientation plots:
  do iy = 1,ny
     if(changeVar.ge.1) then
        if((parID(p1).eq.31.and.parID(p2).eq.32) .or. (parID(p1).eq.52.and.parID(p2).eq.51)) then  
           !Then: RA-Dec or (phi/theta_Jo)/(psi/i) plot, convert lon -> lon * 15 * cos(lat)
           y = tr(4) + tr(6)*iy
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
     end if !if(changeVar.ge.1)
  end do !iy
  
  
  nn = nx*ny
  x1 = reshape(z,(/nn/))                 ! x1 is an 1D array with the same data as the 2D array z
  call rindexx(nn,-x1(1:nn),indx(1:nn))  ! -x1: sort the 1D array to descending value
  
  np = sum(z)
  tot = 0.
  full = 0
  do b=1,nn !Loop over bins in 1D array
     ib = indx(b)
     x2(ib) = 0.
     if(x1(ib).eq.0.) cycle
     tot = tot + x1(ib)
     do i=ni,1,-1 !Loop over intervals
        if(tot.le.np*ivals(i)) then
           x2(ib) = real(ni-i+1)  !e.g. x2(b) = ni if within 68%, ni-1 if within 95%, etc, and 1 if within 99.7%
        else
           if(prProgress.ge.3.and.full(i).eq.0) then !Report the number of points in the lastly selected bin
              if(i.eq.1) write(stdOut,'(A)',advance="no")'Last bin:'
              !write(stdOut,'(F6.3,I5)',advance="no")ivals(i),nint(x1(ib))
              write(stdOut,'(I5)',advance="no")nint(x1(ib))
              full(i) = 1
           end if
        end if
        !write(stdOut,'(2I4, F6.2, 3F20.5)')b,i, ivals(i), np,tot,np*ivals(i)
     end do
  end do
  
  z = reshape(x2, (/nx,ny/))  ! z lies between 1 and ni
  
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
  use constants, only: rd2r
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
              y = tr(4) + tr(6)*iy
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

