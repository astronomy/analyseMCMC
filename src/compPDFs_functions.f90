!> \file compPDFs_functions.f90  Plot functions for compPDFs  (1D/2D PDFs and waveforms)

! 
! LICENCE:
! 
! Copyright (c) 2007-2014  Marc van der Sluys
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



!> \todo put file reading in separate routines


!***********************************************************************************************************************************
module comp_pdfs_settings
  implicit none
  save
  integer, parameter :: nfmax=9       ! nfmax:    maximum number of input files
  integer, parameter :: nParDB=99     ! nParDB: size of the parameter database
  
  integer :: nf,file,type,dim,clr,fillstyle,frames(2),plpars(20),plpars2d(2),clrs(nfmax),lss(nfmax)
  integer :: pltrue,plmedian,plrange,fonttype
  real :: fontsize
  character :: fnames(nfmax)*(299),dirnames(nfmax)*(199),outnamebase*(99),settingsfile*(99)
  
end module comp_pdfs_settings
!***********************************************************************************************************************************


!***********************************************************************************************************************************
!> \brief  Plot a 1D marginalised PDF
!!
!! \param pp   Parameter ID
!! \param lbl  Parameter/plot label

subroutine plotpdf1d(pp,lbl)
  use SUFR_numerics, only: sne
  use comp_pdfs_settings, only: nf,fnames,fontsize,clrs,clr,fillstyle,plrange,pltrue,lss,plmedian
  use general_data, only: parNames,pgunits,pgParNs
  
  implicit none
  integer, intent(in) :: pp
  character, intent(in) :: lbl*(99)
  
  integer, parameter :: np=99,nbin1=500  !np: number of parameters (currently 11-87 are defined)
  integer :: b,p1,io,f,nplvar,nchains,nbin(nf),pID,parIDs(np),pp1,ic,wrap(nf,np),lw,detnan(nf,np),identical
  real :: x,startval(nf,np,2),stats(nf,np,6),ranges(nf,np,5),xmin1(nf,np),xmax1(nf,np),plshift
  real :: xbin1(nf,np,nbin1),ybin1(nf,np,nbin1),xbin(nbin1),ybin(nbin1),xmin,xmax,ymin,ymax,dx,yrange(2),xpeak
  character :: fname*(299),str*(99),tmpstr
  
  
  detnan = 0  ! Used to detect NaNs
  pp1 = 0
  xpeak = -huge(xpeak)
  
  do f=1,nf
     fname = trim(fnames(f))
     if(fname(1:3).eq.'   ') cycle
     open(unit=10,action='read',form='formatted',status='old',position='rewind',file=trim(fname),iostat=io)
     if(io.ne.0) then
        write(*,'(A,I5)')'  Error opening file '//trim(fname)//', error:',io
        write(*,'(A,/)')'  Aborting.'
        stop
     end if
     read(10,'(3I6)') nplvar,nchains,nbin(f)  ! Total number of plot parameters, total number of chains, number of bins
     pp1 = 0
     do p1=1,nplvar
        read(10,*,iostat=io) tmpstr  ! Read empty line
        if(io.ne.0) then  ! Wrong number of parameters listed in top of file?
           !write(*,'(A,2(I0,A))') '  The number of parameters listed at the top of the file '//trim(fname)//' (', &
           !     nplvar,') seems to be incorrect. I found only ',p1-1,'.'
           exit
        end if
        read(10,'(3I6)')    ic,parIDs(p1),wrap(f,p1) !'Chain number, parameter number, and wrap(type)'
        read(10,'(2E15.7)') startval(f,p1,1:2) !'True and starting value'
        read(10,'(6E15.7)') stats(f,p1,1:6) !'Stats: median, mean, absvar1, absvar2, stdev1, stdev2'
        read(10,'(5E15.7)') ranges(f,p1,1:5) !'Ranges: lower,upper limit, centre, width, relative width'
        read(10,'(2E15.7)') xmin1(f,p1),xmax1(f,p1) !'Xmin and Xmax of PDF'
        do b=1,nbin(f)+1
           !read(10,'(2E15.7)',iostat=io)xbin1(f,p1,b),ybin1(f,p1,b)
           read(10,*,iostat=io) xbin1(f,p1,b),ybin1(f,p1,b)  ! Formatted read doesn't work for gfortran when NaNs are present
           if(io.ne.0) then
              if(detnan(f,p1).eq.0) write(*,'(A,I4,A1,I4,A1)') '  Warning while reading file '//trim(fname)//', parameter '// &
                   trim(parNames(parIDs(p1)))//', bin',b,'/',nbin(f),'.'
              ybin1(f,p1,b) = 0.0
              detnan(f,p1) = 1
           end if
           ! Then it's probably a NaN:
           !if(parIDs(p1).eq.pp .and. (.not.ybin1(f,p1,b).gt.-1.e30) .and. (.not.ybin1(f,p1,b).lt.1.e30)) then
           if(parIDs(p1).eq.pp .and. sne(ybin1(f,p1,b),ybin1(f,p1,b)) .and. sne(ybin1(f,p1,b),ybin1(f,p1,b))) then
              detnan(f,p1) = 1
              ybin1(f,p1,b) = 0.
           end if
        end do
        if(parIDs(p1).eq.pp) then
           pp1 = p1
        end if
     end do  ! p1
     
     close(10)
  end do  ! f
  
  
  if(pp1.eq.0) then
     write(0,'(A,I4,A)')'  Variable',pp,' not found!'
     return
  end if
  
  if(sum(detnan).gt.0) write(*,'(A)')'  Warning:  I think I detected NaNs and set them to zero in the PDF for ' &
       //trim(parNames(parIDs(pp)))//'.'
  
  
  
  pID = parIDs(pp1)  !pp1 gives the order in the file (e.g. pp1=1 is first parameter), p1 gives parameter ID (e.g. parID=61: Mc)
  xmin = minval(xmin1(1:nf,pp1))
  xmax = maxval(xmax1(1:nf,pp1))
  
  
  dx = abs(xmax-xmin)
  if(dx.lt.1.e-30) dx = xmin !xmin=xmax
  if(dx.lt.1.e-30) dx = 0.1  !If it's still zero (i.e. xmin=xmax=0)
  xmin = xmin - 0.1*dx
  xmax = xmax + 0.1*dx
  ymin = 0.
  ymax = -1.e30
  do f=1,nf
     !ymax = max(ymax,maxval(ybin1(f,pp1,1:nbin(f)+1)))
     do b=1,nbin(f)+1
        if(ybin1(f,pp1,b).gt.ymax) then
           ymax = ybin1(f,pp1,b)
           xpeak = xbin1(f,pp1,b)
        end if
     end do
  end do
  !print*,f,pp1!,xpeak
  
  ymax = ymax*1.1
  if(nf.eq.1) ymax = ymax*1.2 !Make room for numbers
  if(ymax.lt.1.e-30) ymax = 1.
  
  yrange = (/-1.e30,1.e30/) !Used to plot true value, probability range
  if(nf.eq.1) yrange = (/-1.e30,ymax*0.9/)
  
  
  lw = 3 !Line width for pdf contours, not hatches
  call pgslw(lw)
  call pgsch(fontsize*0.5)
  call pgswin(xmin,xmax,ymin,ymax)
  
  !if(nf.gt.1) call pgsfs(3) !Hatches
  call pgsfs(fillstyle)
  do f=1,nf
     !if(nf.gt.1.and.(clr.eq.0.or.clr.eq.2).and.fillstyle.eq.2) call pgsls(mod(f,4))
     xbin = xbin1(f,pp1,:)
     ybin = ybin1(f,pp1,:)
     call pgshs(45.,1.5,0.)                ! Hatches slanted up
     if(f.eq.2) call pgshs(-45., 1.5, 0.)  ! Hatches slanted down
     if(f.eq.3) call pgshs(  0., 1.5, 0.)  ! Hatches horizontal
     if(f.eq.4) call pgshs( 90., 1.5, 0.)  ! Hatches vertical
     
     call pgsci(clrs(f))
     if(nf.eq.1) call pgsci(15)
     
     if(wrap(f,pp1).eq.0) then
        !if(clr.eq.0) call pgsci(15)
        call pgslw(1)
        call pgpoly(nbin(f)+2,(/xbin(1),xbin(1:nbin(f)+1)/),(/0.,ybin(1:nbin(f)+1)/))
        
        !Plot pdf contour
        if(nf.eq.1.and.clr.eq.1) call pgsci(2)
        if(clr.eq.0.or.clr.eq.2) call pgsci(1)
        call pgslw(lw)
        call pgline(nbin(f)+1,xbin(1:nbin(f)+1),ybin(1:nbin(f)+1)) !:nbin(f)) ?
        
        !Fix the loose ends
        call pgline(2,(/xbin(1),xbin(1)/),(/0.,ybin(1)/))
        call pgline(2,(/xbin(nbin(f)+1),xbin(nbin(f)+1)/),(/ybin(nbin(f)+1),0./))
     else
        !plshift = real(2*pi)
        plshift = 360.
        !if(clr.eq.0) call pgsci(15)
        call pgslw(1)
        call pgpoly(nbin(f)+3,(/xbin(1),xbin(1:nbin(f)),xbin(1)+plshift,xbin(1)+plshift/),(/0.,ybin(1:nbin(f)),ybin(1),0./))
        
        !Plot pdf contour
        if(nf.eq.1.and.clr.eq.1) call pgsci(2)
        if(clr.eq.0.or.clr.eq.2) call pgsci(1)
        !if(clr.eq.0) call pgsci(1)
        call pgslw(lw)
        call pgline(nbin(f),xbin(1:nbin(f)),ybin(1:nbin(f)))
        
        !Plot dotted lines outside the pdf for wrapped periodic parameters
        call pgsls(4)
        call pgline(nbin(f)+1,(/xbin(1:nbin(f))-plshift,xbin(1)/),(/ybin(1:nbin(f)),ybin(1)/))
        call pgline(nbin(f),xbin+plshift,ybin)
        
        !Fix the loose end
        call pgsls(1)
        call pgslw(lw)
        call pgline(2,(/xbin(nbin(f)),xbin(1)+plshift/),(/ybin(nbin(f)),ybin(1)/))
     end if
     
     !plot probability ranges
     if(plrange.ge.1.and.nf.eq.1) then
        !Plot limits
        call pgsls(4)
        if(clr.eq.1) call pgsci(2)
        if(clr.eq.0.or.clr.eq.2) call pgsci(14)
        !if(clr.eq.0) call pgsci(1)
        call pgline(2,(/ranges(f,pp1,1),ranges(f,pp1,1)/),yrange)
        call pgline(2,(/ranges(f,pp1,2),ranges(f,pp1,2)/),yrange)
        
        !Print number
        if(plrange.ge.2) then
           x = ranges(f,pp1,5)
           if(pId.eq.61.or.pId.eq.63.or.pId.eq.64.or.pId.eq.21.or.pId.eq.22.or.pId.eq.71.or.pId.eq.81) x = x*100
           write(str,'(F10.3)')x
           if(pId.eq.61.or.pId.eq.63.or.pId.eq.64.or.pId.eq.21.or.pId.eq.22.or.pId.eq.71.or.pId.eq.81) &
                write(str,'(A)')trim(str)//'%'
           
           if(x.lt.0.01) write(str,'(F6.4)')x
           if(x.ge.0.01.and.x.lt.0.1) write(str,'(F5.3)')x
           if(x.ge.0.1.and.x.lt.1) write(str,'(F4.2)')x
           if(x.ge.1.and.x.lt.10) write(str,'(F3.1)')x
           if(x.ge.10.and.x.lt.100) write(str,'(I2)')nint(x)
           if(x.ge.100) write(str,'(I3)')nint(x)
           write(str,'(A)')'\(2030): '//trim(str)
           if(pId.eq.61.or.pId.eq.63.or.pId.eq.64.or.pId.eq.21.or.pId.eq.22.or.pId.eq.71.or.pId.eq.81) then
              write(str,'(A)')trim(str)//'%'
           else
              write(str,'(A)')trim(str)//trim(pgunits(pp1))
           end if
           !call pgsch(fontsize*1.2)
           !call pgptxt(xmin+0.05*dx,ymax*0.9,0.,0.,trim(pgParNs(pID)))
           call pgptxt(ranges(f,pp1,3),yrange(2)*1.05,0.,0.5,trim(str))
           call pgsls(1)
           call pgline(2,(/ranges(f,pp1,1),ranges(f,pp1,2)/),(/yrange(2),yrange(2)/))
        end if
     end if
     call pgsls(1)
  end do !f
  
  !Plot lines again over surface of overlapping distributions
  if(1.eq.2) then
     call pgsls(4)
     !if(clr.eq.0) call pgsci(1)
     do f=1,nf
        xbin = xbin1(f,pp1,:)
        ybin = ybin1(f,pp1,:)
        if(wrap(f,pp1).eq.0) then
           call pgline(nbin(f)+1,xbin(1:nbin(f)+1),ybin(1:nbin(f)+1))
        else
           call pgline(nbin(f),xbin(1:nbin(f)),ybin(1:nbin(f)))
        end if
     end do
     call pgsls(1)
  end if
  
  
  
  !Plot true values
  if(pltrue.ge.1) then
     identical = 1
     if(nf.gt.1) then
        do f=2,nf
           if(abs(startval(f,pp1,1)-startval(1,pp1,1)).gt.1.e-10) identical = 0
        end do
     end if
     
     if(identical.eq.1) then
        call pgsci(1); call pgslw(lw);call pgsls(2) !; if(clr.eq.0) call pgsci(1)
        !if(nf.eq.1) call pgsci(2)
        ! If a=0, don't plot theta_SL (pp1=6):
        if(.not.(pp1.eq.6.and.abs(startval(1,5,1)).lt.0.001)) call pgline(2,(/startval(1,pp1,1),startval(1,pp1,1)/),yrange)
     else
        do f=1,nf
           if(pp1.eq.6.and.abs(startval(f,5,1)).lt.0.001) cycle  !If a=0, don't plot theta_SL (pp1=6)
           call pgslw(lw); call pgsls(2); call pgsci(1)
           if(clr.eq.1) call pgsci(clrs(f))
           if(pltrue.eq.2) call pgsls(lss(f))
           !call pgsci(clrs(f))  !it seems clearer when all true values are white
           call pgline(2,(/startval(f,pp1,1),startval(f,pp1,1)/),yrange)
        end do
     end if
  end if !if(pltrue.ge.1)
  
  
  
  
  !Plot median (if one file)
  if(plmedian.eq.1.and.nf.eq.1) then
     if(clr.eq.1) call pgsci(2); call pgslw(lw);call pgsls(2)  !; if(clr.eq.0) call pgsci(1)
     call pgline(2,(/stats(1,pp1,1),stats(1,pp1,1)/),yrange)
  end if !if(plmedian.eq.1.and.nf.eq.1)
  
  
  
  
  call pgsls(1)
  call pgsci(1)
  !if(nf.eq.1) then
  call pgbox('BNTS',0.0,0,'',0.0,0)
  !else
  !   call pgbox('BCNTS',0.0,0,'BC',0.0,0)
  !end if
  
  !call pgmtxt('B',2.4,0.5,0.5,trim(pgParNs(pID)))  !Plot label under x axis
  if(abs(xpeak-xmin).gt.abs(xpeak-xmax)) then  !Peak is right, plot varname left
     call pgmtxt('T',-1.,0.05,0.,trim(pgParNs(pID)))  !Plot label in upper-left corner
     call pgmtxt('T',-1.5,0.95,1.,trim(lbl))
  else   !Peak is left, plot varname right
     call pgmtxt('T',-1.,0.95,1.,trim(pgParNs(pID)))  !Plot label in upper-left corner
     call pgmtxt('T',-1.5,0.05,0.,trim(lbl))
  end if
  
  
end subroutine plotpdf1d
!***********************************************************************************************************************************











!***********************************************************************************************************************************
!> \brief  Plot a 2D marginalised PDF
!!
!! \param pID1  Parameter ID 1
!! \param pID2  Parameter ID 2
!! \param lbl   Plot label

subroutine plotpdf2d(pID1,pID2,lbl)
  use comp_pdfs_settings, only: nf,fnames,fontsize,clrs,clr,fillstyle,plrange,pltrue,plmedian
  use general_data, only: parNames,pgunits,pgParNs
  
  implicit none
  integer, intent(in) :: pID1,pID2
  character, intent(in) :: lbl*(99)
  
  integer, parameter :: np=15,nbinx1=500,nbiny1=500
  integer :: bx,by,pID1a,pID2a,p11,p22,pp11,pp22,pp12,p12,io,f,nplvar,nplvar1,nplvar2
  integer :: nchains,nbinx(nf),nbiny(nf),ic,lw,c,foundit
  integer :: identical
  real :: startval(nf,np,2,2),stats(nf,np,2,6),ranges(nf,np,2,5)
  real :: xmin1(nf,np,2),xmax1(nf,np,2),ymin1(nf,np,2),ymax1(nf,np,2),x
  real :: xmin,xmax,ymin,ymax,dx,dy,z(nf,nbinx1,nbiny1),z1(nbinx1,nbiny1),tr(nf,np*np,6),cont(11)
  character :: fname*(299),str*(99),tmpstr
  
  
  pp12 = 0
  dof: do f=1,nf
     foundit = 0
     fname = trim(fnames(f))
     if(fname(1:3).eq.'   ') cycle
     open(unit=10,action='read',form='formatted',status='old',position='rewind',file=trim(fname),iostat=io)
     if(io.ne.0) then
        write(*,'(A,I5)')'  Error opening file '//trim(fname)//'.'
        cycle
     end if
     read(10,'(5I6)') nplvar1,nplvar2,nchains,nbinx(f),nbiny(f)  ! Plot parameter 1,2, total number of chains, number of bins x,y
     nplvar = nplvar2-nplvar1+1
     p12 = 0
     
     dop11: do p11=nplvar1,nplvar2-1
        do p22=p11+1,nplvar2
           p12 = p12+1
           read(10,*,iostat=io) tmpstr
           if(io.ne.0) then
              write(0,'(A)') '  End of file '//trim(fname)//' reached, parameter combination '//trim(parNames(pID1))//'-'// &
                   trim(parNames(pID2))//' not found,  skipping...'
              cycle dof
           end if
           read(10,'(3I6)')    ic,pID1a,pID2a         ! 'Chain number and parameter ID 1,2'
           read(10,'(2E15.7)') startval(f,p12,1,1:2)  ! 'True and starting value p1'
           read(10,'(2E15.7)') startval(f,p12,2,1:2)  ! True and starting value p2'
           read(10,'(6E15.7)') stats(f,p12,1,1:6)     ! Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p1'
           read(10,'(6E15.7)') stats(f,p12,2,1:6)     ! Stats: median, mean, absvar1, absvar2, stdev1, stdev2 for p2'
           
           ! The limits are NOT used as plot limits!!!:
           read(10,'(5E15.7)') ranges(f,p12,1,1:5)    ! Ranges: lower,upper limit, centre, width, relative width for p1'
           read(10,'(5E15.7)') ranges(f,p12,2,1:5)    ! Ranges: lower,upper limit, centre, width, relative width for p2'
           read(10,'(4E15.7)') xmin1(f,p12,1),xmax1(f,p12,1),ymin1(f,p12,2),ymax1(f,p12,2)  ! 'Xmin,Xmax,Ymin,Ymax of PDF'
           read(10,'(6E15.7)') tr(f,p12,1:6)          ! 'Tr'              
           !print*,tr
           read(10,*) tmpstr
           do bx=1,nbinx(f)+1
              do by=1,nbiny(f)
                 read(10,'(E15.7)',advance='no')z1(bx,by)
              end do
              read(10,'(E15.7)')z1(bx,by)
           end do
           
           if(pID1a.eq.pID1.and.pID2a.eq.pID2) then  ! Found the desired parameter combination
              pp12 = p12
              pp11 = pID1a
              pp22 = pID2a
              z(f,:,:) = z1(:,:)
              exit dop11
           end if
           
        end do !p22
     end do dop11 !p11
     
     close(10)
  end do dof !f
  
  
  if(pp12.eq.0) then
     write(0,'(/,A,/)')'  No parameter combinations found,  aborting...'
     stop
  end if
  
  
  
  xmin = minval(xmin1(1:nf,pp12,1))
  xmax = maxval(xmax1(1:nf,pp12,1))
  ymin = minval(ymin1(1:nf,pp12,2))
  ymax = maxval(ymax1(1:nf,pp12,2))
  
  if(pID1.eq.31) then  !RA
     xmin = maxval(xmin1(1:nf,pp12,1))
     xmax = minval(xmax1(1:nf,pp12,1))
  end if
  
  if(1.eq.2.and.pID1.eq.31.and.pID2.eq.31) then !Get whole sky for sky plot (done in analysemcmc?)
     xmin = 24.
     xmax = 0.
     ymin = -90.
     ymax = 90.
  end if
  
  dx = xmax - xmin
  dy = ymax - ymin
  if(abs(dx).lt.1.e-30) dx = xmax  !If xmin=xmax
  if(abs(dx).lt.1.e-30) dx = 0.1   !If xmin=xmax=0
  if(abs(dy).lt.1.e-30) dy = ymax  !If ymin=ymax
  if(abs(dy).lt.1.e-30) dy = 0.1   !If ymin=ymax=0
  
  lw = 3
  call pgsls(1)
  call pgsci(1)
  call pgslw(lw)
  call pgsch(fontsize)
  call pgswin(xmin,xmax,ymin,ymax)
  call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
  call pgmtxt('B',2.4,0.5,0.5,trim(pgParNs(pID1)))  !Cheat a bit
  call pgmtxt('L',2.2,0.5,0.5,trim(pgParNs(pID2)))  !Cheat a bit
  
  !select contours
  do c=1,11
     cont(c) = 0.01 + 2*real(c-1)/10.
  end do
  
  do f=1,nf
     ! Fill contours with hatches
     call pgsci(1)
     !call pgsfs(3)  ! Hatches
     call pgsfs(fillstyle)
     call pgshs(45.,1.5,0.)                ! Hatches slanted up
     if(f.eq.2) call pgshs(-45., 1.5, 0.)  ! Hatches slanted down
     if(f.eq.3) call pgshs(  0., 1.5, 0.)  ! Hatches horizontal
     if(f.eq.4) call pgshs( 90., 1.5, 0.)  ! Hatches vertical
     
     !Fill contours with solid colour, override hatches
     !call pgsci(15) !Light grey
     !call pgsfs(1)
     
     call pgslw(1)
     if(f.eq.2) call pgsci(14) !Dark grey
     !if(f.eq.2) call pgshs(-45.,0.7,0.) !Hatches down
     call pgsci(clrs(f))
     call pgconf(z(f,1:nbinx(f)+1,1:nbiny(f)+1),nbinx(f)+1,nbiny(f)+1,1,nbinx(f)+1,1,nbiny(f),cont(1),cont(11),tr(f,pp12,1:6))
     !call pgsci(1) !When not using hatches
     call pgslw(lw)
     call pgcont(z(f,1:nbinx(f)+1,1:nbiny(f)+1),nbinx(f)+1,nbiny(f)+1,1,nbinx(f)+1,1,nbiny(f),cont(1:1),1,tr(f,pp12,1:6))
     
     
     !plot probability ranges
     !if(plrange.eq.1) then
     if(plrange.ge.1.and.nf.eq.1) then
        !Plot limits
        call pgsls(1)
        if(clr.eq.1) call pgsci(2)
        call pgsci(clrs(f))
        x = 0.05+f*0.02
        call pgarro(ranges(f,pp12,1,3),ymin+dy*x,ranges(f,pp12,1,1),ymin+dy*x)
        call pgarro(ranges(f,pp12,1,3),ymin+dy*x,ranges(f,pp12,1,2),ymin+dy*x)
        call pgarro(xmin+dx*x*0.7,ranges(f,pp12,2,3),xmin+dx*x*0.7,ranges(f,pp12,2,1))
        call pgarro(xmin+dx*x*0.7,ranges(f,pp12,2,3),xmin+dx*x*0.7,ranges(f,pp12,2,2))
        
        !Print number
        !if(nf.eq.1) then
        if(plrange.ge.2.and.nf.eq.1) then
           x = ranges(f,pp12,1,5)
           if(pId1.eq.61.or.pId1.eq.63.or.pId1.eq.64.or.pId1.eq.21.or.pId1.eq.22.or.pId1.eq.71.or.pId1.eq.81) x = x*100
           write(str,'(F10.3)')x
           if(pId1.eq.61.or.pId1.eq.63.or.pId1.eq.64.or.pId1.eq.21.or.pId1.eq.22.or.pId1.eq.71.or.pId1.eq.81) &
                write(str,'(A)')trim(str)//'%'
           
           if(x.lt.0.01) write(str,'(F6.4)')x
           if(x.ge.0.01.and.x.lt.0.1) write(str,'(F5.3)')x
           if(x.ge.0.1.and.x.lt.1) write(str,'(F4.2)')x
           if(x.ge.1.and.x.lt.10) write(str,'(F3.1)')x
           if(x.ge.10.and.x.lt.100) write(str,'(I2)')nint(x)
           if(x.ge.100) write(str,'(I3)')nint(x)
           write(str,'(A)')'\(2030): '//trim(str)
           if(pId1.eq.61.or.pId1.eq.63.or.pId1.eq.64.or.pId1.eq.21.or.pId1.eq.22.or.pId1.eq.71.or.pId1.eq.81) then
              write(str,'(A)')trim(str)//'%'
           else
              write(str,'(A)')trim(str)//trim(pgUnits(pID1))
           end if
           !call pgptxt(ranges(f,pp12,1,3),yrange(2)*1.05,0.,0.5,trim(str))
           call pgsls(1)
           !call pgline(2,(/ranges(f,pp12,1,1),ranges(f,pp12,1,2)/),(/yrange(2),yrange(2)/))
        end if
     end if !if(plrange.ge.1.and.nf.eq.1)
     call pgsls(1)
     
  end do !do f=1,nf
  
  
  !Plot true values 
  if(pltrue.ge.1) then
     identical = 1
     if(nf.gt.1) then
        do f=2,nf
           if(abs(startval(f,pp12,1,1)-startval(1,pp12,1,1)).gt.1.e-10) identical = 0
        end do
     end if
     ! Then all have the same value:
     !if(abs((startval(1,pp12,1,1) - sum(startval(1:nf,pp12,1,1))/(real(nf)))/startval(1,pp12,1,1)).lt.0.001) then
     if(identical.eq.1) then
        call pgslw(lw);call pgsls(2); call pgsci(1)
        call pgline(2,(/startval(1,pp12,1,1),startval(1,pp12,1,1)/),(/-1.e20,1.e20/))
     else
        do f=1,nf
           call pgslw(lw);call pgsls(2)
           call pgsci(1)
           !call pgsci(clrs(f))  !it seems clearer when all true values are white
           call pgline(2,(/startval(f,pp12,1,1),startval(f,pp12,1,1)/),(/-1.e20,1.e20/))
        end do
     end if
     
     identical = 1
     do f=2,nf
        if(abs(startval(f,pp12,2,1)-startval(1,pp12,2,1)).gt.1.e-10) identical = 0
     end do
     ! Then all have the same value:
     !if(abs((startval(1,pp12,2,1) - sum(startval(1:nf,pp12,2,1))/(real(nf)))/startval(1,pp12,2,1)).lt.0.001) then
     if(identical.eq.1) then
        call pgslw(lw);call pgsls(2); call pgsci(1)
        call pgline(2,(/-1.e20,1.e20/),(/startval(1,pp12,2,1),startval(1,pp12,2,1)/))
     else
        do f=1,nf
           call pgslw(lw);call pgsls(2)
           call pgsci(1)
           !call pgsci(clrs(f))  !it seems clearer when all true values are white
           call pgline(2,(/-1.e20,1.e20/),(/startval(f,pp12,2,1),startval(f,pp12,2,1)/))
        end do
     end if
  end if !if(pltrue.ge.1)
  
  
  !Plot median (if one file)
  if(plmedian.eq.1.and.nf.eq.1) then
     if(clr.eq.1) call pgsci(2); call pgslw(lw);call pgsls(2)  !; if(clr.eq.0) call pgsci(1)
     call pgline(2,(/stats(1,pp12,1,1),stats(1,pp12,1,1)/),(/-1.e20,1.e20/))
     call pgline(2,(/-1.e20,1.e20/), (/stats(1,pp12,2,1),stats(1,pp12,2,1)/) )
  end if !if(plmedian.eq.1.and.nf.eq.1)
  
  
  call pgmtxt('T',-1.5,0.05,0.,trim(lbl))
  
end subroutine plotpdf2d
!***********************************************************************************************************************************










!***********************************************************************************************************************************
!> \brief  Plot a GW
!!
!! \param fname1  Input file name
!! \param thingy  Thingy (something with the axes)
!! \param lbl     Plot label

subroutine plotwave(fname1,thingy,lbl)
  use SUFR_kinds, only: double, dbl
  
  implicit none
  character, intent(in) :: fname1*(99), lbl*(99)
  integer, intent(in) :: thingy
  
  integer, parameter :: nf=1, n1=1000000
  integer :: i,f,n(nf),io,lw
  integer :: nfrx,nfry
  real :: xwinmin,xwinmax,ywinmin,ywinmax
  real :: t(nf,n1),h(nf,n1),dx,dy,xmin,xmax,ymin,ymax
  real(double) :: t1,t0,m1,m2,mc,eta,tc,dl,lat,lon,phase,spin,kappa,thJ0,phJ0,alpha
  character :: fname*(299),bla
  
  t0 = 0.0_dbl
  
  do f=1,nf
     fname = trim(fname1)
     !if(f.eq.2) fname = fname2
     open(unit=10,form='formatted',status='old',file=trim(fname),iostat=io)
     rewind(10)
     
     !write(6,'(A)', advance='no')'Reading input file '//trim(fname)//'...     '
     read(10,*)bla
     read(10,*)m1,m2,mc,eta,tc,dl,lat,lon,phase,spin,kappa,thJ0,phJ0,alpha
     read(10,*)bla
     do i=1,n1
        read(10,*,err=195,end=199)t1,h(f,i)
        if(i.eq.1) t0 = (nint(t1*1.d-3)-1)*1.d3
        t(f,i) = real(t1-t0)
     end do
195  write(6,'(A,I8)')'error reading file '//trim(fname)//' line',i+1
199  close(10)
     n(f) = i-1
     !write(6,'(I,A12)')n(f),'lines read.'
  end do
  
  h = h*1.e21
  
  xmin =  1.e30
  xmax = -1.e30
  ymin =  1.e30
  ymax = -1.e30
  do f=1,nf
     xmin = min(minval(t(f,1:n(f))),xmin)
     xmax = max(maxval(t(f,1:n(f))),xmax)
     ymin = min(minval(h(f,1:n(f))),ymin)
     ymax = max(maxval(h(f,1:n(f))),ymax)
  end do
  
  ymax = maxval(abs(h))
  ymin = -ymax
  
  xmin = 1006.65
  xmin = 1010.55
  xmin = 1011.6
  xmax = 1012.32
  !xmin = 1012.3
  xmax = 1012.36
  
  
  
  dx = abs(xmax-xmin)
  dy = abs(ymax-ymin)
  
  lw = 3
  call pgslw(lw)
  nfrx = 1      !Number of frames
  nfry = nf      
  xwinmin = 0.10    !Absolute limits
  xwinmax = 0.95
  ywinmin = 0.20
  ywinmax = 0.95
  if(nf.eq.1) ywinmin = 0.25    
  if(nf.eq.1) ywinmax = 0.8
  if(nf.eq.2) ywinmax = 0.9
  if(nf.eq.3) ywinmax = 0.95
  !call pgswin(xmin,xmax,2*ymin-dy*0.05,2*ymax+dy*0.05)
  call pgswin(xmin,xmax,ymin-dy*0.05,ymax+dy*0.05)
  if(thingy.eq.1) call pgbox('BCTS',0.0,0,'BCNTS',0.0,0)
  if(thingy.eq.2) call pgbox('BCNTS',0.0,0,'BCNTS',0.0,0)
  !call pgbox('BCNTS',0.0,0,'BC',0.0,0)
  if(thingy.eq.2) call pgmtxt('B',2.4,0.5,0.5,'time (s)')
  call pgmtxt('L',2.2,0.5,0.5,'h (10\u-21\d)')
  do f=1,nf
     call pgslw(1)
     !call pgline(n(f),t(f,1:n(f)),h(f,1:n(f))-(f-1.5)*dy)
     call pgline(n(f), t(f,1:n(f)), h(f,1:n(f)) )
     !call pgsch(fontsize)
     call pgsci(1)
     call pgslw(lw)
  end do  !f
  call pgsch(1.)
  call pgmtxt('T',-1.5,0.03,0.,trim(lbl))
  
  
end subroutine plotwave
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Read the compPDF settings file

subroutine read_inputfile()
  use comp_pdfs_settings, only: nf,fnames,fontsize,clr,fillstyle,plrange,pltrue,plmedian, plpars,frames,dirnames
  use comp_pdfs_settings, only: settingsfile,file,type,dim,fonttype,plpars2d,outnamebase
  
  implicit none
  integer :: u,io
  character :: bla
  
  u = 14
  open(unit=u,form='formatted',status='old',action='read',file=trim(settingsfile),iostat=io)
  if(io.ne.0) write(*,'(A)')'  Error opening input file.'
  
  
  read(u,*)bla
  
  read(u,*)bla
  read(u,*)nf
  read(u,*)file
  read(u,*)type
  read(u,*)dim
  read(u,*)pltrue
  read(u,*)plmedian
  read(u,*)plrange
  read(u,*)clr
  read(u,*)fillstyle
  read(u,*)fonttype
  read(u,*)fontsize
  
  read(u,*)bla
  read(u,*)bla
  read(u,*)frames
  read(u,*)bla
  plpars = 0
  !read(u,*)plpars
  read(u,*)plpars(1:frames(1)*frames(2))  !Don't always read 20 parameters
  
  read(u,*)bla
  read(u,*)bla
  read(u,*)plpars2d
  
  read(u,*)bla
  read(u,*)bla
  read(u,*)fnames(1:nf)
  read(u,*)bla
  read(u,*)dirnames(1:nf)
  read(u,*)bla
  read(u,*)outnamebase
  
  close(u)
end subroutine read_inputfile
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Write the compPDF settings file

subroutine write_inputfile()
  use comp_pdfs_settings, only: nf,fnames,fontsize,clr,fillstyle,plrange,pltrue,plmedian, plpars,frames,dirnames
  use comp_pdfs_settings, only: file,type,dim,fonttype,plpars2d,outnamebase
  
  implicit none
  integer :: u,i,io
  
  u = 14
  open(unit=u,form='formatted',status='replace',file='comp_pdfs.used',iostat=io)
  if(io.ne.0) write(*,'(A)')'  Error writing input file.'
  
11 format(I10,1x,A9,5x,A)
12 format(1x,2I5)
13 format(3x,20I3)
  
21 format(F10.5,1x,A9,5x,A)
  
  
  write(u,'(A,/)')' Input file for comp_pdfs.f'
  
  
  write(u,'(/,A)')' General plot options:'
  write(u,11)nf, 'nf',   'Number of input files'
  write(u,11)file, 'file',   'Output file type: 0:screen, 1: png, 2:eps, 3:pdf'
  write(u,11)type, 'type',   'Type of plot:  1: default, 2: poster, 3: talk'
  write(u,11)dim, 'dim',   'Dimension of PDF: 1 or 2'
  write(u,11)pltrue, 'pltrue',   &
       'Plot true values for 1D or 2D plots:  0-no, 1-yes, 2-use different line styles for multiple sets of true values'
  write(u,11)plmedian, 'plmedian',   'Plot medians for 1D or 2D plots, only if reading 1 file'
  write(u,11)plrange, 'plrange',  'Plot ranges  for 1D or 2D plots, only if reading 1 file: 0-no, 1:plot lines, 2:plot numbers too'
  write(u,11)clr, 'clr',   'Use colour: 0-B/W, 1-colour, 2-grey scales'
  write(u,11)fillstyle, 'fillstyle',   'Fill style for PDFs: 1-solid, 2-outline, 3-hatched, 4-cross-hatched'
  write(u,11)fonttype, 'fonttype',   'Font type: 1: arial, 2: woman... no no, rrrroman, 3: italic, 4: script'
  write(u,21)fontsize, 'fontsize',   'Scale of font, etc.  Default: 2.'
  
  
  write(u,'(//,A)')' 1D plot options:'
  write(u,'(A)')'   Number of frames in 1D plot: (w, h  or  x, y)  (frames):'
  write(u,12)frames
  write(u,'(A)') &
       "   Labels of the plot parameters for the 1d frames (plpars)  (need at least one for each frame.  0-don't plot panel):  "
  write(u,13)plpars !Array of size 20
  
  
  write(u,'(//,A)')' 2D plot options:'
  write(u,'(A)')'   Labels of the plot parameters for the 2d plot (plpars2d):'
  write(u,12)plpars2d
  
  
  write(u,'(//,A)')' Input files:'
  write(u,'(A)')'   Base of input files names (without extension) (fnames):'
  write(u,'(A)', advance='no')'     '
  do i=1,nf
     write(u,'(A)', advance='no')trim(fnames(i))//'   '
  end do
  write(u,*)
  
  write(u,'(/,A)')'   Directories of the input files (dirnames):'
  write(u,'(A)', advance='no')'     '
  do i=1,nf
     write(u,'(A)', advance='no')trim(dirnames(i))//'   '
  end do
  write(u,*)
  
  write(u,'(/,A)')'   Base of the output file name (without extention, can start with dir) (outnamebase)'
  write(u,'(A)', advance='no')'     '//trim(outnamebase)
  
  
  write(u,*)
  close(u)
end subroutine write_inputfile
!***********************************************************************************************************************************



!***********************************************************************************************************************************
!> \brief  Set the title in a Postscript file generated by PGPlot
!! 
!! \param PSfile   The file to adapt
!! \param PStitle  The title to give to the PS file
!!
!! \note  Requires sed installed

subroutine set_PGPS_title(PSfile,PStitle)
  implicit none
  character, intent(in) :: PSfile*(*),PStitle*(*)
  integer :: status,system
  character :: tempfile*(99)
  
  tempfile = 'temp_PGPS_file.eps'
  
  status = system("sed -e 's/Title: PGPLOT PostScript plot/Title: "//trim(PStitle)//"/' "//trim(PSfile)//" > "//trim(tempfile))
  if(status.eq.0) then
     status = system('mv -f '//trim(tempfile)//' '//trim(PSfile))
  else
     status = system('rm -f '//trim(tempfile))
  end if
  
end subroutine set_PGPS_title
!***********************************************************************************************************************************


